/////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (C) Statoil ASA
//  Copyright (C) Ceetron Solutions AS
//
//  ResInsight is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  ResInsight is distributed in the hope that it will be useful, but WITHOUT ANY
//  WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.
//
//  See the GNU General Public License at <http://www.gnu.org/licenses/gpl.html>
//  for more details.
//
/////////////////////////////////////////////////////////////////////////////////

// This file is based on the following files of ResInsight:
//     ApplicationLibCode\ReservoirDataModel\RigEclipseWellLogExtractor.cpp
//     ApplicationLibCode\ReservoirDataModel\RigCellGeometryTools.cpp
//     ApplicationLibCode\ReservoirDataModel\RigMainGrid.cpp
//     ApplicationLibCode\ReservoirDataModel\RigWellPathIntersectionTools.cpp


#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/Schedule/ScheduleGrid.hpp>

#include <opm/input/eclipse/Schedule/WellTraj/RigEclipseWellLogExtractor.hpp>
#include <external/resinsight/ReservoirDataModel/RigWellLogExtractionTools.h>
#include <external/resinsight/ReservoirDataModel/RigWellPath.h>
#include <external/resinsight/ReservoirDataModel/cvfGeometryTools.h>
#include <external/resinsight/ReservoirDataModel/RigWellLogExtractor.h>
#include <external/resinsight/ReservoirDataModel/RigCellGeometryTools.h>
#include <external/resinsight/CommonCode/cvfStructGrid.h>
#include <external/resinsight/LibGeometry/cvfBoundingBox.h>

#include <map>

#include <opm/grid/common/CartesianIndexMapper.hpp>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/gridenums.hh>
#include <opm/grid/CpGrid.hpp>


namespace external {

 RigEclipseWellLogExtractorGrid::RigEclipseWellLogExtractorGrid( const RigWellPath* wellpath, 
                                                         const Dune::CpGrid& grid, 
                                                         cvf::ref<cvf::BoundingBoxTree>& cellSearchTree )
    : RigWellLogExtractor( wellpath, "" )
      ,m_grid(grid)
      ,m_cellSearchTree(cellSearchTree)
{
  // could be done when making bounding box
        using Grid = ::Dune::CpGrid;
        //using GridView = typename Grid::LeafGridView;
    using GridView = typename Grid::LeafGridView;
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    const auto& gv = m_grid.leafGridView();
    ElementMapper mapper(gv, Dune::mcmgElementLayout()); // used id sets interall
    for(const auto& elem : Dune::elements(gv)){
      auto geom = elem.geometry();
      int index = mapper.index(elem);
      assert(geom.corners() == 8);
      cvf::BoundingBox cellBB;
      //cvf::Vec3d cornerPoint;
      std::array<cvf::Vec3d,8> cornerPointArray;
      // dune 0 1 2 3 4 5 6 7 is resinsight 0 1 3 2 4 5 7 6 (i think)
      for (std::size_t l = 0; l < 8; l++) {
	auto duneCornerPointArray = geom.corner(l);
	cvf::Vec3d cornerPoint = {duneCornerPointArray[0], duneCornerPointArray[1], duneCornerPointArray[2]};     
	cornerPointArray[l] = cornerPoint;
      }
      all_hexCorners_.push_back(cornerPointArray);
    }
    calculateIntersection();
}


// Modified version of ApplicationLibCode\ReservoirDataModel\RigEclipseWellLogExtractor.cpp
void RigEclipseWellLogExtractorGrid::calculateIntersection()
{
    std::map<RigMDCellIdxEnterLeaveKey, HexIntersectionInfo> uniqueIntersections;

    if ( m_wellPathGeometry->wellPathPoints().empty() ) return;

    buildCellSearchTree();
    

    for ( size_t wpp = 0; wpp < m_wellPathGeometry->wellPathPoints().size() - 1; ++wpp )
    {
        std::vector<HexIntersectionInfo> intersections;
        cvf::Vec3d                       p1 = m_wellPathGeometry->wellPathPoints()[wpp];
        cvf::Vec3d                       p2 = m_wellPathGeometry->wellPathPoints()[wpp + 1];

        cvf::BoundingBox bb;

        bb.add( p1 );
        bb.add( p2 );

        std::vector<size_t> closeCellIndices = findCloseCellIndices( bb );

        for ( const auto& globalCellIndex : closeCellIndices )
        {
            cvf::Vec3d    hexCorners[8];  //resinsight numbering, see RigCellGeometryTools.cpp
            RigEclipseWellLogExtractorGrid::hexCornersOpmToResinsight( hexCorners, globalCellIndex);
	    // for (std::size_t l = 0; l < 8; l++) {
	    //   hexCorners[l] = all_hexCorners_[globalCellIndex][l];
	    // }
            RigHexIntersectionTools::lineHexCellIntersection( p1, p2, hexCorners, globalCellIndex, &intersections );
        }

        // Now, with all the intersections of this piece of line, we need to
        // sort them in order, and set the measured depth and corresponding cell index

        // Inserting the intersections in this map will remove identical intersections
        // and sort them according to MD, CellIdx, Leave/enter

        double md1 = m_wellPathGeometry->measuredDepths()[wpp];
        double md2 = m_wellPathGeometry->measuredDepths()[wpp + 1];

        insertIntersectionsInMap( intersections, p1, md1, p2, md2, &uniqueIntersections );
    }

    this->populateReturnArrays( uniqueIntersections );
}

// Modified version of ApplicationLibCode\ReservoirDataModel\RigEclipseWellLogExtractor.cpp
cvf::Vec3d RigEclipseWellLogExtractorGrid::calculateLengthInCell( size_t            cellIndex,
                                                              const cvf::Vec3d& startPoint,
                                                              const cvf::Vec3d& endPoint ) const
{
    cvf::Vec3d hexCorners[8];  //resinsight numbering, see RigCellGeometryTools.cpp
    RigEclipseWellLogExtractorGrid::hexCornersOpmToResinsight( hexCorners, cellIndex );

    std::array<cvf::Vec3d, 8> hexCorners2;
    for (size_t i = 0; i < 8; i++)
           hexCorners2[i] =  hexCorners[i];

    return RigEclipseWellLogExtractorGrid::calculateLengthInCell( hexCorners2, startPoint, endPoint );
}

// Modified version of ApplicationLibCode\ReservoirDataModel\RigWellPathIntersectionTools.cpp
cvf::Vec3d RigEclipseWellLogExtractorGrid::calculateLengthInCell( const std::array<cvf::Vec3d, 8>& hexCorners,
                                                              const cvf::Vec3d&                startPoint,
                                                              const cvf::Vec3d&                endPoint ) const
{

    cvf::Vec3d vec = endPoint - startPoint;
    cvf::Vec3d iAxisDirection;
    cvf::Vec3d jAxisDirection;
    cvf::Vec3d kAxisDirection;

    RigEclipseWellLogExtractorGrid::findCellLocalXYZ( hexCorners, iAxisDirection, jAxisDirection, kAxisDirection );

    cvf::Mat3d localCellCoordinateSystem( iAxisDirection.x(),
                                          jAxisDirection.x(),
                                          kAxisDirection.x(),
                                          iAxisDirection.y(),
                                          jAxisDirection.y(),
                                          kAxisDirection.y(),
                                          iAxisDirection.z(),
                                          jAxisDirection.z(),
                                          kAxisDirection.z() );

    auto signedVector = vec.getTransformedVector( localCellCoordinateSystem.getInverted() );

    return { std::fabs( signedVector.x() ), std::fabs( signedVector.y() ), std::fabs( signedVector.z() ) };
}

// Modified version of ApplicationLibCode\ReservoirDataModel\RigCellGeometryTools.cpp
void RigEclipseWellLogExtractorGrid::findCellLocalXYZ( const std::array<cvf::Vec3d, 8>& hexCorners,
                                                   cvf::Vec3d&                      localXdirection,
                                                   cvf::Vec3d&                      localYdirection,
                                                   cvf::Vec3d&                      localZdirection ) const
{

    cvf::Vec3d faceCenterNegI = cvf::GeometryTools::computeFaceCenter( hexCorners[0],
                                                                       hexCorners[4],
                                                                       hexCorners[7],
                                                                       hexCorners[3] );

    cvf::Vec3d faceCenterPosI = cvf::GeometryTools::computeFaceCenter( hexCorners[1],
                                                                       hexCorners[2],
                                                                       hexCorners[6],
                                                                       hexCorners[5] );

    cvf::Vec3d faceCenterNegJ = cvf::GeometryTools::computeFaceCenter( hexCorners[0],
                                                                       hexCorners[1],
                                                                       hexCorners[5],
                                                                       hexCorners[4] );
    
    cvf::Vec3d faceCenterPosJ = cvf::GeometryTools::computeFaceCenter( hexCorners[3],
                                                                       hexCorners[7],
                                                                       hexCorners[6],
                                                                       hexCorners[2] );

    cvf::Vec3d faceCenterCenterVectorI = faceCenterPosI - faceCenterNegI;
    cvf::Vec3d faceCenterCenterVectorJ = faceCenterPosJ - faceCenterNegJ;

    localZdirection.cross( faceCenterCenterVectorI, faceCenterCenterVectorJ );
    localZdirection.normalize();

    cvf::Vec3d crossPoductJZ;
    crossPoductJZ.cross( faceCenterCenterVectorJ, localZdirection );
    localXdirection = faceCenterCenterVectorI + crossPoductJZ;
    localXdirection.normalize();

    cvf::Vec3d crossPoductIZ;
    crossPoductIZ.cross( faceCenterCenterVectorI, localZdirection );
    localYdirection = faceCenterCenterVectorJ - crossPoductIZ;
    localYdirection.normalize();
}

// Convert opm to resinsight numbering of cornerpoints, see RigCellGeometryTools.cpp
void RigEclipseWellLogExtractorGrid::hexCornersOpmToResinsight( cvf::Vec3d hexCorners[8], size_t cellIndex ) const
{
  //NB NB this should depend on grid type
    std::array<std::size_t, 8> opm2resinsight = {0, 1, 3, 2, 4, 5, 7, 6};
    for (std::size_t l = 0; l < 8; l++) {
      //std::array<double, 3> cornerPointArray;
         cvf::Vec3d cornerPointArray = all_hexCorners_[cellIndex][l];
         hexCorners[opm2resinsight[l]]= cvf::Vec3d(cornerPointArray[0], cornerPointArray[1], cornerPointArray[2]);
    } 
}

// Modified version of ApplicationLibCode\ReservoirDataModel\RigMainGrid.cpp
void RigEclipseWellLogExtractorGrid::buildCellSearchTree()
{
    
    if (m_cellSearchTree.isNull()) {
      assert(false);
    }
}

// From ApplicationLibCode\ReservoirDataModel\RigMainGrid.cpp
void RigEclipseWellLogExtractorGrid::findIntersectingCells( const cvf::BoundingBox& inputBB, std::vector<size_t>* cellIndices ) const
{
    CVF_ASSERT( m_cellSearchTree.notNull() );

    m_cellSearchTree->findIntersections( inputBB, cellIndices );
}

// Modified version of ApplicationLibCode\ReservoirDataModel\RigEclipseWellLogExtractorGrid.cpp 
std::vector<size_t> RigEclipseWellLogExtractorGrid::findCloseCellIndices( const cvf::BoundingBox& bb )
{
    std::vector<size_t> closeCells;
    this->findIntersectingCells( bb, &closeCells );
    return closeCells;
}

 cvf::ref<cvf::BoundingBoxTree> RigEclipseWellLogExtractorGrid::getCellSearchTree()
 {
    return m_cellSearchTree;
 }
 } //namespace external
