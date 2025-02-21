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
o#include <external/resinsight/ReservoirDataModel/cvfGeometryTools.h>
#include <external/resinsight/ReservoirDataModel/RigWellLogExtractor.h>
#include <external/resinsight/ReservoirDataModel/RigCellGeometryTools.h>
#include <external/resinsight/CommonCode/cvfStructGrid.h>
#include <external/resinsight/LibGeometry/cvfBoundingBox.h>

#include <map>


namespace external {

 RigEclipseWellLogExtractorGrid::RigEclipseWellLogExtractorGrid( const RigWellPath* wellpath, 
                                                         const Opm::CpGrid& grid, 
                                                         cvf::ref<cvf::BoundingBoxTree>& cellSearchTree )
    : RigWellLogExtractor( wellpath, "" )
      ,m_grid(grid)
      ,m_cellSearchTree(cellSearchTree)
{
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
            RigEclipseWellLogExtractor::hexCornersOpmToResinsight( hexCorners, globalCellIndex);
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
cvf::Vec3d RigEclipseWellLogExtractor::calculateLengthInCell( size_t            cellIndex,
                                                              const cvf::Vec3d& startPoint,
                                                              const cvf::Vec3d& endPoint ) const
{
    cvf::Vec3d hexCorners[8];  //resinsight numbering, see RigCellGeometryTools.cpp
    RigEclipseWellLogExtractor::hexCornersOpmToResinsight( hexCorners, cellIndex );

    std::array<cvf::Vec3d, 8> hexCorners2;
    for (size_t i = 0; i < 8; i++)
           hexCorners2[i] =  hexCorners[i];

    return RigEclipseWellLogExtractor::calculateLengthInCell( hexCorners2, startPoint, endPoint );
}

// Modified version of ApplicationLibCode\ReservoirDataModel\RigWellPathIntersectionTools.cpp
cvf::Vec3d RigEclipseWellLogExtractor::calculateLengthInCell( const std::array<cvf::Vec3d, 8>& hexCorners,
                                                              const cvf::Vec3d&                startPoint,
                                                              const cvf::Vec3d&                endPoint ) const
{

    cvf::Vec3d vec = endPoint - startPoint;
    cvf::Vec3d iAxisDirection;
    cvf::Vec3d jAxisDirection;
    cvf::Vec3d kAxisDirection;

    RigEclipseWellLogExtractor::findCellLocalXYZ( hexCorners, iAxisDirection, jAxisDirection, kAxisDirection );

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
void RigEclipseWellLogExtractor::findCellLocalXYZ( const std::array<cvf::Vec3d, 8>& hexCorners,
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
void RigEclipseWellLogExtractor::hexCornersOpmToResinsight( cvf::Vec3d hexCorners[8], size_t cellIndex ) const
{
    const auto[i,j,k] = m_grid.getIJK(cellIndex);
    std::array<std::size_t, 8> opm2resinsight = {0, 1, 3, 2, 4, 5, 7, 6};
   
    for (std::size_t l = 0; l < 8; l++) {
         std::array<double, 3> cornerPointArray;
         cornerPointArray = m_grid.getCornerPos(i,j,k,l);
         hexCorners[opm2resinsight[l]]= cvf::Vec3d(cornerPointArray[0], cornerPointArray[1], cornerPointArray[2]);
    } 
}

// Modified version of ApplicationLibCode\ReservoirDataModel\RigMainGrid.cpp
void RigEclipseWellLogExtractor::buildCellSearchTree()
{
    if (m_cellSearchTree.isNull()) {

        auto nx = m_grid.getNX();
        auto ny = m_grid.getNY();
        auto nz = m_grid.getNZ();

        size_t cellCount = nx * ny * nz;
        std::vector<size_t>           cellIndicesForBoundingBoxes;
        std::vector<cvf::BoundingBox> cellBoundingBoxes;

// #pragma omp parallel
//         {
            size_t threadCellCount = cellCount;

            std::vector<size_t>           threadIndicesForBoundingBoxes;
            std::vector<cvf::BoundingBox> threadBoundingBoxes;

            threadIndicesForBoundingBoxes.reserve( threadCellCount );
            threadBoundingBoxes.reserve( threadCellCount );

            std::array<double, 3> cornerPointArray;
            cvf::Vec3d cornerPoint;

// #pragma omp for
            for (int cIdx = 0; cIdx < (int)cellCount; ++cIdx) {
                const auto[i,j,k] = m_grid.getIJK(cIdx);
                cvf::BoundingBox cellBB;
                for (std::size_t l = 0; l < 8; l++) {
                     cornerPointArray = m_grid.getCornerPos(i,j,k,l);
                     cornerPoint = cvf::Vec3d(cornerPointArray[0], cornerPointArray[1], cornerPointArray[2]);
                     cellBB.add(cornerPoint);
                }

                if (cellBB.isValid()) {
                    threadIndicesForBoundingBoxes.emplace_back(cIdx);
                    threadBoundingBoxes.emplace_back(cellBB);
                }
            }

            threadIndicesForBoundingBoxes.shrink_to_fit();
            threadBoundingBoxes.shrink_to_fit();
// #pragma omp critical
//          {
             cellIndicesForBoundingBoxes.insert( cellIndicesForBoundingBoxes.end(),
                                                 threadIndicesForBoundingBoxes.begin(),
                                                 threadIndicesForBoundingBoxes.end() );

             cellBoundingBoxes.insert( cellBoundingBoxes.end(), threadBoundingBoxes.begin(), threadBoundingBoxes.end() );
        //  }
    // } #pragma omp parallel
         m_cellSearchTree = new cvf::BoundingBoxTree;
         m_cellSearchTree->buildTreeFromBoundingBoxes( cellBoundingBoxes, &cellIndicesForBoundingBoxes );
    }
}

// From ApplicationLibCode\ReservoirDataModel\RigMainGrid.cpp
void RigEclipseWellLogExtractor::findIntersectingCells( const cvf::BoundingBox& inputBB, std::vector<size_t>* cellIndices ) const
{
    CVF_ASSERT( m_cellSearchTree.notNull() );

    m_cellSearchTree->findIntersections( inputBB, cellIndices );
}

// Modified version of ApplicationLibCode\ReservoirDataModel\RigEclipseWellLogExtractor.cpp 
std::vector<size_t> RigEclipseWellLogExtractor::findCloseCellIndices( const cvf::BoundingBox& bb )
{
    std::vector<size_t> closeCells;
    this->findIntersectingCells( bb, &closeCells );
    return closeCells;
}

 cvf::ref<cvf::BoundingBoxTree> RigEclipseWellLogExtractor::getCellSearchTree()
 {
    return m_cellSearchTree;
 }
 } //namespace external
