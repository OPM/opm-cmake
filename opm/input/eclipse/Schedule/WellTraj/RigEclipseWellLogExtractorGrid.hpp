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

#pragma once

#include <external/resinsight/ReservoirDataModel/RigWellLogExtractor.h>
#include <external/resinsight/LibGeometry/cvfBoundingBoxTree.h>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/Schedule/ScheduleGrid.hpp>

// from origina cp file
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



namespace Opm {
class EclipseGrid;
class ScheduleGrid;
}
namespace Dune{
  class CpGrid;
}

namespace external {

class RigEclipseCaseData;
class RigWellPath;
class RigResultAccessor;

namespace cvf
{
class BoundingBox;
}

//==================================================================================================
///
//==================================================================================================
class RigEclipseWellLogExtractorGrid : public RigWellLogExtractor
{
public:
  RigEclipseWellLogExtractorGrid( const RigWellPath* wellpath, const Dune::CpGrid& grid, cvf::ref<cvf::BoundingBoxTree>& cellSearchTree);

    cvf::ref<cvf::BoundingBoxTree> getCellSearchTree();
private:
    void                calculateIntersection();
    std::vector<size_t> findCloseCellIndices( const cvf::BoundingBox& bb );
    cvf::Vec3d
        calculateLengthInCell( size_t cellIndex, const cvf::Vec3d& startPoint, const cvf::Vec3d& endPoint ) const override;

    cvf::Vec3d calculateLengthInCell( const std::array<cvf::Vec3d, 8>& hexCorners,
                                      const cvf::Vec3d&                startPoint,
                                      const cvf::Vec3d&                endPoint ) const;

    void hexCornersOpmToResinsight( cvf::Vec3d    hexCorners[8],
                                    size_t        cellIndex ) const;

    void findCellLocalXYZ( const std::array<cvf::Vec3d, 8>& hexCorners,
                           cvf::Vec3d&                      localXdirection,
                           cvf::Vec3d&                      localYdirection,
                           cvf::Vec3d&                      localZdirection ) const;
    void buildCellSearchTree();
    void findIntersectingCells( const cvf::BoundingBox& inputBB, std::vector<size_t>* cellIndices ) const;
    void computeCachedData();

    const Dune::CpGrid& m_grid;
    cvf::ref<cvf::BoundingBoxTree> m_cellSearchTree;
    std::vector<std::array<cvf::Vec3d,8>> all_hexCorners_;
};
} //namespace external
//#include "RigEclipseWellLogExtractorGrid_impl.hpp"
