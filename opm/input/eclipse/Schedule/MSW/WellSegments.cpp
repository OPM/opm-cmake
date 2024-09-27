/*
  Copyright 2016 Statoil ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <opm/input/eclipse/Schedule/MSW/WellSegments.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/input/eclipse/Schedule/MSW/Segment.hpp>
#include <opm/input/eclipse/Schedule/MSW/SICD.hpp>
#include <opm/input/eclipse/Schedule/MSW/Valve.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>

#include <opm/input/eclipse/Deck/DeckItem.hpp>
#include <opm/input/eclipse/Deck/DeckKeyword.hpp>
#include <opm/input/eclipse/Deck/DeckRecord.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#ifdef _WIN32
#define _USE_MATH_DEFINES
#include <math.h>
#endif

#include <fmt/format.h>

namespace Opm {

    WellSegments::WellSegments(CompPressureDrop compDrop,
                               const std::vector<Segment>& segments)
       : m_comp_pressure_drop(compDrop)
    {
        for (const auto& segment : segments)
            this->addSegment(segment);
    }


    WellSegments::WellSegments(const DeckKeyword& keyword) {
        this->loadWELSEGS(keyword);
    }


    WellSegments WellSegments::serializationTestObject()
    {
        WellSegments result;
        result.m_comp_pressure_drop = CompPressureDrop::HF_;
        result.m_segments = {Opm::Segment::serializationTestObject()};
        result.segment_number_to_index = {{1, 2}};

        return result;
    }


    std::size_t WellSegments::size() const {
        return m_segments.size();
    }

    bool WellSegments::empty() const {
        return this->m_segments.empty();
    }

    int WellSegments::maxSegmentID() const
    {
        return std::accumulate(this->m_segments.begin(),
                               this->m_segments.end(), 0,
            [](const int maxID, const Segment& seg)
        {
            return std::max(maxID, seg.segmentNumber());
        });
    }

    int WellSegments::maxBranchID() const
    {
        return std::accumulate(this->m_segments.begin(),
                               this->m_segments.end(), 0,
            [](const int maxID, const Segment& seg)
        {
            return std::max(maxID, seg.branchNumber());
        });
    }

    const Segment& WellSegments::topSegment() const {
        return this->m_segments[0];
    }

    double WellSegments::depthTopSegment() const {
        return this->topSegment().depth();
    }

    double WellSegments::lengthTopSegment() const {
        return this->topSegment().totalLength();
    }

    double WellSegments::volumeTopSegment() const {
        return this->topSegment().volume();
    }


    WellSegments::CompPressureDrop WellSegments::compPressureDrop() const {
        return m_comp_pressure_drop;
    }

    const std::vector<Segment>::const_iterator WellSegments::begin() const {
        return this->m_segments.begin();
    }

    const std::vector<Segment>::const_iterator WellSegments::end() const {
        return this->m_segments.end();
    }

    const Segment& WellSegments::operator[](size_t idx) const {
        return m_segments[idx];
    }

    int WellSegments::segmentNumberToIndex(const int segment_number) const {
        const auto it = segment_number_to_index.find(segment_number);
        if (it != segment_number_to_index.end()) {
            return it->second;
        } else {
            return -1;
        }
    }

    void WellSegments::addSegment(const Segment& new_segment)
    {
        // Decide whether to push_back or insert
        const auto segment_number = new_segment.segmentNumber();
        const auto segment_index = segmentNumberToIndex(segment_number);
        if (segment_index < 0) {
            // New segment object.
            const auto new_index = static_cast<int>(this->size());
            this->segment_number_to_index.insert_or_assign(segment_number, new_index);
            this->m_segments.push_back(new_segment);
        }
        else {
            // Update to existing segment object.
            this->m_segments[segment_index] = new_segment;
        }
    }

    void WellSegments::addSegment(const int segment_number,
                                  const int branch,
                                  const int outlet_segment,
                                  const double length,
                                  const double depth,
                                  const double internal_diameter,
                                  const double roughness,
                                  const double cross_area,
                                  const double volume,
                                  const bool data_ready,
                                  const double node_x,
                                  const double node_y)
    {
        const auto segment = Segment {
            segment_number, branch, outlet_segment,
            length, depth, internal_diameter, roughness,
            cross_area, volume,
            data_ready, node_x, node_y
        };

        this->addSegment(segment);
    }


    void WellSegments::loadWELSEGS(const DeckKeyword& welsegsKeyword)
    {
        // For the first record, which provides the information for the top
        // segment and information for the whole segment set.
        const auto& record1 = welsegsKeyword.getRecord(0);

        // Meaningless value to indicate unspecified values.
        const double invalid_value = Segment::invalidValue();

        const double depth_top = record1.getItem("TOP_DEPTH").getSIDouble(0);
        const double length_top = record1.getItem("TOP_LENGTH").getSIDouble(0);
        const double volume_top = record1.getItem("WELLBORE_VOLUME").getSIDouble(0);
        const LengthDepth length_depth_type = LengthDepthFromString(record1.getItem("INFO_TYPE").getTrimmedString(0));
        m_comp_pressure_drop = CompPressureDropFromString(record1.getItem("PRESSURE_COMPONENTS").getTrimmedString(0));

        const auto nodeX_top = record1.getItem("TOP_X").getSIDouble(0);
        const auto nodeY_top = record1.getItem("TOP_Y").getSIDouble(0);

        // The main branch is 1 instead of 0.  The segment number for top
        // segment is also 1.
        if (length_depth_type == LengthDepth::INC) {
            const auto segmentID = 1;
            const auto branchID = 1;
            const auto outletSegment = 0;
            const auto length = 0.0;
            const auto depth = 0.0;
            const auto internal_diameter = invalid_value;
            const auto roughness = invalid_value;
            const auto cross_area = invalid_value;
            const auto data_ready = false;

            this->addSegment(segmentID, branchID, outletSegment, length, depth,
                             internal_diameter, roughness, cross_area,
                             volume_top, data_ready, nodeX_top, nodeY_top);
        }
        else if (length_depth_type == LengthDepth::ABS) {
            const auto segmentID = 1;
            const auto branchID = 1;
            const auto outletSegment = 0;
            const auto internal_diameter = invalid_value;
            const auto roughness = invalid_value;
            const auto cross_area = invalid_value;
            const auto data_ready = true;

            this->addSegment(segmentID, branchID, outletSegment,
                             length_top, depth_top,
                             internal_diameter, roughness, cross_area,
                             volume_top, data_ready, nodeX_top, nodeY_top);
        }

        // Read all the information out from the DECK first then process to
        // get all the requisite information.
        for (size_t recordIndex = 1; recordIndex < welsegsKeyword.size(); ++recordIndex) {
            const auto& record = welsegsKeyword.getRecord(recordIndex);
            const int segment1 = record.getItem("SEGMENT1").get<int>(0);
            const int segment2 = record.getItem("SEGMENT2").get<int>(0);
            if (segment1 < 2) {
                throw std::logic_error {
                    fmt::format("Illegal segment 1 number in WELSEGS\n"
                                "Expected 2..NSEGMX, but got {}", segment1)
                };
            }

            if (segment2 < segment1) {
                throw std::logic_error {
                    fmt::format("Illegal segment 2 number in WELSEGS\n"
                                "Expected {}..NSEGMX, but got {}",
                                segment1, segment2)
                };
            }

            if ((segment1 != segment2) && (length_depth_type == LengthDepth::ABS)) {
                throw std::logic_error{
                    fmt::format("In WELSEGS, it is not supported to enter multiple segments in one record "
                                     "with ABS type of tubing length and depth information")
                };
            }

            const int branch = record.getItem("BRANCH").get<int>(0);
            if (branch < 1) {
                throw std::logic_error {
                    fmt::format("Illegal branch number input "
                                "({}) is found in WELSEGS!", branch)
                };
            }

            const double diameter = record.getItem("DIAMETER").getSIDouble(0);
            double area = M_PI * diameter * diameter / 4.0;
            {
                const auto& itemArea = record.getItem("AREA");
                if (itemArea.hasValue(0)) {
                    area = itemArea.getSIDouble(0);
                }
            }

            // If the length_depth_type is INC, then the length is the length of the segment,
            // If the length_depth_type is ABS, then the length is the length of the last segment node in the range.
            const double length = record.getItem("LENGTH").getSIDouble(0);

            // If the length_depth_type is INC, then the depth is the depth change of the segment from the outlet segment.
            // If the length_depth_type is ABS, then the depth is the absolute depth of last the segment node in the range.
            const double depth = record.getItem("DEPTH").getSIDouble(0);

            double volume = invalid_value;
            {
                const auto& itemVolume = record.getItem("VOLUME");
                if (itemVolume.hasValue(0)) {
                    volume = itemVolume.getSIDouble(0);
                }
                else if (length_depth_type == LengthDepth::INC) {
                    volume = area * length;
                }
            }

            const double roughness = record.getItem("ROUGHNESS").getSIDouble(0);

            const auto node_X = record.getItem("LENGTH_X").getSIDouble(0);
            const auto node_Y = record.getItem("LENGTH_Y").getSIDouble(0);

            for (int segment_number = segment1; segment_number <= segment2; ++segment_number) {
                // For the first or the only segment in the range is the one
                // specified in WELSEGS.  From the second segment in the
                // range, the outlet segment is the previous segment in the
                // range.
                const int outlet_segment = (segment_number == segment1)
                    ? record.getItem("JOIN_SEGMENT").get<int>(0)
                    : segment_number - 1;

                const auto data_ready = (length_depth_type != LengthDepth::INC)
                    && (segment_number == segment2);

                this->addSegment(segment_number, branch, outlet_segment,
                                 length, depth, diameter,
                                 roughness, area, volume, data_ready,
                                 node_X, node_Y);
            }
        }

        for (const auto& segment : this->m_segments) {
            const int outlet_segment = segment.outletSegment();
            if (outlet_segment <= 0) { // no outlet segment
                continue;
            }

            const int outlet_segment_index = segment_number_to_index[outlet_segment];
            m_segments[outlet_segment_index].addInletSegment(segment.segmentNumber());
        }

        this->process(length_depth_type, depth_top, length_top);
    }

    const Segment& WellSegments::getFromSegmentNumber(const int segment_number) const {
        // the index of segment in the vector of segments
        const int segment_index = segmentNumberToIndex(segment_number);
        if (segment_index < 0) {
            throw std::runtime_error {
                fmt::format("Unknown segment number {}", segment_number)
            };
        }
        return m_segments[segment_index];
    }

    void WellSegments::process(const LengthDepth length_depth,
                               const double      depth_top,
                               const double      length_top)
    {
        if (length_depth == LengthDepth::ABS) {
            this->processABS();
        }
        else if (length_depth == LengthDepth::INC) {
            this->processINC(depth_top, length_top);
        }
        else {
            throw std::logic_error {
                fmt::format("Invalid length/depth type "
                            "{} in segment data structure",
                            static_cast<int>(length_depth))
            };
        }
    }

    void WellSegments::processABS()
    {
        // Meaningless value to indicate unspecified/uncompleted values
        const double invalid_value = Segment::invalidValue();

        orderSegments();

        std::size_t current_index = 1;
        while (current_index < size()) {
            if (m_segments[current_index].dataReady()) {
                current_index++;
                continue;
            }

            const int range_begin = current_index;
            const int outlet_segment = m_segments[range_begin].outletSegment();
            const int outlet_index = segmentNumberToIndex(outlet_segment);

            assert(m_segments[outlet_index].dataReady() == true);

            std::size_t range_end = range_begin + 1;
            for (; range_end < size(); ++range_end) {
                if (m_segments[range_end].dataReady() == true) {
                    break;
                }
            }

            if (range_end >= size()) {
                throw std::logic_error(" One range records in WELSEGS is wrong. ");
            }

            // set the length and depth values in the range.
            int number_segments = range_end - range_begin + 1;
            assert(number_segments > 1); //if only 1, the information should be complete

            const double length_outlet = m_segments[outlet_index].totalLength();
            const double depth_outlet = m_segments[outlet_index].depth();

            const double length_last = m_segments[range_end].totalLength();
            const double depth_last = m_segments[range_end].depth();

            // incremental length and depth for the segments within the range
            const double length_inc = (length_last - length_outlet) / number_segments;
            const double depth_inc = (depth_last - depth_outlet) / number_segments;
            const double volume_segment = m_segments[range_end].crossArea() * length_inc;

            const auto x_outlet = m_segments[outlet_index].node_X();
            const auto y_outlet = m_segments[outlet_index].node_Y();

            const auto dx = (m_segments[range_end].node_X() - x_outlet) / number_segments;
            const auto dy = (m_segments[range_end].node_Y() - y_outlet) / number_segments;

            for (std::size_t k = range_begin; k <= range_end; ++k) {
                const auto& old_segment = this->m_segments[k];

                double new_length, new_depth, new_x, new_y;
                if (k == range_end) {
                    new_length = old_segment.totalLength();
                    new_depth  = old_segment.depth();
                    new_x      = old_segment.node_X();
                    new_y      = old_segment.node_Y();
                }
                else {
                    const auto num_inc = k - range_begin + 1;

                    new_length = length_outlet + num_inc*length_inc;
                    new_depth  = depth_outlet  + num_inc*depth_inc;
                    new_x      = x_outlet      + num_inc*dx;
                    new_y      = y_outlet      + num_inc*dy;
                }

                const auto new_volume = (old_segment.volume() < 0.5 * invalid_value)
                    ? volume_segment
                    : old_segment.volume();

                this->addSegment(Segment {
                    old_segment, new_length, new_depth,
                    new_volume, new_x, new_y
                });
            }

            current_index = range_end + 1;
        }

        // Then update the volume for all the segments except the top
        // segment.  This is for the segments specified individually while
        // the volume is not specified.
        for (std::size_t i = 1; i < size(); ++i) {
            assert(m_segments[i].dataReady());
            if (m_segments[i].volume() == invalid_value) {
                const auto& old_segment = this->m_segments[i];
                const int outlet_segment = m_segments[i].outletSegment();
                const int outlet_index = segmentNumberToIndex(outlet_segment);
                const double segment_length = m_segments[i].totalLength() - m_segments[outlet_index].totalLength();
                const double segment_volume = m_segments[i].crossArea() * segment_length;

                Segment new_segment(old_segment, segment_volume);
                addSegment(new_segment);
            }
        }
    }

    void WellSegments::processINC(double depth_top, double length_top)
    {
        // Update the information inside the WellSegments to be in ABS way
        Segment new_top_segment(this->m_segments[0], depth_top, length_top);
        this->addSegment(new_top_segment);

        orderSegments();

        // Begin with the second segment
        for (std::size_t i_index = 1; i_index < size(); ++i_index) {
            if (m_segments[i_index].dataReady()) {
                continue;
            }

            // Find its outlet segment
            const int outlet_segment = m_segments[i_index].outletSegment();
            const int outlet_index = segmentNumberToIndex(outlet_segment);

            // assert some information of the outlet_segment
            assert(outlet_index >= 0);
            assert(m_segments[outlet_index].dataReady());

            const double outlet_depth = m_segments[outlet_index].depth();
            const double outlet_length = m_segments[outlet_index].totalLength();
            const double temp_depth = outlet_depth + m_segments[i_index].depth();
            const double temp_length = outlet_length + m_segments[i_index].totalLength();
            const auto new_x = m_segments[outlet_index].node_X() + m_segments[i_index].node_X();
            const auto new_y = m_segments[outlet_index].node_Y() + m_segments[i_index].node_Y();

            // Applying the calculated length and depth to the current segment
            this->addSegment(Segment {
                this->m_segments[i_index],
                temp_depth, temp_length, new_x, new_y
            });
        }
    }

    void WellSegments::orderSegments()
    {
        // Re-ordering the segments to make later use easier.
        //
        // Two principles
        //
        //   1. A segment's outlet segment is ordered before the segment
        //      itself.
        //
        //   2. Segments on the same branch are stored consecutively.
        //
        // Top segment always at index zero so we only reorder the segments
        // in the index range [1 .. size()).
        std::size_t current_index = 1;

        // Clear the mapping from segment number to store index.
        segment_number_to_index.clear();

        // For the top segment
        segment_number_to_index.insert_or_assign(1, 0);

        while (current_index < size()) {
            // the branch number of the last segment that is done re-ordering
            const int last_branch_number = m_segments[current_index-1].branchNumber();
            // the one need to be swapped to the current_index.
            int target_segment_index= -1;

            // looking for target_segment_index
            for (std::size_t i_index= current_index; i_index< size(); ++i_index) {
                const int outlet_segment_number = m_segments[i_index].outletSegment();
                const int outlet_segment_index = segmentNumberToIndex(outlet_segment_number);
                if (outlet_segment_index < 0) { // not found the outlet_segment in the done re-ordering segments
                    continue;
                }
                if (target_segment_index< 0) { // first time found a candidate
                    target_segment_index= i_index;
                } else { // there is already a candidate, chosing the one with the same branch number with last_branch_number
                    const int old_target_segment_index_branch = m_segments[target_segment_index].branchNumber();
                    const int new_target_segment_index_branch = m_segments[i_index].branchNumber();
                    if (new_target_segment_index_branch == last_branch_number) {
                        if (old_target_segment_index_branch != last_branch_number) {
                            target_segment_index= i_index;
                        } else {
                            throw std::logic_error("two segments in the same branch share the same outlet segment !!\n");
                        }
                    }
                }
            }

            if (target_segment_index< 0) {
                throw std::logic_error("could not find candidate segment to swap in before the re-odering process get done !!\n");
            }

            assert(target_segment_index >= static_cast<int>(current_index));
            if (target_segment_index > static_cast<int>(current_index)) {
                std::swap(m_segments[current_index], m_segments[target_segment_index]);
            }
            const int segment_number = m_segments[current_index].segmentNumber();
            segment_number_to_index[segment_number] = current_index;
            current_index++;
        }
    }

    bool WellSegments::operator==( const WellSegments& rhs ) const {
        return this->m_comp_pressure_drop == rhs.m_comp_pressure_drop
            && this->m_segments.size() == rhs.m_segments.size()
            && this->segment_number_to_index.size() == rhs.segment_number_to_index.size()
            && std::equal( this->m_segments.begin(),
                           this->m_segments.end(),
                           rhs.m_segments.begin() )
            && std::equal( this->segment_number_to_index.begin(),
                           this->segment_number_to_index.end(),
                           rhs.segment_number_to_index.begin() );
    }

    double WellSegments::segmentLength(const int segment_number) const {
        const Segment& segment = this->getFromSegmentNumber(segment_number);
        if (segment_number == 1) // top segment
            return segment.totalLength();

        // other segments
        const Segment &outlet_segment = getFromSegmentNumber(segment.outletSegment());
        const double segment_length = segment.totalLength() - outlet_segment.totalLength();
        if (segment_length <= 0.)
            throw std::runtime_error(" non positive segment length is obtained for segment "
                                     + std::to_string(segment_number));

        return segment_length;
    }


    double WellSegments::segmentDepthChange(const int segment_number) const {
        const Segment& segment = getFromSegmentNumber(segment_number);
        if (segment_number == 1) // top segment
            return segment.depth();

        // other segments
        const Segment &outlet_segment = this->getFromSegmentNumber(segment.outletSegment());
        return segment.depth() - outlet_segment.depth();
    }

    void WellSegments::checkSegmentDepthConsistency(const std::string& well_name) const {
        constexpr double limit = 1. + std::numeric_limits<double>::epsilon();
        for (const auto& segment : this->m_segments) {
            const int segment_number = segment.segmentNumber();
            if (segment_number == 1) {
                continue; // not totally sure how to handle the top segment yet
            }
            const double segment_length = this->segmentLength(segment_number);
            const double segment_depth_change = this->segmentDepthChange(segment_number);
            if (std::abs(segment_depth_change) > limit * segment_length) {
                // TODO: output considering UNIT
                const std::string msg = fmt::format(" Segment {} of well {} has a depth change of {} meters,"
                                                    " while it has a length of {} meters, which is unphysical.",
                                                    segment_number, well_name, segment_depth_change, segment_length);
                OpmLog::warning(msg);
            }
        }
    }

    std::set<int> WellSegments::branches() const {
        std::set<int> bset;
        for (const auto& segment : this->m_segments)
            bset.insert( segment.branchNumber() );
        return bset;
    }


    std::vector<Segment> WellSegments::branchSegments(int branch) const {
        std::vector<Segment> segments;
        std::unordered_set<int> segment_set;
        for (const auto& segment : this->m_segments) {
            if (segment.branchNumber() == branch) {
                segments.push_back(segment);
                segment_set.insert( segment.segmentNumber() );
            }
        }


        std::size_t head_index = 0;
        while (head_index < segments.size()) {
            const auto& head_segment = segments[head_index];
            if (segment_set.count(head_segment.outletSegment()) != 0) {
                auto head_iter = std::find_if(std::next(segments.begin(), head_index), segments.end(),
                                              [&segment_set] (const Segment& segment) { return (segment_set.count(segment.outletSegment()) == 0); });

                if (head_iter == segments.end())
                    throw std::logic_error("Loop detected in branch/segment structure");
                std::swap( segments[head_index], *head_iter);
            }
            segment_set.erase( segments[head_index].segmentNumber() );
            head_index++;
        }

        return segments;
    }

    bool WellSegments::updateWSEGSICD(const std::vector<std::pair<int, SICD> >& sicd_pairs) {
        if (m_comp_pressure_drop == CompPressureDrop::H__) {
            const std::string msg = "to use spiral ICD segment you have to activate the frictional pressure drop calculation";
            throw std::runtime_error(msg);
        }

        for (const auto& pair_elem : sicd_pairs) {
            const int segment_number = pair_elem.first;
            const SICD& spiral_icd = pair_elem.second;
            Segment segment = this->getFromSegmentNumber(segment_number);
            segment.updateSpiralICD(spiral_icd);
            this->addSegment(segment);
        }

        return true;
    }

    bool WellSegments::updateWSEGVALV(const std::vector<std::pair<int, Valve> >& valve_pairs) {

        if (m_comp_pressure_drop == CompPressureDrop::H__) {
            const std::string msg = "to use WSEGVALV segment you have to activate the frictional pressure drop calculation";
            throw std::runtime_error(msg);
        }

        for (const auto& pair : valve_pairs) {
            const int segment_number = pair.first;
            const Valve& valve = pair.second;
            Segment segment = this->getFromSegmentNumber(segment_number);
            const double segment_length = this->segmentLength(segment_number);
            // this function can return bool
            segment.updateValve(valve, segment_length);
            this->addSegment(segment);
        }

        return true;
    }

    bool WellSegments::updateWSEGAICD(const std::vector<std::pair<int, AutoICD> >& aicd_pairs, const KeywordLocation& location) {
        if (m_comp_pressure_drop == CompPressureDrop::H__) {
            const std::string msg = fmt::format("to use Autonomous ICD segment with keyword {} "
                                                "at line {} in file {},\n"
                                                "you have to activate frictional pressure drop calculation in WELSEGS",
                                                location.keyword, location.lineno, location.filename);
            throw std::runtime_error(msg);
        }

        for (const auto& pair_elem : aicd_pairs) {
            const int segment_number = pair_elem.first;
            const AutoICD& auto_icd = pair_elem.second;
            Segment segment = this->getFromSegmentNumber(segment_number);
            segment.updateAutoICD(auto_icd);
            this->addSegment(segment);
        }

        return true;
    }


    bool WellSegments::operator!=( const WellSegments& rhs ) const {
        return !( *this == rhs );
    }


const std::string WellSegments::LengthDepthToString(LengthDepth enumValue) {
    switch (enumValue) {
    case LengthDepth::INC:
        return "INC";
    case LengthDepth::ABS:
        return "ABS";
    default:
        throw std::invalid_argument("unhandled LengthDepth value");
    }
}


WellSegments::LengthDepth WellSegments::LengthDepthFromString(const std::string& string_value ) {
    if (string_value == "INC") {
        return LengthDepth::INC;
    } else if (string_value == "ABS") {
        return LengthDepth::ABS;
    } else {
        throw std::invalid_argument("Unknown enum string_value: " + string_value + " for LengthDepth");
    }
}


const std::string WellSegments::CompPressureDropToString(CompPressureDrop enumValue) {
    switch (enumValue) {
    case CompPressureDrop::HFA:
        return "HFA";
    case CompPressureDrop::HF_:
        return "HF-";
    case CompPressureDrop::H__:
        return "H--";
    default:
        throw std::invalid_argument("unhandled CompPressureDrop value");
    }
}

WellSegments::CompPressureDrop WellSegments::CompPressureDropFromString( const std::string& string_value ) {

    if (string_value == "HFA") {
        return CompPressureDrop::HFA;
    } else if (string_value == "HF-") {
        return CompPressureDrop::HF_;
    } else if (string_value == "H--") {
        return CompPressureDrop::H__;
    } else {
        throw std::invalid_argument("Unknown enum string_value: " + string_value + " for CompPressureDrop");
    }
}

const std::string WellSegments::MultiPhaseModelToString(MultiPhaseModel enumValue) {
    switch (enumValue) {
    case MultiPhaseModel::HO:
        return "HO";
    case MultiPhaseModel::DF:
        return "DF";
    default:
        throw std::invalid_argument("unhandled MultiPhaseModel value");
    }
}

WellSegments::MultiPhaseModel WellSegments::MultiPhaseModelFromString(const std::string& string_value ) {

    if ((string_value == "HO") || (string_value == "H0")) {
        return MultiPhaseModel::HO;
    } else if (string_value == "DF") {
        return MultiPhaseModel::DF;
    } else {
        throw std::invalid_argument("Unknown enum string_value: " + string_value + " for MultiPhaseModel");
    }
}


void WellSegments::updatePerfLength(const WellConnections& connections) {
    for (auto& segment : this->m_segments) {
        auto perf_length = connections.segment_perf_length( segment.segmentNumber() );
        segment.updatePerfLength(perf_length);
    }
}

}
