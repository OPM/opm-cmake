/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.

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

#ifndef SEGMENTSET_HPP_HEADER_INCLUDED
#define SEGMENTSET_HPP_HEADER_INCLUDED

#include <map>
#include <set>
#include <vector>

#include <opm/input/eclipse/Schedule/MSW/Segment.hpp>

namespace Opm {
    class SICD;
    class AutoICD;
    class Valve;
    class WellConnections;
}

namespace Opm {

    enum class WellSegmentCompPressureDrop {
          HFA = 0,
          HF_ = 1,
          H__ = 2
    };

    class DeckKeyword;
    class KeywordLocation;

    class WellSegments {
    public:
        enum class LengthDepth {
            INC = 0,
            ABS = 1
        };
        static const std::string LengthDepthToString(LengthDepth enumValue);
        static LengthDepth LengthDepthFromString(const std::string& stringValue);

        using CompPressureDrop = WellSegmentCompPressureDrop;

        static const std::string CompPressureDropToString(CompPressureDrop enumValue);
        static CompPressureDrop CompPressureDropFromString(const std::string& stringValue);


        enum class MultiPhaseModel {
            HO = 0,
            DF = 1
        };
        static const std::string MultiPhaseModelToString(MultiPhaseModel enumValue);
        static MultiPhaseModel MultiPhaseModelFromString(const std::string& stringValue);

        WellSegments() = default;
        WellSegments(CompPressureDrop compDrop,
                     const std::vector<Segment>& segments);
        explicit WellSegments(const DeckKeyword& keyword);
        void loadWELSEGS( const DeckKeyword& welsegsKeyword);

        static WellSegments serializationTestObject();

        std::size_t size() const;
        bool empty() const;
        int maxSegmentID() const;
        int maxBranchID() const;
        double depthTopSegment() const;
        double lengthTopSegment() const;
        double volumeTopSegment() const;

        CompPressureDrop compPressureDrop() const;

        // mapping the segment number to the index in the vector of segments
        int segmentNumberToIndex(const int segment_number) const;



        const Segment& getFromSegmentNumber(const int segment_number) const;

        const Segment& operator[](size_t idx) const;
        void orderSegments();
        void updatePerfLength(const WellConnections& connections);

        bool operator==( const WellSegments& ) const;
        bool operator!=( const WellSegments& ) const;

        double segmentLength(const int segment_number) const;
        double segmentDepthChange(const int segment_number) const;
        std::vector<Segment> branchSegments(int branch) const;
        std::set<int> branches() const;

        // it returns true if there is no error encountered during the update
        bool updateWSEGSICD(const std::vector<std::pair<int, SICD> >& sicd_pairs);

        bool updateWSEGVALV(const std::vector<std::pair<int, Valve> >& valve_pairs);
        bool updateWSEGAICD(const std::vector<std::pair<int, AutoICD> >& aicd_pairs, const KeywordLocation& location);
        const std::vector<Segment>::const_iterator begin() const;
        const std::vector<Segment>::const_iterator end() const;

        void checkSegmentDepthConsistency(const std::string& well_name) const;

        template<class Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(m_comp_pressure_drop);
            serializer(m_segments);
            serializer(segment_number_to_index);
        }

    private:
        void processABS();
        void processINC(double depth_top, double length_top);
        void process(LengthDepth length_depth, double depth_top, double length_top);
        void addSegment(const Segment& new_segment);
        void addSegment(const int segment_number,
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
                        const double node_y);
        const Segment& topSegment() const;

        // components of the pressure drop to be included
        CompPressureDrop m_comp_pressure_drop{CompPressureDrop::HFA};
        // There are other three properties for segment related to thermal conduction,
        // while they are not supported by the keyword at the moment.

        std::vector< Segment > m_segments{};
        // the mapping from the segment number to the
        // storage index in the vector
        std::map<int, int> segment_number_to_index{};
    };
}

#endif
