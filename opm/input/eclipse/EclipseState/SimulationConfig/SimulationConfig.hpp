/*
  Copyright 2015 Statoil ASA.

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

#ifndef OPM_SIMULATION_CONFIG_HPP
#define OPM_SIMULATION_CONFIG_HPP

#include <opm/input/eclipse/EclipseState/SimulationConfig/BCConfig.hpp>
#include <opm/input/eclipse/EclipseState/SimulationConfig/DatumDepth.hpp>
#include <opm/input/eclipse/EclipseState/SimulationConfig/RockConfig.hpp>
#include <opm/input/eclipse/EclipseState/SimulationConfig/ThresholdPressure.hpp>

namespace Opm {

    class Deck;
    class FieldPropsManager;

} // namespace Opm

namespace Opm {

    class SimulationConfig
    {
    public:
        SimulationConfig() = default;
        SimulationConfig(bool restart,
                         const Deck& deck,
                         const FieldPropsManager& fp);

        static SimulationConfig serializationTestObject();

        const RockConfig& rock_config() const;
        const ThresholdPressure& getThresholdPressure() const;
        const BCConfig& bcconfig() const;
        const DatumDepth& datumDepths() const;

        bool useThresholdPressure() const;
        bool useCPR() const;
        bool useNONNC() const;
        bool hasDISGAS() const;
        bool hasDISGASW() const;
        bool hasVAPOIL() const;
        bool hasVAPWAT() const;
        bool isThermal() const;
        bool isTemp() const;
        bool isDiffusive() const;
        bool hasPRECSALT() const;

        bool operator==(const SimulationConfig& data) const;
        static bool rst_cmp(const SimulationConfig& full_config,
                            const SimulationConfig& rst_config);

        template<class Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(m_ThresholdPressure);
            serializer(m_bcconfig);
            serializer(m_rock_config);
            serializer(this->m_datum_depth);
            serializer(m_useCPR);
            serializer(m_useNONNC);
            serializer(m_DISGAS);
            serializer(m_DISGASW);
            serializer(m_VAPOIL);
            serializer(m_VAPWAT);
            serializer(m_isTemp);
            serializer(m_isThermal);
            serializer(m_diffuse);
            serializer(m_PRECSALT);
        }

    private:
        friend class EclipseState;

        ThresholdPressure m_ThresholdPressure{};
        BCConfig m_bcconfig{};
        RockConfig m_rock_config{};
        DatumDepth m_datum_depth{};
        bool m_useCPR{false};
        bool m_useNONNC{false};
        bool m_DISGAS{false};
        bool m_DISGASW{false};
        bool m_VAPOIL{false};
        bool m_VAPWAT{false};
        bool m_isTemp{false};
        bool m_isThermal{false};
        bool m_diffuse{false};
        bool m_PRECSALT{false};
    };

} // namespace Opm

#endif // OPM_SIMULATION_CONFIG_HPP
