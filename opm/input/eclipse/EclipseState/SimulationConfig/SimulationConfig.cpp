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

#include <opm/input/eclipse/EclipseState/SimulationConfig/SimulationConfig.hpp>

#include <opm/input/eclipse/EclipseState/Grid/FieldPropsManager.hpp>

#include <opm/input/eclipse/EclipseState/SimulationConfig/BCConfig.hpp>
#include <opm/input/eclipse/EclipseState/SimulationConfig/DatumDepth.hpp>
#include <opm/input/eclipse/EclipseState/SimulationConfig/RockConfig.hpp>
#include <opm/input/eclipse/EclipseState/SimulationConfig/ThresholdPressure.hpp>

#include <opm/input/eclipse/Deck/DeckSection.hpp>

#include <opm/input/eclipse/Parser/ParserKeywords/C.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/D.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/G.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/N.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/P.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/T.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/V.hpp>

#include <stdexcept>

// The internalization of the CPR keyword has been temporarily disabled,
// suddenly decks with 'CPR' in the summary section turned up. Keywords with
// section aware keyword semantics is currently not handled by the parser.
//
// When the CPR is added again the following keyword configuration must be
// added:
//
//   {"name" : "CPR" , "sections" : ["RUNSPEC"], "size": 1 }

namespace Opm {

    SimulationConfig::SimulationConfig(bool restart,
                                       const Deck& deck,
                                       const FieldPropsManager& fp)
        : m_ThresholdPressure { restart, deck, fp }
        , m_bcconfig          { deck }
        , m_rock_config       { deck, fp }
        , m_datum_depth       { SOLUTIONSection { deck } }
    {
        if (DeckSection::hasRUNSPEC(deck)) {
            const RUNSPECSection runspec(deck);

            if (runspec.hasKeyword<ParserKeywords::CPR>()) {
                const auto& cpr = runspec.get<ParserKeywords::CPR>().back();
                if (cpr.size() > 0) {
                    throw std::invalid_argument {
                        "ERROR: In the RUNSPEC section the CPR keyword "
                        "should contain EXACTLY one empty record."
                    };
                }

                m_useCPR = true;
            }

            if (runspec.hasKeyword<ParserKeywords::NONNC>()) {
                const auto& nonnc = runspec.get<ParserKeywords::NONNC>().back();
                if (nonnc.size() > 0) {
                    throw std::invalid_argument {
                        "ERROR: In the RUNSPEC section the NONNC keyword "
                        "should contain EXACTLY one empty record."
                    };
                }

                m_useNONNC = true;
            }

            if (runspec.hasKeyword<ParserKeywords::DISGAS>()) {
                m_DISGAS = true;
            }

            if (runspec.hasKeyword<ParserKeywords::DISGASW>()) {
                m_DISGASW = true;
            }

            if (runspec.hasKeyword<ParserKeywords::VAPOIL>()) {
                m_VAPOIL = true;
            }

            if (runspec.hasKeyword<ParserKeywords::VAPWAT>()) {
                m_VAPWAT = true;
            }

            if (runspec.hasKeyword<ParserKeywords::DIFFUSE>()) {
                m_diffuse = true;
            }

            this->m_isThermal = runspec.hasKeyword<ParserKeywords::THERMAL>();
            this->m_isTemp = runspec.hasKeyword<ParserKeywords::TEMP>();

            this->m_useEnthalpy = runspec.hasKeyword<ParserKeywords::THERMAL>();

            if(runspec.hasKeyword<ParserKeywords::TEMP>()){
                this->m_useEnthalpy = false;
                if(runspec.hasKeyword<ParserKeywords::THERMAL>()){
                    throw std::invalid_argument {
                        "ERROR: In the RUNSPEC section the BOTH TEMP and THERMAL keyword "
                        "only one should be spesified"
                    };
                }
            }

            if (runspec.hasKeyword<ParserKeywords::PRECSALT>()) {
                m_PRECSALT = true;
            }
        }
    }

    SimulationConfig SimulationConfig::serializationTestObject()
    {
        SimulationConfig result;

        result.m_ThresholdPressure = ThresholdPressure::serializationTestObject();
        result.m_bcconfig = BCConfig::serializationTestObject();
        result.m_rock_config = RockConfig::serializationTestObject();
        result.m_datum_depth = DatumDepth::serializationTestObjectUserDefined();
        result.m_useCPR = true;
        result.m_useNONNC = true;
        result.m_DISGAS = true;
        result.m_DISGASW = true;
        result.m_VAPOIL = true;
        result.m_VAPWAT = true;
        result.m_isThermal = true;
        result.m_isTemp = true;
        result.m_diffuse = true;
        result.m_PRECSALT = true;

        return result;
    }

    const RockConfig& SimulationConfig::rock_config() const
    {
        return this->m_rock_config;
    }

    const ThresholdPressure& SimulationConfig::getThresholdPressure() const
    {
        return m_ThresholdPressure;
    }

    const BCConfig& SimulationConfig::bcconfig() const
    {
        return m_bcconfig;
    }

    const DatumDepth& SimulationConfig::datumDepths() const
    {
        return this->m_datum_depth;
    }

    bool SimulationConfig::useThresholdPressure() const
    {
        return m_ThresholdPressure.active();
    }

    bool SimulationConfig::useCPR() const
    {
        return m_useCPR;
    }

    bool SimulationConfig::useNONNC() const
    {
        return m_useNONNC;
    }

    bool SimulationConfig::hasDISGAS() const
    {
        return m_DISGAS;
    }

    bool SimulationConfig::hasDISGASW() const
    {
        return m_DISGASW;
    }

    bool SimulationConfig::hasVAPOIL() const
    {
        return m_VAPOIL;
    }

    bool SimulationConfig::hasVAPWAT() const
    {
        return m_VAPWAT;
    }

    bool SimulationConfig::isThermal() const
    {
        return this->m_isThermal;
    }

    bool SimulationConfig::isTemp() const
    {
        return this->m_isTemp;
    }

    bool SimulationConfig::useEnthalpy() const
    {
        return this->m_useEnthalpy;
    }

    bool SimulationConfig::isDiffusive() const
    {
        return this->m_diffuse;
    }

    bool SimulationConfig::hasPRECSALT() const
    {
        return m_PRECSALT;
    }

    bool SimulationConfig::operator==(const SimulationConfig& data) const
    {
        return (this->getThresholdPressure() == data.getThresholdPressure())
            && (this->bcconfig() == data.bcconfig())
            && (this->rock_config() == data.rock_config())
            && (this->datumDepths() == data.datumDepths())
            && (this->useCPR() == data.useCPR())
            && (this->useNONNC() == data.useNONNC())
            && (this->hasDISGAS() == data.hasDISGAS())
            && (this->hasDISGASW() == data.hasDISGASW())
            && (this->hasVAPOIL() == data.hasVAPOIL())
            && (this->hasVAPWAT() == data.hasVAPWAT())
            && (this->isTemp() == data.isTemp())
            && (this->isThermal() == data.isThermal())
            && (this->isDiffusive() == data.isDiffusive())
            && (this->hasPRECSALT() == data.hasPRECSALT())
            ;
    }

    bool SimulationConfig::rst_cmp(const SimulationConfig& full_config,
                                   const SimulationConfig& rst_config)
    {
        return ThresholdPressure::rst_cmp(full_config.getThresholdPressure(), rst_config.getThresholdPressure())
            && (full_config.bcconfig() == rst_config.bcconfig())
            && (full_config.rock_config() == rst_config.rock_config())
            && (full_config.datumDepths() == rst_config.datumDepths())
            && (full_config.useCPR() == rst_config.useCPR())
            && (full_config.useNONNC() == rst_config.useNONNC())
            && (full_config.hasDISGAS() == rst_config.hasDISGAS())
            && (full_config.hasDISGASW() == rst_config.hasDISGASW())
            && (full_config.hasVAPOIL() == rst_config.hasVAPOIL())
            && (full_config.hasVAPWAT() == rst_config.hasVAPWAT())
            && (full_config.isTemp() == rst_config.isTemp())
            && (full_config.isThermal() == rst_config.isThermal())
            && (full_config.isDiffusive() == rst_config.isDiffusive())
            && (full_config.hasPRECSALT() == rst_config.hasPRECSALT())
            ;
    }

} // namespace Opm
