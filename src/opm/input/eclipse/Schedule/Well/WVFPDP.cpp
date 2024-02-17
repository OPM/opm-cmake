/*
  Copyright 2023 Equinor.

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

#include <opm/input/eclipse/Schedule/Well/WVFPDP.hpp>
#include <opm/io/eclipse/rst/well.hpp>
#include <opm/output/eclipse/VectorItems/well.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/W.hpp>
#include <opm/input/eclipse/Deck/DeckRecord.hpp>

#include <string>
#include <vector>

namespace Opm {

    WVFPDP WVFPDP::serializationTestObject()
    {
        WVFPDP result;
        result.m_dp = 1.23;
        result.m_fp = 0.456;
        result.m_active = false;

        return result;
    }

    bool WVFPDP::operator==(const WVFPDP& other) const {
        return (m_dp == other.m_dp)
            && (m_fp == other.m_fp)
            && (m_active == other.m_active);
    }

    void WVFPDP::update(const DeckRecord& record) {
        m_dp = record.getItem<ParserKeywords::WVFPDP::DELTA_P>().getSIDouble(0);
        m_fp = record.getItem<ParserKeywords::WVFPDP::LOSS_SCALING_FACTOR>().get<double>(0);
        this->m_active = true;
    }

    void WVFPDP::update(const RestartIO::RstWell& rst_well) {
        this->m_dp = rst_well.vfp_bhp_adjustment;
        this->m_fp = rst_well.vfp_bhp_scaling_factor;
        // \Note: not totally sure we should do this
        this->m_active = true;
    }

    bool WVFPDP::active() const {
        return this->m_active;
    }

    double WVFPDP::getPressureLoss(double bhp_tab, double thp_limit) const {
        auto dpt = bhp_tab - thp_limit;
        return m_dp + (m_fp-1)*dpt;
    }

    bool WVFPDP::operator!=(const WVFPDP& other) const {
        return !(*this == other);
    }

}
