// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/

#include <config.h>
#include <opm/material/fluidsystems/blackoilpvt/WaterPvtMultiplexer.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Runspec.hpp>

namespace Opm {

template<class Scalar, bool enableThermal, bool enableBrine>
WaterPvtMultiplexer<Scalar,enableThermal,enableBrine>::
~WaterPvtMultiplexer()
{
    switch (approach_) {
    case WaterPvtApproach::ConstantCompressibilityWater: {
        delete &getRealPvt<WaterPvtApproach::ConstantCompressibilityWater>();
        break;
    }
    case WaterPvtApproach::ConstantCompressibilityBrine: {
        delete &getRealPvt<WaterPvtApproach::ConstantCompressibilityBrine>();
        break;
    }
    case WaterPvtApproach::ThermalWater: {
        delete &getRealPvt<WaterPvtApproach::ThermalWater>();
        break;
    }
    case WaterPvtApproach::BrineCo2: {
        delete &getRealPvt<WaterPvtApproach::BrineCo2>();
        break;
    }
    case WaterPvtApproach::BrineH2: {
        delete &getRealPvt<WaterPvtApproach::BrineH2>();
        break;
    }
    case WaterPvtApproach::NoWater:
        break;
    }
}

#if HAVE_ECL_INPUT
template<class Scalar, bool enableThermal, bool enableBrine>
void WaterPvtMultiplexer<Scalar,enableThermal,enableBrine>::
initFromState(const EclipseState& eclState, const Schedule& schedule)
{
    if (!eclState.runspec().phases().active(Phase::WATER))
        return;

    // The co2Storage option both works with oil + gas
    // and water/brine + gas
    if (eclState.runspec().co2Storage() || eclState.runspec().co2Sol())
        setApproach(WaterPvtApproach::BrineCo2);
    else if (eclState.runspec().h2Storage() || eclState.runspec().h2Sol())
        setApproach(WaterPvtApproach::BrineH2);
    else if (enableThermal && (eclState.getSimulationConfig().isThermal() || eclState.getSimulationConfig().isTemp()))
        setApproach(WaterPvtApproach::ThermalWater);
    else if (!eclState.getTableManager().getPvtwTable().empty())
        setApproach(WaterPvtApproach::ConstantCompressibilityWater);
    else if (enableBrine && !eclState.getTableManager().getPvtwSaltTables().empty())
        setApproach(WaterPvtApproach::ConstantCompressibilityBrine);

    OPM_WATER_PVT_MULTIPLEXER_CALL(pvtImpl.initFromState(eclState, schedule), break);
}
#endif

template<class Scalar, bool enableThermal, bool enableBrine>
void WaterPvtMultiplexer<Scalar,enableThermal,enableBrine>::
initEnd()
{
    OPM_WATER_PVT_MULTIPLEXER_CALL(pvtImpl.initEnd(), break);
}


template<class Scalar, bool enableThermal, bool enableBrine>
unsigned WaterPvtMultiplexer<Scalar,enableThermal,enableBrine>::
numRegions() const
{
    OPM_WATER_PVT_MULTIPLEXER_CALL(return pvtImpl.numRegions());
}

template<class Scalar, bool enableThermal, bool enableBrine>
void WaterPvtMultiplexer<Scalar,enableThermal,enableBrine>::
setVapPars(const Scalar par1, const Scalar par2)
{
    OPM_WATER_PVT_MULTIPLEXER_CALL(pvtImpl.setVapPars(par1, par2), break);
}

template<class Scalar, bool enableThermal, bool enableBrine>
Scalar WaterPvtMultiplexer<Scalar,enableThermal,enableBrine>::
waterReferenceDensity(unsigned regionIdx) const
{
    OPM_WATER_PVT_MULTIPLEXER_CALL(return pvtImpl.waterReferenceDensity(regionIdx));
}

template<class Scalar, bool enableThermal, bool enableBrine>
Scalar WaterPvtMultiplexer<Scalar,enableThermal,enableBrine>::
hVap(unsigned regionIdx) const
{
    OPM_WATER_PVT_MULTIPLEXER_CALL(return pvtImpl.hVap(regionIdx));
}

template<class Scalar, bool enableThermal, bool enableBrine>
void WaterPvtMultiplexer<Scalar,enableThermal,enableBrine>::
setApproach(WaterPvtApproach appr)
{
    switch (appr) {
    case WaterPvtApproach::ConstantCompressibilityWater:
        realWaterPvt_ = new ConstantCompressibilityWaterPvt<Scalar>;
        break;

    case WaterPvtApproach::ConstantCompressibilityBrine:
        realWaterPvt_ = new ConstantCompressibilityBrinePvt<Scalar>;
        break;

    case WaterPvtApproach::ThermalWater:
        realWaterPvt_ = new WaterPvtThermal<Scalar, enableBrine>;
        break;

    case WaterPvtApproach::BrineCo2:
        realWaterPvt_ = new BrineCo2Pvt<Scalar>;
        break;

    case WaterPvtApproach::BrineH2:
        realWaterPvt_ = new BrineH2Pvt<Scalar>;
        break;

    case WaterPvtApproach::NoWater:
        throw std::logic_error("Not implemented: Water PVT of this deck!");
    }

    approach_ = appr;
}

template<class Scalar, bool enableThermal, bool enableBrine>
WaterPvtMultiplexer<Scalar,enableThermal,enableBrine>&
WaterPvtMultiplexer<Scalar,enableThermal,enableBrine>::
operator=(const WaterPvtMultiplexer<Scalar,enableThermal,enableBrine>& data)
{
    approach_ = data.approach_;
    switch (approach_) {
    case WaterPvtApproach::ConstantCompressibilityWater:
        realWaterPvt_ = new ConstantCompressibilityWaterPvt<Scalar>(*static_cast<const ConstantCompressibilityWaterPvt<Scalar>*>(data.realWaterPvt_));
        break;
    case WaterPvtApproach::ConstantCompressibilityBrine:
        realWaterPvt_ = new ConstantCompressibilityBrinePvt<Scalar>(*static_cast<const ConstantCompressibilityBrinePvt<Scalar>*>(data.realWaterPvt_));
        break;
    case WaterPvtApproach::ThermalWater:
        realWaterPvt_ = new WaterPvtThermal<Scalar, enableBrine>(*static_cast<const WaterPvtThermal<Scalar, enableBrine>*>(data.realWaterPvt_));
        break;
    case WaterPvtApproach::BrineCo2:
        realWaterPvt_ = new BrineCo2Pvt<Scalar>(*static_cast<const BrineCo2Pvt<Scalar>*>(data.realWaterPvt_));
        break;
    case WaterPvtApproach::BrineH2:
        realWaterPvt_ = new BrineH2Pvt<Scalar>(*static_cast<const BrineH2Pvt<Scalar>*>(data.realWaterPvt_));
        break;
    default:
        break;
    }

    return *this;
}

template class WaterPvtMultiplexer<double,false,false>;
template class WaterPvtMultiplexer<double,true,false>;
template class WaterPvtMultiplexer<double,false,true>;
template class WaterPvtMultiplexer<double,true,true>;
template class WaterPvtMultiplexer<float,false,false>;
template class WaterPvtMultiplexer<float,true,false>;
template class WaterPvtMultiplexer<float,false,true>;
template class WaterPvtMultiplexer<float,true,true>;

} // namespace Opm
