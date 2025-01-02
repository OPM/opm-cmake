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
#include <opm/material/fluidsystems/blackoilpvt/GasPvtMultiplexer.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>

namespace Opm {

template <class Scalar, bool enableThermal>
GasPvtMultiplexer<Scalar,enableThermal>::
~GasPvtMultiplexer()
{
    switch (gasPvtApproach_) {
    case GasPvtApproach::DryGas: {
        delete &getRealPvt<GasPvtApproach::DryGas>();
        break;
    }
    case GasPvtApproach::DryHumidGas: {
        delete &getRealPvt<GasPvtApproach::DryHumidGas>();
        break;
    }
    case GasPvtApproach::WetHumidGas: {
        delete &getRealPvt<GasPvtApproach::WetHumidGas>();
        break;
    }
    case GasPvtApproach::WetGas: {
        delete &getRealPvt<GasPvtApproach::WetGas>();
        break;
    }
    case GasPvtApproach::ThermalGas: {
        delete &getRealPvt<GasPvtApproach::ThermalGas>();
        break;
    }
    case GasPvtApproach::Co2Gas: {
        delete &getRealPvt<GasPvtApproach::Co2Gas>();
        break;
    }
    case GasPvtApproach::H2Gas: {
        delete &getRealPvt<GasPvtApproach::H2Gas>();
        break;
    }
    case GasPvtApproach::NoGas:
        break;
    }
}

template <class Scalar, bool enableThermal>
void GasPvtMultiplexer<Scalar,enableThermal>::
initEnd()
{
    OPM_GAS_PVT_MULTIPLEXER_CALL(pvtImpl.initEnd(), break);
}

template <class Scalar, bool enableThermal>
unsigned GasPvtMultiplexer<Scalar,enableThermal>::
numRegions() const
{
    OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.numRegions());
}


template <class Scalar, bool enableThermal>
void GasPvtMultiplexer<Scalar,enableThermal>::
setVapPars(const Scalar par1, const Scalar par2)
{
    OPM_GAS_PVT_MULTIPLEXER_CALL(pvtImpl.setVapPars(par1, par2), break);
}


template <class Scalar, bool enableThermal>
Scalar GasPvtMultiplexer<Scalar,enableThermal>::
gasReferenceDensity(unsigned regionIdx)
{
    OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.gasReferenceDensity(regionIdx));
}


template <class Scalar, bool enableThermal>
Scalar GasPvtMultiplexer<Scalar,enableThermal>::
hVap(unsigned regionIdx) const
{
    OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.hVap(regionIdx));
}

#if HAVE_ECL_INPUT
template <class Scalar, bool enableThermal>
void GasPvtMultiplexer<Scalar,enableThermal>::
initFromState(const EclipseState& eclState, const Schedule& schedule)
{
    if (!eclState.runspec().phases().active(Phase::GAS))
        return;

    if (eclState.runspec().co2Storage())
        setApproach(GasPvtApproach::Co2Gas);
    else if (eclState.runspec().h2Storage())
        setApproach(GasPvtApproach::H2Gas);
    else if (enableThermal && (eclState.getSimulationConfig().isThermal() || eclState.getSimulationConfig().isTemp()))
        setApproach(GasPvtApproach::ThermalGas);
    else if (!eclState.getTableManager().getPvtgwTables().empty() &&
            !eclState.getTableManager().getPvtgTables().empty())
        setApproach(GasPvtApproach::WetHumidGas);
    else if (!eclState.getTableManager().getPvtgTables().empty())
        setApproach(GasPvtApproach::WetGas);
    else if (eclState.getTableManager().hasTables("PVDG"))
        setApproach(GasPvtApproach::DryGas);
    else if (!eclState.getTableManager().getPvtgwTables().empty())
        setApproach(GasPvtApproach::DryHumidGas);


    OPM_GAS_PVT_MULTIPLEXER_CALL(pvtImpl.initFromState(eclState, schedule), break);
}
#endif

template <class Scalar, bool enableThermal>
void GasPvtMultiplexer<Scalar,enableThermal>::
setApproach(GasPvtApproach gasPvtAppr)
{
    switch (gasPvtAppr) {
    case GasPvtApproach::DryGas:
        realGasPvt_ = new DryGasPvt<Scalar>;
        break;

    case GasPvtApproach::DryHumidGas:
        realGasPvt_ = new DryHumidGasPvt<Scalar>;
        break;

    case GasPvtApproach::WetHumidGas:
        realGasPvt_ = new WetHumidGasPvt<Scalar>;
        break;

    case GasPvtApproach::WetGas:
        realGasPvt_ = new WetGasPvt<Scalar>;
        break;

    case GasPvtApproach::ThermalGas:
        realGasPvt_ = new GasPvtThermal<Scalar>;
        break;

    case GasPvtApproach::Co2Gas:
        realGasPvt_ = new Co2GasPvt<Scalar>;
        break;

    case GasPvtApproach::H2Gas:
        realGasPvt_ = new H2GasPvt<Scalar>;
        break;

    case GasPvtApproach::NoGas:
        throw std::logic_error("Not implemented: Gas PVT of this deck!");
    }

    gasPvtApproach_ = gasPvtAppr;
}

template <class Scalar, bool enableThermal>
GasPvtMultiplexer<Scalar,enableThermal>&
GasPvtMultiplexer<Scalar,enableThermal>::
operator=(const GasPvtMultiplexer<Scalar,enableThermal>& data)
{
    gasPvtApproach_ = data.gasPvtApproach_;
    switch (gasPvtApproach_) {
    case GasPvtApproach::DryGas:
        realGasPvt_ = new DryGasPvt<Scalar>(*static_cast<const DryGasPvt<Scalar>*>(data.realGasPvt_));
        break;
    case GasPvtApproach::DryHumidGas:
        realGasPvt_ = new DryHumidGasPvt<Scalar>(*static_cast<const DryHumidGasPvt<Scalar>*>(data.realGasPvt_));
        break;
    case GasPvtApproach::WetHumidGas:
        realGasPvt_ = new WetHumidGasPvt<Scalar>(*static_cast<const WetHumidGasPvt<Scalar>*>(data.realGasPvt_));
        break;
    case GasPvtApproach::WetGas:
        realGasPvt_ = new WetGasPvt<Scalar>(*static_cast<const WetGasPvt<Scalar>*>(data.realGasPvt_));
        break;
    case GasPvtApproach::ThermalGas:
        realGasPvt_ = new GasPvtThermal<Scalar>(*static_cast<const GasPvtThermal<Scalar>*>(data.realGasPvt_));
        break;
    case GasPvtApproach::Co2Gas:
        realGasPvt_ = new Co2GasPvt<Scalar>(*static_cast<const Co2GasPvt<Scalar>*>(data.realGasPvt_));
        break;
case GasPvtApproach::H2Gas:
        realGasPvt_ = new H2GasPvt<Scalar>(*static_cast<const H2GasPvt<Scalar>*>(data.realGasPvt_));
        break;
    default:
        break;
    }

    return *this;
}

template class GasPvtMultiplexer<double,false>;
template class GasPvtMultiplexer<double,true>;
template class GasPvtMultiplexer<float,false>;
template class GasPvtMultiplexer<float,true>;

} // namespace Opm
