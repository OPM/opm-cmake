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
#include <opm/material/fluidsystems/blackoilpvt/OilPvtMultiplexer.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Runspec.hpp>

namespace Opm {

template <class Scalar, bool enableThermal>
OilPvtMultiplexer<Scalar,enableThermal>::
~OilPvtMultiplexer()
{
    switch (approach_) {
    case OilPvtApproach::LiveOil: {
        delete &getRealPvt<OilPvtApproach::LiveOil>();
        break;
    }
    case OilPvtApproach::DeadOil: {
        delete &getRealPvt<OilPvtApproach::DeadOil>();
        break;
    }
    case OilPvtApproach::ConstantCompressibilityOil: {
        delete &getRealPvt<OilPvtApproach::ConstantCompressibilityOil>();
        break;
    }
    case OilPvtApproach::ThermalOil: {
        delete &getRealPvt<OilPvtApproach::ThermalOil>();
        break;
    }
    case OilPvtApproach::BrineCo2: {
        delete &getRealPvt<OilPvtApproach::BrineCo2>();
        break;
    }
    case OilPvtApproach::BrineH2: {
        delete &getRealPvt<OilPvtApproach::BrineH2>();
        break;
    }
    case OilPvtApproach::NoOil:
        break;
    }
}

#if HAVE_ECL_INPUT
template <class Scalar, bool enableThermal>
void OilPvtMultiplexer<Scalar,enableThermal>::
initFromState(const EclipseState& eclState, const Schedule& schedule)
{
    if (!eclState.runspec().phases().active(Phase::OIL))
        return;

    // The co2Storage option both works with oil + gas
    // and water/brine + gas
    if (eclState.runspec().co2Storage())
        setApproach(OilPvtApproach::BrineCo2);
    else if (eclState.runspec().h2Storage())
        setApproach(OilPvtApproach::BrineH2);
    else if (enableThermal && (eclState.getSimulationConfig().isThermal() || eclState.getSimulationConfig().isTemp()))
        setApproach(OilPvtApproach::ThermalOil);
    else if (!eclState.getTableManager().getPvcdoTable().empty())
        setApproach(OilPvtApproach::ConstantCompressibilityOil);
    else if (eclState.getTableManager().hasTables("PVDO"))
        setApproach(OilPvtApproach::DeadOil);
    else if (!eclState.getTableManager().getPvtoTables().empty())
        setApproach(OilPvtApproach::LiveOil);

    OPM_OIL_PVT_MULTIPLEXER_CALL(pvtImpl.initFromState(eclState, schedule), break);
}
#endif

template <class Scalar, bool enableThermal>
void OilPvtMultiplexer<Scalar,enableThermal>::
initEnd()
{
    OPM_OIL_PVT_MULTIPLEXER_CALL(pvtImpl.initEnd(), break);
}

template <class Scalar, bool enableThermal>
unsigned OilPvtMultiplexer<Scalar,enableThermal>::
numRegions() const
{
    OPM_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.numRegions());
}

template <class Scalar, bool enableThermal>
void OilPvtMultiplexer<Scalar,enableThermal>::
setVapPars(const Scalar par1, const Scalar par2)
{
    OPM_OIL_PVT_MULTIPLEXER_CALL(pvtImpl.setVapPars(par1, par2), break);
}

template <class Scalar, bool enableThermal>
Scalar OilPvtMultiplexer<Scalar,enableThermal>::
oilReferenceDensity(unsigned regionIdx) const
{
    OPM_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.oilReferenceDensity(regionIdx));
}

template <class Scalar, bool enableThermal>
void OilPvtMultiplexer<Scalar,enableThermal>::
setApproach(OilPvtApproach appr)
{
    switch (appr) {
    case OilPvtApproach::LiveOil:
        realOilPvt_ = new LiveOilPvt<Scalar>;
        break;

    case OilPvtApproach::DeadOil:
        realOilPvt_ = new DeadOilPvt<Scalar>;
        break;

    case OilPvtApproach::ConstantCompressibilityOil:
        realOilPvt_ = new ConstantCompressibilityOilPvt<Scalar>;
        break;

    case OilPvtApproach::ThermalOil:
        realOilPvt_ = new OilPvtThermal<Scalar>;
        break;

    case OilPvtApproach::BrineCo2:
        realOilPvt_ = new BrineCo2Pvt<Scalar>;
        break;

    case OilPvtApproach::BrineH2:
        realOilPvt_ = new BrineH2Pvt<Scalar>;
        break;

    case OilPvtApproach::NoOil:
        throw std::logic_error("Not implemented: Oil PVT of this deck!");
    }

    approach_ = appr;
}

template <class Scalar, bool enableThermal>
OilPvtMultiplexer<Scalar,enableThermal>&
OilPvtMultiplexer<Scalar,enableThermal>::
operator=(const OilPvtMultiplexer<Scalar,enableThermal>& data)
{
    approach_ = data.approach_;
    switch (approach_) {
    case OilPvtApproach::ConstantCompressibilityOil:
        realOilPvt_ = new ConstantCompressibilityOilPvt<Scalar>(*static_cast<const ConstantCompressibilityOilPvt<Scalar>*>(data.realOilPvt_));
        break;
    case OilPvtApproach::DeadOil:
        realOilPvt_ = new DeadOilPvt<Scalar>(*static_cast<const DeadOilPvt<Scalar>*>(data.realOilPvt_));
        break;
    case OilPvtApproach::LiveOil:
        realOilPvt_ = new LiveOilPvt<Scalar>(*static_cast<const LiveOilPvt<Scalar>*>(data.realOilPvt_));
        break;
    case OilPvtApproach::ThermalOil:
        realOilPvt_ = new OilPvtThermal<Scalar>(*static_cast<const OilPvtThermal<Scalar>*>(data.realOilPvt_));
        break;
    case OilPvtApproach::BrineCo2:
        realOilPvt_ = new BrineCo2Pvt<Scalar>(*static_cast<const BrineCo2Pvt<Scalar>*>(data.realOilPvt_));
        break;
    case OilPvtApproach::BrineH2:
        realOilPvt_ = new BrineH2Pvt<Scalar>(*static_cast<const BrineH2Pvt<Scalar>*>(data.realOilPvt_));
        break;
    default:
        break;
    }

    return *this;
}

template class OilPvtMultiplexer<double,false>;
template class OilPvtMultiplexer<double,true>;
template class OilPvtMultiplexer<float,false>;
template class OilPvtMultiplexer<float,true>;

} // namespace Opm
