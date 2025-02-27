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

namespace Opm
{

template <class Scalar, bool enableThermal>
void
GasPvtMultiplexer<Scalar, enableThermal>::initEnd()
{
    OPM_GAS_PVT_MULTIPLEXER_CALL(pvtImpl.initEnd(), break);
}

template <class Scalar, bool enableThermal>
unsigned
GasPvtMultiplexer<Scalar, enableThermal>::numRegions() const
{
    OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.numRegions());
}


template <class Scalar, bool enableThermal>
void
GasPvtMultiplexer<Scalar, enableThermal>::setVapPars(const Scalar par1, const Scalar par2)
{
    OPM_GAS_PVT_MULTIPLEXER_CALL(pvtImpl.setVapPars(par1, par2), break);
}


template <class Scalar, bool enableThermal>
Scalar
GasPvtMultiplexer<Scalar, enableThermal>::gasReferenceDensity(unsigned regionIdx)
{
    OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.gasReferenceDensity(regionIdx));
}


template <class Scalar, bool enableThermal>
Scalar
GasPvtMultiplexer<Scalar, enableThermal>::hVap(unsigned regionIdx) const
{
    OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.hVap(regionIdx));
}

#if HAVE_ECL_INPUT
template <class Scalar, bool enableThermal>
void
GasPvtMultiplexer<Scalar, enableThermal>::initFromState(const EclipseState& eclState, const Schedule& schedule)
{
    if (!eclState.runspec().phases().active(Phase::GAS))
        return;

    if (eclState.runspec().co2Storage())
        setApproach(GasPvtApproach::Co2Gas);
    else if (eclState.runspec().h2Storage())
        setApproach(GasPvtApproach::H2Gas);
    else if (enableThermal && (eclState.getSimulationConfig().isThermal() || eclState.getSimulationConfig().isTemp()))
        setApproach(GasPvtApproach::ThermalGas);
    else if (!eclState.getTableManager().getPvtgwTables().empty()
             && !eclState.getTableManager().getPvtgTables().empty())
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

// Helper function to keep the switch case tidy when constructing different pvts
template <class Scalar, bool enableThermal>
template <class ConcreteGasPvt>
std::unique_ptr<void, std::function<void(void*)>>
GasPvtMultiplexer<Scalar, enableThermal>::makeGasPvt()
{
    return UniqueVoidPtrWithDeleter(
        new ConcreteGasPvt,
        [this](void* ptr) { deleter(ptr); }
    );
}

template <class Scalar, bool enableThermal>
void GasPvtMultiplexer<Scalar, enableThermal>::
setApproach(GasPvtApproach gasPvtAppr)
{
    switch (gasPvtAppr) {
    case GasPvtApproach::DryGas:
        realGasPvt_ = makeGasPvt<DryGasPvt<Scalar>>();
        break;

    case GasPvtApproach::DryHumidGas:
        realGasPvt_ = makeGasPvt<DryHumidGasPvt<Scalar>>();
        break;

    case GasPvtApproach::WetHumidGas:
        realGasPvt_ = makeGasPvt<WetHumidGasPvt<Scalar>>();
        break;

    case GasPvtApproach::WetGas:
        realGasPvt_ = makeGasPvt<WetGasPvt<Scalar>>();
        break;

    case GasPvtApproach::ThermalGas:
        realGasPvt_ = makeGasPvt<GasPvtThermal<Scalar>>();
        break;

    case GasPvtApproach::Co2Gas:
        realGasPvt_ = makeGasPvt<Co2GasPvt<Scalar>>();
        break;

    case GasPvtApproach::H2Gas:
        realGasPvt_ = makeGasPvt<H2GasPvt<Scalar>>();
        break;

    case GasPvtApproach::NoGas:
        throw std::logic_error("Not implemented: Gas PVT of this deck!");
    }

    gasPvtApproach_ = gasPvtAppr;
}

// Helper template to create copies of PVT objects
template <class Scalar, bool enableThermal>
template <class ConcretePvt>
std::unique_ptr<void, std::function<void(void*)>>
GasPvtMultiplexer<Scalar, enableThermal>::
copyPvt(const std::unique_ptr<void, std::function<void(void*)>>& sourcePvt) {
    return UniqueVoidPtrWithDeleter(
        new ConcretePvt(*static_cast<const ConcretePvt*>(sourcePvt.get())),
        [this](void* ptr) { deleter(ptr); }
    );
}

template <class Scalar, bool enableThermal>
GasPvtMultiplexer<Scalar, enableThermal>&
GasPvtMultiplexer<Scalar, enableThermal>::operator=(const GasPvtMultiplexer<Scalar, enableThermal>& data)
{
    gasPvtApproach_ = data.gasPvtApproach_;
    switch (gasPvtApproach_) {
    case GasPvtApproach::DryGas:
        realGasPvt_ = copyPvt<DryGasPvt<Scalar>>(data.realGasPvt_);
        break;
    case GasPvtApproach::DryHumidGas:
        realGasPvt_ = copyPvt<DryHumidGasPvt<Scalar>>(data.realGasPvt_);
        break;
    case GasPvtApproach::WetHumidGas:
        realGasPvt_ = copyPvt<WetHumidGasPvt<Scalar>>(data.realGasPvt_);
        break;
    case GasPvtApproach::WetGas:
        realGasPvt_ = copyPvt<WetGasPvt<Scalar>>(data.realGasPvt_);
        break;
    case GasPvtApproach::ThermalGas:
        realGasPvt_ = copyPvt<GasPvtThermal<Scalar>>(data.realGasPvt_);
        break;
    case GasPvtApproach::Co2Gas:
        realGasPvt_ = copyPvt<Co2GasPvt<Scalar>>(data.realGasPvt_);
        break;
    case GasPvtApproach::H2Gas:
        realGasPvt_ = copyPvt<H2GasPvt<Scalar>>(data.realGasPvt_);
        break;
    default:
        break;
    }

    return *this;
}

template <class Scalar, bool enableThermal>
void GasPvtMultiplexer<Scalar, enableThermal>::deleter(void* ptr)
{
    switch (gasPvtApproach_) {
        case GasPvtApproach::DryGas: {
            delete static_cast<DryGasPvt<Scalar>*>(ptr);
            break;
        }
        case GasPvtApproach::DryHumidGas: {
            delete static_cast<DryHumidGasPvt<Scalar>*>(ptr);
            break;
        }
        case GasPvtApproach::WetHumidGas: {
            delete static_cast<WetHumidGasPvt<Scalar>*>(ptr);
            break;
        }
        case GasPvtApproach::WetGas: {
            delete static_cast<WetGasPvt<Scalar>*>(ptr);
            break;
        }
        case GasPvtApproach::ThermalGas: {
            delete static_cast<GasPvtThermal<Scalar>*>(ptr);
            break;
        }
        case GasPvtApproach::Co2Gas: {
            delete static_cast<Co2GasPvt<Scalar>*>(ptr);
            break;
        }
        case GasPvtApproach::H2Gas: {
            delete static_cast<H2GasPvt<Scalar>*>(ptr);
            break;
        }
        case GasPvtApproach::NoGas:
            break;
    }
}

template class GasPvtMultiplexer<double, false>;
template class GasPvtMultiplexer<double, true>;
template class GasPvtMultiplexer<float, false>;
template class GasPvtMultiplexer<float, true>;

} // namespace Opm
