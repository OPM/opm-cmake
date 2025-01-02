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
#include <opm/material/fluidsystems/blackoilpvt/DryGasPvt.hpp>

#include <opm/common/ErrorMacros.hpp>

#if HAVE_ECL_INPUT
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PvdgTable.hpp>
#endif

#include <fmt/format.h>
namespace Opm {

#if HAVE_ECL_INPUT
template<class Scalar>
void DryGasPvt<Scalar>::
initFromState(const EclipseState& eclState, const Schedule&)
{
    const auto& pvdgTables = eclState.getTableManager().getPvdgTables();
    const auto& densityTable = eclState.getTableManager().getDensityTable();

    if (pvdgTables.size() != densityTable.size()) {
        OPM_THROW(std::runtime_error,
                  fmt::format("Table sizes mismatch. PVDG: {}, DensityTable: {}\n",
                              pvdgTables.size(), densityTable.size()));
    }

    std::size_t regions = pvdgTables.size();
    setNumRegions(regions);

    for (unsigned regionIdx = 0; regionIdx < regions; ++regionIdx) {
        Scalar rhoRefO = densityTable[regionIdx].oil;
        Scalar rhoRefG = densityTable[regionIdx].gas;
        Scalar rhoRefW = densityTable[regionIdx].water;

        setReferenceDensities(regionIdx, rhoRefO, rhoRefG, rhoRefW);

        // determine the molar masses of the components
        constexpr Scalar p = 1.01325e5; // surface pressure, [Pa]
        constexpr Scalar T = 273.15 + 15.56; // surface temperature, [K]
        constexpr Scalar MO = 175e-3; // [kg/mol]
        Scalar MG = Constants<Scalar>::R*T*rhoRefG / p; // [kg/mol], consequence of the ideal gas law
        constexpr Scalar MW = 18.0e-3; // [kg/mol]
        // TODO (?): the molar mass of the components can possibly specified
        // explicitly in the deck.
        setMolarMasses(regionIdx, MO, MG, MW);

        const auto& pvdgTable = pvdgTables.getTable<PvdgTable>(regionIdx);

        // say 99.97% of all time: "premature optimization is the root of all
        // evil". Eclipse does this "optimization" for apparently no good reason!
        std::vector<Scalar> invB(pvdgTable.numRows());
        const auto& Bg = pvdgTable.getFormationFactorColumn();
        for (unsigned i = 0; i < Bg.size(); ++i) {
            invB[i] = 1.0 / Bg[i];
        }

        std::size_t numSamples = invB.size();
        inverseGasB_[regionIdx].setXYArrays(numSamples, pvdgTable.getPressureColumn(), invB);
        gasMu_[regionIdx].setXYArrays(numSamples, pvdgTable.getPressureColumn(), pvdgTable.getViscosityColumn());
    }

    initEnd();
}
#endif

template<class Scalar>
void DryGasPvt<Scalar>::setNumRegions(std::size_t numRegions)
{
    gasReferenceDensity_.resize(numRegions);
    inverseGasB_.resize(numRegions);
    inverseGasBMu_.resize(numRegions);
    gasMu_.resize(numRegions);
}

template<class Scalar>
void DryGasPvt<Scalar>::
setGasFormationVolumeFactor(unsigned regionIdx,
                            const SamplingPoints& samplePoints)
{
    SamplingPoints tmp(samplePoints);
    auto it = tmp.begin();
    const auto& endIt = tmp.end();
    for (; it != endIt; ++ it)
        std::get<1>(*it) = 1.0/std::get<1>(*it);

    inverseGasB_[regionIdx].setContainerOfTuples(tmp);
    assert(inverseGasB_[regionIdx].monotonic());
}

template<class Scalar>
void DryGasPvt<Scalar>::initEnd()
{
    // calculate the final 2D functions which are used for interpolation.
    std::size_t regions = gasMu_.size();
    for (unsigned regionIdx = 0; regionIdx < regions; ++ regionIdx) {
        // calculate the table which stores the inverse of the product of the gas
        // formation volume factor and the gas viscosity
        const auto& gasmu = gasMu_[regionIdx];
        const auto& invGasB = inverseGasB_[regionIdx];
        assert(gasmu.numSamples() == invGasB.numSamples());

        std::vector<Scalar> pressureValues(gasmu.numSamples());
        std::vector<Scalar> invGasBMuValues(gasmu.numSamples());
        for (unsigned pIdx = 0; pIdx < gasmu.numSamples(); ++pIdx) {
            pressureValues[pIdx] = invGasB.xAt(pIdx);
            invGasBMuValues[pIdx] = invGasB.valueAt(pIdx) * (1.0 / gasmu.valueAt(pIdx));
        }

        inverseGasBMu_[regionIdx].setXYContainers(pressureValues, invGasBMuValues);
    }
}

template class DryGasPvt<double>;
template class DryGasPvt<float>;

} // namespace Opm
