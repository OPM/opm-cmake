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
#include <opm/material/fluidsystems/blackoilpvt/Co2GasPvt.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>

namespace Opm {

template<class Scalar>
void Co2GasPvt<Scalar>::
initFromState(const EclipseState& eclState, const Schedule&)
{

    setEnableVaporizationWater(eclState.getSimulationConfig().hasVAPOIL() || eclState.getSimulationConfig().hasVAPWAT());
    setActivityModelSalt(eclState.getTableManager().actco2s());

    bool co2sol = eclState.runspec().co2Sol();
    if (!co2sol && (eclState.getTableManager().hasTables("PVDG") ||
        !eclState.getTableManager().getPvtgTables().empty())) {
        OpmLog::warning("CO2STORE is enabled but PVDG or PVTG is in the deck. \n"
                        "CO2 PVT properties are computed based on the Span-Wagner "
                        "pvt model and PVDG/PVTG input is ignored.");
    }

    // We only supported single pvt region for the co2-brine module
    size_t numRegions = 1;
    setNumRegions(numRegions);
    size_t regionIdx = 0;

    // Currently we only support constant salinity
    const Scalar molality = eclState.getTableManager().salinity(); // mol/kg
    const Scalar MmNaCl = 58.44e-3; // molar mass of NaCl [kg/mol]
    // convert to mass fraction
    salinity_[regionIdx] = 1 / ( 1 + 1 / (molality*MmNaCl));

    if (co2sol) {
        brineReferenceDensity_[regionIdx] = eclState.getTableManager().getDensityTable()[regionIdx].water;
        Scalar T_ref = eclState.getTableManager().stCond().temperature;
        Scalar P_ref = eclState.getTableManager().stCond().pressure;
        // Throw an error if STCOND is not (T, p) = (15.56 C, 1 atm) = (288.71 K, 1.01325e5 Pa)
        if (T_ref != Scalar(288.71) || P_ref != Scalar(1.01325e5)) {
            OPM_THROW(std::runtime_error, "CO2SOL can only be used with default values for STCOND!");
        }
        gasReferenceDensity_[regionIdx] = CO2::gasDensity(T_ref, P_ref, extrapolate);
    } else if (!eclState.getTableManager().getDensityTable().empty()) {
        OpmLog::info("CO2STORE is enabled and DENSITY is given in the deck. \n"
                     "The surface density from DENSITY is used.");
        brineReferenceDensity_[regionIdx] = eclState.getTableManager().getDensityTable()[regionIdx].water;
        gasReferenceDensity_[regionIdx] = eclState.getTableManager().getDensityTable()[regionIdx].gas;
    } else {
        Scalar T_ref = eclState.getTableManager().stCond().temperature;
        Scalar P_ref = eclState.getTableManager().stCond().pressure;
        // Throw an error if STCOND is not (T, p) = (15.56 C, 1 atm) = (288.71 K, 1.01325e5 Pa)
        if (T_ref != Scalar(288.71) || P_ref != Scalar(1.01325e5)) {
            OPM_THROW(std::runtime_error, "CO2STORE can only be used with default values for STCOND!");
        }
        gasReferenceDensity_[regionIdx] = CO2::gasDensity(T_ref, P_ref, extrapolate);
        brineReferenceDensity_[regionIdx] = Brine::liquidDensity(T_ref, P_ref, salinity_[regionIdx], extrapolate);
        OpmLog::info("CO2STORE is enabled and DENSITY is not given in the deck. \n "
            "The surface densities are computed using the reference pressure ( " + std::to_string(P_ref) 
            + "Pa) and the reference temperature (" +  std::to_string(T_ref) + "K)." 
            + "The surface density of CO2 is  " + std::to_string(gasReferenceDensity_[regionIdx])
            + "kg/m3 \n The surface density of Brine is  " + std::to_string(brineReferenceDensity_[regionIdx])
            + "kg/m3"
            );
    }
    initEnd();
}

template class Co2GasPvt<double>;
template class Co2GasPvt<float>;

} // namespace Opm
