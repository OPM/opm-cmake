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
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/common/ErrorMacros.hpp>

#if HAVE_ECL_INPUT
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Tables/FlatTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>
#endif

#include <string_view>

#include <fmt/format.h>

namespace Opm {

#if HAVE_ECL_INPUT
template <class Scalar, class IndexTraits, template<typename> typename Storage, template<typename> typename SmartPointer>
void BlackOilFluidSystem<Scalar,IndexTraits, Storage, SmartPointer>::
initFromState(const EclipseState& eclState, const Schedule& schedule)
{
    std::size_t num_regions = eclState.runspec().tabdims().getNumPVTTables();
    initBegin(num_regions);

    numActivePhases_ = 0;
    std::fill_n(&phaseIsActive_[0], numPhases, false);

    if (eclState.runspec().phases().active(Phase::OIL)) {
        phaseIsActive_[oilPhaseIdx] = true;
        ++numActivePhases_;
    }

    if (eclState.runspec().phases().active(Phase::GAS)) {
        phaseIsActive_[gasPhaseIdx] = true;
        ++numActivePhases_;
    }

    if (eclState.runspec().phases().active(Phase::WATER)) {
        phaseIsActive_[waterPhaseIdx] = true;
        ++numActivePhases_;
    }

    // this fluidsystem only supports one, two or three phases
    if (numActivePhases_ < 1 || numActivePhases_ > 3) {
        OPM_THROW(std::runtime_error,
                  fmt::format("Fluidsystem supports 1-3 phases, but {} is active\n",
                              numActivePhases_));
    }

    // set the surface conditions using the STCOND keyword
    surfaceTemperature = eclState.getTableManager().stCond().temperature;
    surfacePressure = eclState.getTableManager().stCond().pressure;

    // The reservoir temperature does not really belong into the table manager. TODO:
    // change this in opm-parser
    setReservoirTemperature(eclState.getTableManager().rtemp());

    setEnableDissolvedGas(eclState.getSimulationConfig().hasDISGAS());
    setEnableVaporizedOil(eclState.getSimulationConfig().hasVAPOIL());
    setEnableVaporizedWater(eclState.getSimulationConfig().hasVAPWAT());

    if (eclState.getSimulationConfig().hasDISGASW()) {
        if (eclState.runspec().co2Storage() || eclState.runspec().h2Storage())
            setEnableDissolvedGasInWater(eclState.getSimulationConfig().hasDISGASW());
        else if (eclState.runspec().co2Sol() || eclState.runspec().h2Sol()) {
            // For CO2SOL and H2SOL the dissolved gas in water is added in the solvent model
            // The HC gas is not allowed to dissolved into water.
            // For most HC gasses this is a resonable assumption.
            OpmLog::info("CO2SOL/H2SOL is activated together with DISGASW. \n"
                         "Only CO2/H2 is allowed to dissolve into water");
        } else
            OPM_THROW(std::runtime_error,
                      "DISGASW only supported in combination with CO2STORE/H2STORE or CO2SOL/H2SOL");
    }

    if (phaseIsActive(gasPhaseIdx)) {
        gasPvt_ = std::make_shared<GasPvt>();
        gasPvt_->initFromState(eclState, schedule);
    }

    if (phaseIsActive(oilPhaseIdx)) {
        oilPvt_ = std::make_shared<OilPvt>();
        oilPvt_->initFromState(eclState, schedule);
    }

    if (phaseIsActive(waterPhaseIdx)) {
        waterPvt_ = std::make_shared<WaterPvt>();
        waterPvt_->initFromState(eclState, schedule);
    }

    // set the reference densities of all PVT regions
    for (unsigned regionIdx = 0; regionIdx < num_regions; ++regionIdx) {
        setReferenceDensities(oilPvt_ ? oilPvt_->oilReferenceDensity(regionIdx) : 700.0,
                              waterPvt_ ? waterPvt_->waterReferenceDensity(regionIdx) : 1000.0,
                              gasPvt_ ? gasPvt_->gasReferenceDensity(regionIdx) : 2.0,
                              regionIdx);
    }

    // set default molarMass and mappings
    initEnd();

    // use molarMass of CO2 and Brine as default
    // when we are using the the CO2STORE option
    if (eclState.runspec().co2Storage()) {
        const Scalar salinity = eclState.getCo2StoreConfig().salinity();  // mass fraction
        for (unsigned regionIdx = 0; regionIdx < num_regions; ++regionIdx) {
            if (phaseIsActive(oilPhaseIdx)) // The oil component is used for the brine if OIL is active
                molarMass_[regionIdx][oilCompIdx] = BrineCo2Pvt<Scalar>::Brine::molarMass(salinity);
            if (phaseIsActive(waterPhaseIdx))
                molarMass_[regionIdx][waterCompIdx] = BrineCo2Pvt<Scalar>::Brine::molarMass(salinity);
            if (!phaseIsActive(gasPhaseIdx)) {
                OPM_THROW(std::runtime_error,
                          "CO2STORE requires gas phase\n");
            }
            molarMass_[regionIdx][gasCompIdx] = BrineCo2Pvt<Scalar>::CO2::molarMass();
        }
    }

    // Use molar mass of H2 and Brine as default in H2STORE keyword
    if (eclState.runspec().h2Storage()) {
        // Salinity in mass fraction
        const Scalar molality = eclState.getTableManager().salinity(); // mol/kg
        const Scalar MmNaCl = 58.44e-3; // molar mass of NaCl [kg/mol]
        const Scalar salinity = 1 / ( 1 + 1 / (molality*MmNaCl));
        for (unsigned regionIdx = 0; regionIdx < num_regions; ++regionIdx) {
            if (phaseIsActive(oilPhaseIdx)) // The oil component is used for the brine if OIL is active
                molarMass_[regionIdx][oilCompIdx] = BrineH2Pvt<Scalar>::Brine::molarMass(salinity);
            if (phaseIsActive(waterPhaseIdx))
                molarMass_[regionIdx][waterCompIdx] = BrineH2Pvt<Scalar>::Brine::molarMass(salinity);
            if (!phaseIsActive(gasPhaseIdx)) {
                OPM_THROW(std::runtime_error,
                          "H2STORE requires gas phase\n");
            }
            molarMass_[regionIdx][gasCompIdx] = BrineH2Pvt<Scalar>::H2::molarMass();
        }
    }

    // For co2storage and h2storage we dont have a concept of tables and should not spend time on
    // checking if we are at the saturated front
    setUseSaturatedTables(!(eclState.runspec().h2Storage() || eclState.runspec().co2Storage()));

    setEnableDiffusion(eclState.getSimulationConfig().isDiffusive());
    if (enableDiffusion()) {
        const auto& diffCoeffTables = eclState.getTableManager().getDiffusionCoefficientTable();
        if (!diffCoeffTables.empty()) {
            // if diffusion coefficient table is empty we relay on the PVT model to
            // to give us the coefficients.
            diffusionCoefficients_.resize(num_regions,{0,0,0,0,0,0,0,0,0});
            if (diffCoeffTables.size() != num_regions) {
                OPM_THROW(std::runtime_error,
                          fmt::format("Table sizes mismatch. DiffCoeffs: {}, NumRegions: {}\n",
                                      diffCoeffTables.size(), num_regions));
            }
            for (unsigned regionIdx = 0; regionIdx < num_regions; ++regionIdx) {
                const auto& diffCoeffTable = diffCoeffTables[regionIdx];
                molarMass_[regionIdx][oilCompIdx] = diffCoeffTable.oil_mw;
                molarMass_[regionIdx][gasCompIdx] = diffCoeffTable.gas_mw;
                setDiffusionCoefficient(diffCoeffTable.gas_in_gas, gasCompIdx, gasPhaseIdx, regionIdx);
                setDiffusionCoefficient(diffCoeffTable.oil_in_gas, oilCompIdx, gasPhaseIdx, regionIdx);
                setDiffusionCoefficient(diffCoeffTable.gas_in_oil, gasCompIdx, oilPhaseIdx, regionIdx);
                setDiffusionCoefficient(diffCoeffTable.oil_in_oil, oilCompIdx, oilPhaseIdx, regionIdx);
                if (diffCoeffTable.gas_in_oil_cross_phase > 0 || diffCoeffTable.oil_in_oil_cross_phase > 0) {
                    OPM_THROW(std::runtime_error,
                              "Cross phase diffusion is set in the deck, "
                              "but not implemented in Flow. "
                              "Please default DIFFC item 7 and item 8 "
                              "or set it to zero.");
                }
            }
        } else if ( (eclState.runspec().co2Storage() || eclState.runspec().h2Storage())
                && eclState.runspec().phases().active(Phase::GAS)
                && eclState.runspec().phases().active(Phase::WATER))
        {
            diffusionCoefficients_.resize(num_regions, {0,0,0,0,0,0,0,0,0});
            // diffusion coefficients can be set using DIFFCGAS and DIFFCWAT
            // for CO2STORE and H2STORE cases with gas + water
            const auto& diffCoeffWatTables = eclState.getTableManager().getDiffusionCoefficientWaterTable();
            if (!diffCoeffWatTables.empty()) {
                for (unsigned regionIdx = 0; regionIdx < num_regions; ++regionIdx) {
                    const auto& diffCoeffWatTable = diffCoeffWatTables[regionIdx];
                    setDiffusionCoefficient(diffCoeffWatTable.co2_in_water, gasCompIdx, waterPhaseIdx, regionIdx);
                    setDiffusionCoefficient(diffCoeffWatTable.h2o_in_water, waterCompIdx, waterPhaseIdx, regionIdx);
                }
            }
            const auto& diffCoeffGasTables = eclState.getTableManager().getDiffusionCoefficientGasTable();
            if (!diffCoeffGasTables.empty()) {
                for (unsigned regionIdx = 0; regionIdx < num_regions; ++regionIdx) {
                    const auto& diffCoeffGasTable = diffCoeffGasTables[regionIdx];
                    setDiffusionCoefficient(diffCoeffGasTable.co2_in_gas, gasCompIdx, gasPhaseIdx, regionIdx);
                    setDiffusionCoefficient(diffCoeffGasTable.h2o_in_gas, waterCompIdx, gasPhaseIdx, regionIdx);
                }
            }
        }
    }
}
#endif

#define INSTANTIATE_TYPE(T) \
    template<> unsigned char BlackOilFluidSystem<T, BlackOilDefaultIndexTraits, VectorWithDefaultAllocator, std::shared_ptr>::numActivePhases_ = 0; \
    template<> std::array<bool, BlackOilFluidSystem<T, BlackOilDefaultIndexTraits, VectorWithDefaultAllocator, std::shared_ptr>::numPhases> \
        BlackOilFluidSystem<T, BlackOilDefaultIndexTraits, VectorWithDefaultAllocator, std::shared_ptr>::phaseIsActive_ = {false, false, false}; \
    template<> std::array<short, BlackOilFluidSystem<T, BlackOilDefaultIndexTraits, VectorWithDefaultAllocator, std::shared_ptr>::numPhases> \
        BlackOilFluidSystem<T, BlackOilDefaultIndexTraits, VectorWithDefaultAllocator, std::shared_ptr>::activeToCanonicalPhaseIdx_ = {0, 1, 2}; \
    template<> std::array<short, BlackOilFluidSystem<T, BlackOilDefaultIndexTraits, VectorWithDefaultAllocator, std::shared_ptr>::numPhases> \
        BlackOilFluidSystem<T , BlackOilDefaultIndexTraits, VectorWithDefaultAllocator, std::shared_ptr>::canonicalToActivePhaseIdx_ = {0, 1, 2}; \
    template<> T BlackOilFluidSystem<T, BlackOilDefaultIndexTraits, VectorWithDefaultAllocator, std::shared_ptr>::surfaceTemperature = 0.0; \
    template<> T BlackOilFluidSystem<T, BlackOilDefaultIndexTraits, VectorWithDefaultAllocator, std::shared_ptr>::surfacePressure = 0.0; \
    template<> T BlackOilFluidSystem<T, BlackOilDefaultIndexTraits, VectorWithDefaultAllocator, std::shared_ptr>::reservoirTemperature_ = 0.0; \
    template<> bool BlackOilFluidSystem<T, BlackOilDefaultIndexTraits, VectorWithDefaultAllocator, std::shared_ptr>::enableDissolvedGas_ = true; \
    template<> bool BlackOilFluidSystem<T, BlackOilDefaultIndexTraits, VectorWithDefaultAllocator, std::shared_ptr>::enableDissolvedGasInWater_ = false; \
    template<> bool BlackOilFluidSystem<T, BlackOilDefaultIndexTraits, VectorWithDefaultAllocator, std::shared_ptr>::enableVaporizedOil_ = false; \
    template<> bool BlackOilFluidSystem<T, BlackOilDefaultIndexTraits, VectorWithDefaultAllocator, std::shared_ptr>::enableVaporizedWater_ = false; \
    template<> bool BlackOilFluidSystem<T, BlackOilDefaultIndexTraits, VectorWithDefaultAllocator, std::shared_ptr>::enableDiffusion_ = false; \
    template<> std::shared_ptr<OilPvtMultiplexer<T>> \
        BlackOilFluidSystem<T, BlackOilDefaultIndexTraits, VectorWithDefaultAllocator, std::shared_ptr>::oilPvt_ = {}; \
    template<> std::shared_ptr<GasPvtMultiplexer<T>> \
        BlackOilFluidSystem<T, BlackOilDefaultIndexTraits, VectorWithDefaultAllocator, std::shared_ptr>::gasPvt_ = {}; \
    template<> std::shared_ptr<WaterPvtMultiplexer<T>> \
        BlackOilFluidSystem<T, BlackOilDefaultIndexTraits, VectorWithDefaultAllocator, std::shared_ptr>::waterPvt_ = {}; \
    template<> std::vector<std::array<T, 3>> \
        BlackOilFluidSystem<T, BlackOilDefaultIndexTraits, VectorWithDefaultAllocator, std::shared_ptr>::referenceDensity_ = {}; \
    template<> std::vector<std::array<T, 3>> \
        BlackOilFluidSystem<T, BlackOilDefaultIndexTraits, VectorWithDefaultAllocator, std::shared_ptr>::molarMass_ = {}; \
    template<> std::vector<std::array<T, 9>> \
        BlackOilFluidSystem<T, BlackOilDefaultIndexTraits, VectorWithDefaultAllocator, std::shared_ptr>::diffusionCoefficients_ = {}; \
    template<> bool BlackOilFluidSystem<T, BlackOilDefaultIndexTraits, VectorWithDefaultAllocator, std::shared_ptr>::isInitialized_ = false; \
    template<> bool BlackOilFluidSystem<T, BlackOilDefaultIndexTraits, VectorWithDefaultAllocator, std::shared_ptr>::useSaturatedTables_ = false; \
    template<> bool BlackOilFluidSystem<T, BlackOilDefaultIndexTraits, VectorWithDefaultAllocator, std::shared_ptr>::enthalpy_eq_energy_ = false; \
    template class BlackOilFluidSystem<T,BlackOilDefaultIndexTraits, VectorWithDefaultAllocator, std::shared_ptr>;
    // IMPORTANT: The class must be instantiated after the template template specializations
    //    or else the static variable above will appear as undefined in the generated object file.

INSTANTIATE_TYPE(double)
INSTANTIATE_TYPE(float)

} // namespace Opm
