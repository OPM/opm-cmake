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
#include <opm/material/fluidmatrixinteractions/EclMaterialLawManagerSimple.hpp>

#include <opm/common/TimingMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Grid/SatfuncPropertyInitializers.hpp>

#include <opm/material/fluidmatrixinteractions/EclEpsGridProperties.hpp>
#include <opm/material/fluidstates/SimpleModularFluidState.hpp>

namespace Opm {

template<class TraitsT>
EclMaterialLawManagerSimple<TraitsT>::EclMaterialLawManagerSimple() = default;

template<class TraitsT>
EclMaterialLawManagerSimple<TraitsT>::~EclMaterialLawManagerSimple() = default;

template<class TraitsT>
void EclMaterialLawManagerSimple<TraitsT>::
initFromState(const EclipseState& eclState)
{
    // get the number of saturation regions and the number of cells in the deck
    const auto&  runspec       = eclState.runspec();
    const size_t numSatRegions = runspec.tabdims().getNumSatTables();

    const auto& ph = runspec.phases();
    this->hasGas = ph.active(Phase::GAS);
    this->hasOil = ph.active(Phase::OIL);
    this->hasWater = ph.active(Phase::WATER);

    // readGlobalEpsOptions_(eclState);
    // readGlobalHysteresisOptions_(eclState);
    readGlobalThreePhaseOptions_(runspec);

    // Read the end point scaling configuration (once per run).
    // gasOilConfig_ = std::make_shared<EclEpsConfig>();
    // oilWaterConfig_ = std::make_shared<EclEpsConfig>();
    // gasWaterConfig_ = std::make_shared<EclEpsConfig>();
    // gasOilConfig_->initFromState(eclState, EclTwoPhaseSystemType::GasOil);
    // oilWaterConfig_->initFromState(eclState, EclTwoPhaseSystemType::OilWater);
    // gasWaterConfig_->initFromState(eclState, EclTwoPhaseSystemType::GasWater);


    const auto& tables = eclState.getTableManager();

    {
        const auto& stone1exTables = tables.getStone1exTable();

        if (! stone1exTables.empty()) {
            stoneEtas_.clear();
            stoneEtas_.reserve(numSatRegions);

            for (const auto& table : stone1exTables) {
                stoneEtas_.push_back(table.eta);
            }
        }
        
        const auto& ppcwmaxTables = tables.getPpcwmax();
        this->enablePpcwmax_ = !ppcwmaxTables.empty();

        if (this->enablePpcwmax_) {
            maxAllowPc_.clear();
            modifySwl_.clear();

            maxAllowPc_.reserve(numSatRegions);
            modifySwl_.reserve(numSatRegions);

            for (const auto& table : ppcwmaxTables) {
                maxAllowPc_.push_back(table.max_cap_pres);
                modifySwl_.push_back(table.option);
            }
        }
    }

    this->unscaledEpsInfo_.resize(numSatRegions);

    if (this->hasGas + this->hasOil + this->hasWater == 1) {
        // Single-phase simulation.  Special case.  Nothing to do here.
        return;
    }

    // Multiphase simulation.  Common case.
    const auto tolcrit = runspec.saturationFunctionControls()
        .minimumRelpermMobilityThreshold();

    const auto rtep  = satfunc::getRawTableEndpoints(tables, ph, tolcrit);
    const auto rfunc = satfunc::getRawFunctionValues(tables, ph, rtep);

    for (unsigned satRegionIdx = 0; satRegionIdx < numSatRegions; ++satRegionIdx) {
        this->unscaledEpsInfo_[satRegionIdx]
            .extractUnscaled(rtep, rfunc, satRegionIdx);
    }

    // WAG hysteresis parameters per SATNUM.
    // if (eclState.runspec().hysterPar().activeWag()) {
    //     if (numSatRegions != eclState.getWagHysteresis().size())
    //         throw std::runtime_error("Inconsistent Wag-hysteresis data");
    //     wagHystersisConfig_.resize(numSatRegions);
    //     for (unsigned satRegionIdx = 0; satRegionIdx < numSatRegions; ++satRegionIdx) {
    //         wagHystersisConfig_[satRegionIdx] = std::make_shared<WagHysteresisConfig::
    //             WagHysteresisConfigRecord >(eclState.getWagHysteresis()[satRegionIdx]);
    //     }
    // }
}

template<class TraitsT>
void EclMaterialLawManagerSimple<TraitsT>::
initParamsForElements(const EclipseState& eclState, size_t numCompressedElems,
                      const std::function<std::vector<int>(const FieldPropsManager&, const std::string&, bool)>& fieldPropIntOnLeafAssigner,
                      const std::function<unsigned(unsigned)>& lookupIdxOnLevelZeroAssigner)
{
    InitParams initParams {*this, eclState, numCompressedElems};
    initParams.run(fieldPropIntOnLeafAssigner, lookupIdxOnLevelZeroAssigner);
}

// TODO: Better (proper?) handling of mixed wettability systems - see ecl kw OPTIONS switch 74
// Note: Without OPTIONS[74] the negative part of the Pcow curve is not scaled
template<class TraitsT>
std::pair<typename TraitsT::Scalar, bool> EclMaterialLawManagerSimple<TraitsT>::
applySwatinit(unsigned elemIdx,
              Scalar pcow,
              Scalar Sw)
{
    // Default is no SWATINIT scaling of the negative part of the Pcow curve, so look up saturation using the input Pcow curve
    if (pcow <= 0.0) {
        return {Sw, /*newSwatInit*/ true};
    }

    auto& elemScaledEpsInfo = oilWaterScaledEpsInfoDrainage_[elemIdx];
    if (Sw <= elemScaledEpsInfo.Swl)
        Sw = elemScaledEpsInfo.Swl;

    // specify a fluid state which only stores the saturations
    using FluidState = SimpleModularFluidState<Scalar,
                                                numPhases,
                                                /*numComponents=*/0,
                                                /*FluidSystem=*/void, /* -> don't care */
                                                /*storePressure=*/false,
                                                /*storeTemperature=*/false,
                                                /*storeComposition=*/false,
                                                /*storeFugacity=*/false,
                                                /*storeSaturation=*/true,
                                                /*storeDensity=*/false,
                                                /*storeViscosity=*/false,
                                                /*storeEnthalpy=*/false>;
    FluidState fs;
    fs.setSaturation(waterPhaseIdx, Sw);
    fs.setSaturation(gasPhaseIdx, 0);
    fs.setSaturation(oilPhaseIdx, 0);
    std::array<Scalar, numPhases> pc = { 0 };
    MaterialLaw::capillaryPressures(pc, materialLawParams(elemIdx), fs);
    Scalar pcowAtSw = pc[oilPhaseIdx] - pc[waterPhaseIdx];
    constexpr const Scalar pcowAtSwThreshold = 1.0e-6; //Pascal

    // avoid division by very small number and avoid negative PCW at connate Sw
    // (look up saturation on input Pcow curve in this case)
    if (pcowAtSw < pcowAtSwThreshold) {
        return {Sw, /*newSwatInit*/ true};
    }

    // Sufficiently positive value, continue with max. capillary pressure (PCW) scaling to honor SWATINIT value
    Scalar newMaxPcow = elemScaledEpsInfo.maxPcow * (pcow/pcowAtSw);

    // Limit max. capillary pressure with PPCWMAX
    bool newSwatInit = false;
    int satRegionIdx = satnumRegionIdx(elemIdx);
    if (enablePpcwmax() && (newMaxPcow > maxAllowPc_[satRegionIdx])) {
        // Two options in PPCWMAX to modify connate Sw or not.  In both cases, init. Sw needs to be
        // re-calculated (done in opm-simulators)
        newSwatInit = true;
        if (modifySwl_[satRegionIdx] == false) {
            // Max. cap. pressure set to PCWO in PPCWMAX
            elemScaledEpsInfo.maxPcow = maxAllowPc_[satRegionIdx];
        }
        else {
            // Max. cap. pressure remains unscaled and connate Sw is set to SWATINIT value
            elemScaledEpsInfo.Swl = Sw;
        }
    }
    // Max. cap. pressure adjusted from SWATINIT data
    else
        elemScaledEpsInfo.maxPcow = newMaxPcow;

    // auto& elemEclEpsScalingPoints = oilWaterScaledEpsPointsDrainage(elemIdx);
    // elemEclEpsScalingPoints.init(elemScaledEpsInfo,
    //                                 *oilWaterEclEpsConfig_,
    //                                 EclTwoPhaseSystemType::OilWater);

    return {Sw, newSwatInit};
}

template<class TraitsT>
void
EclMaterialLawManagerSimple<TraitsT>::applyRestartSwatInit(const unsigned elemIdx,
                                                     const Scalar   maxPcow)
{
    // Maximum capillary pressure adjusted from SWATINIT data.

    auto& elemScaledEpsInfo =
        this->oilWaterScaledEpsInfoDrainage_[elemIdx];

    elemScaledEpsInfo.maxPcow = maxPcow;

    // this->oilWaterScaledEpsPointsDrainage(elemIdx)
    //     .init(elemScaledEpsInfo,
    //           *this->oilWaterEclEpsConfig_,
    //           EclTwoPhaseSystemType::OilWater);
}

template<class TraitsT>
const typename EclMaterialLawManagerSimple<TraitsT>::MaterialLawParams&
EclMaterialLawManagerSimple<TraitsT>::
connectionMaterialLawParams(unsigned satRegionIdx, unsigned elemIdx) const
{
    MaterialLawParams& mlp = const_cast<MaterialLawParams&>(materialLawParams_[elemIdx]);

    if (enableHysteresis())
        OpmLog::warning("Warning: Using non-default satnum regions for connection is not tested in combination with hysteresis");

    auto& realParams = mlp;
    if (realParams.approach() == EclTwoPhaseApproach::GasOil) {
        realParams.setGasOilParams(gasOilEffectiveParamVector_[satRegionIdx]);
    }
    else if (realParams.approach() == EclTwoPhaseApproach::GasWater) {
        realParams.setGasWaterParams(gasWaterEffectiveParamVector_[satRegionIdx]);
    }
    else if (realParams.approach() == EclTwoPhaseApproach::OilWater) {
        realParams.setOilWaterParams(oilWaterEffectiveParamVector_[satRegionIdx]);
    }
    return mlp;
}

template<class TraitsT>
int EclMaterialLawManagerSimple<TraitsT>::
getKrnumSatIdx(unsigned elemIdx, FaceDir::DirEnum facedir) const
{
    using Dir = FaceDir::DirEnum;
    const std::vector<int>* array = nullptr;
    switch(facedir) {
    case Dir::XPlus:
      array = &krnumXArray_;
      break;
    case Dir::YPlus:
      array = &krnumYArray_;
      break;
    case Dir::ZPlus:
      array = &krnumZArray_;
      break;
    default:
      throw std::runtime_error("Unknown face direction");
    }
    if (array->size() > 0) {
      return (*array)[elemIdx];
    }
    else {
      return satnumRegionArray_[elemIdx];
    }
}

template<class TraitsT>
void EclMaterialLawManagerSimple<TraitsT>::
oilWaterHysteresisParams(Scalar& soMax,
                         Scalar& swMax,
                         Scalar& swMin,
                         unsigned elemIdx) const
{
    if (!enableHysteresis())
        throw std::runtime_error("Cannot get hysteresis parameters if hysteresis not enabled.");

    const auto& params = materialLawParams(elemIdx);
    MaterialLaw::oilWaterHysteresisParams(soMax, swMax, swMin, params);
}

template<class TraitsT>
void EclMaterialLawManagerSimple<TraitsT>::
setOilWaterHysteresisParams(const Scalar& soMax,
                            const Scalar& swMax,
                            const Scalar& swMin,
                            unsigned elemIdx)
{
    if (!enableHysteresis())
        throw std::runtime_error("Cannot set hysteresis parameters if hysteresis not enabled.");

    auto& params = materialLawParams(elemIdx);
    MaterialLaw::setOilWaterHysteresisParams(soMax, swMax, swMin, params);
}

template<class TraitsT>
void EclMaterialLawManagerSimple<TraitsT>::
gasOilHysteresisParams(Scalar& sgmax,
                       Scalar& shmax,
                       Scalar& somin,
                       unsigned elemIdx) const
{
    if (!enableHysteresis())
        throw std::runtime_error("Cannot get hysteresis parameters if hysteresis not enabled.");

    const auto& params = materialLawParams(elemIdx);
    MaterialLaw::gasOilHysteresisParams(sgmax, shmax, somin, params);
}

template<class TraitsT>
void EclMaterialLawManagerSimple<TraitsT>::
setGasOilHysteresisParams(const Scalar& sgmax,
                          const Scalar& shmax,
                          const Scalar& somin,
                          unsigned elemIdx)
{
    if (!enableHysteresis())
        throw std::runtime_error("Cannot set hysteresis parameters if hysteresis not enabled.");

    auto& params = materialLawParams(elemIdx);
    MaterialLaw::setGasOilHysteresisParams(sgmax, shmax, somin, params);
}

// template<class TraitsT>
// EclEpsScalingPoints<typename TraitsT::Scalar>&
// EclMaterialLawManagerSimple<TraitsT>::
// oilWaterScaledEpsPointsDrainage(unsigned elemIdx)
// {
//     auto& materialParams = materialLawParams_[elemIdx];
//     auto& realParams = materialParams;
//     return realParams.oilWaterParams().drainageParams().scaledPoints();
// }

template<class TraitsT>
const typename EclMaterialLawManagerSimple<TraitsT>::MaterialLawParams& EclMaterialLawManagerSimple<TraitsT>::
materialLawParamsFunc_(unsigned elemIdx, FaceDir::DirEnum facedir) const
{
    using Dir = FaceDir::DirEnum;
    if (dirMaterialLawParams_) {
        switch(facedir) {
            case Dir::XMinus:
            case Dir::XPlus:
                return dirMaterialLawParams_->materialLawParamsX_[elemIdx];
            case Dir::YMinus:
            case Dir::YPlus:
                return dirMaterialLawParams_->materialLawParamsY_[elemIdx];
            case Dir::ZMinus:
            case Dir::ZPlus:
                return dirMaterialLawParams_->materialLawParamsZ_[elemIdx];
            default:
                throw std::runtime_error("Unexpected face direction");
        }
    }
    else {
        return materialLawParams_[elemIdx];
    }
}

// template<class TraitsT>
// void EclMaterialLawManagerSimple<TraitsT>::
// readGlobalEpsOptions_(const EclipseState& eclState)
// {
//     oilWaterEclEpsConfig_ = std::make_shared<EclEpsConfig>();
//     oilWaterEclEpsConfig_->initFromState(eclState, EclTwoPhaseSystemType::OilWater);

//     enableEndPointScaling_ = eclState.getTableManager().hasTables("ENKRVD");
// }

// template<class TraitsT>
// void EclMaterialLawManagerSimple<TraitsT>::
// readGlobalHysteresisOptions_(const EclipseState& state)
// {
//     hysteresisConfig_ = std::make_shared<EclHysteresisConfig>();
//     hysteresisConfig_->initFromState(state.runspec());
// }

template<class TraitsT>
void EclMaterialLawManagerSimple<TraitsT>::
readGlobalThreePhaseOptions_(const Runspec& runspec)
{
    bool gasEnabled = runspec.phases().active(Phase::GAS);
    bool oilEnabled = runspec.phases().active(Phase::OIL);
    bool waterEnabled = runspec.phases().active(Phase::WATER);

    int numEnabled =
        (gasEnabled?1:0)
        + (oilEnabled?1:0)
        + (waterEnabled?1:0);

    if (numEnabled == 0) {
        throw std::runtime_error("At least one fluid phase must be enabled. (Is: "+std::to_string(numEnabled)+")");
    } else if (numEnabled == 1) {
        threePhaseApproach_ = EclMultiplexerApproach::OnePhase;
    } else if ( numEnabled == 2) {
        threePhaseApproach_ = EclMultiplexerApproach::TwoPhase;
        if (!gasEnabled)
            twoPhaseApproach_ = EclTwoPhaseApproach::OilWater;
        else if (!oilEnabled)
            twoPhaseApproach_ = EclTwoPhaseApproach::GasWater;
        else if (!waterEnabled)
            twoPhaseApproach_ = EclTwoPhaseApproach::GasOil;
    }
    else {
        assert(numEnabled == 3);

        threePhaseApproach_ = EclMultiplexerApproach::Default;
        const auto& satctrls = runspec.saturationFunctionControls();
        if (satctrls.krModel() == SatFuncControls::ThreePhaseOilKrModel::Stone2)
            threePhaseApproach_ = EclMultiplexerApproach::Stone2;
        else if (satctrls.krModel() == SatFuncControls::ThreePhaseOilKrModel::Stone1)
            threePhaseApproach_ = EclMultiplexerApproach::Stone1;
    }
}


template class EclMaterialLawManagerSimple<ThreePhaseMaterialTraits<double,0,1,2>>;
template class EclMaterialLawManagerSimple<ThreePhaseMaterialTraits<float,0,1,2>>;

} // namespace Opm
