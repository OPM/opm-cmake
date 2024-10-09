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
/*!
 * \file
 * \copydoc Opm::EclTwoPhaseMaterialParams
 */
#ifndef OPM_ECL_TWO_PHASE_MATERIAL_PARAMS_HPP
#define OPM_ECL_TWO_PHASE_MATERIAL_PARAMS_HPP

#include <memory>

#include <opm/common/utility/gpuDecorators.hpp>
#include <opm/material/common/EnsureFinalized.hpp>

namespace Opm {

enum class EclTwoPhaseApproach {
    GasOil,
    OilWater,
    GasWater
};

/*!
 * \brief Implementation for the parameters required by the material law for two-phase
 *        simulations.
 *
 * Essentially, this class just stores the two parameter objects for
 * the twophase capillary pressure laws.
 */
template<class Traits, class GasOilParamsT, class OilWaterParamsT, class GasWaterParamsT, template <typename> class PtrT = std::shared_ptr>
class EclTwoPhaseMaterialParams : public EnsureFinalized
{
    using Scalar = typename Traits::Scalar;
    enum { numPhases = 3 };
public:
    using EnsureFinalized :: finalize;

    using GasOilParams = GasOilParamsT;
    using OilWaterParams = OilWaterParamsT;
    using GasWaterParams = GasWaterParamsT;


    // make it possible to opt out of using shared pointers to stay within valid cuda c++
    using GasOilParamsPtr = PtrT<GasOilParams>;
    using OilWaterParamsPtr = PtrT<OilWaterParams>;
    using GasWaterParamsPtr = PtrT<GasWaterParams>;

    /*!
     * \brief The default constructor.
     */
    OPM_HOST_DEVICE EclTwoPhaseMaterialParams()
    {
    }

    OPM_HOST_DEVICE void setApproach(EclTwoPhaseApproach newApproach)
    { approach_ = newApproach; }

    OPM_HOST_DEVICE EclTwoPhaseApproach approach() const
    { return approach_; }

    /*!
     * \brief The parameter object for the gas-oil twophase law.
     */
    OPM_HOST_DEVICE const GasOilParams& gasOilParams() const
    { EnsureFinalized::check(); return *gasOilParams_; }

    /*!
     * \brief The parameter object for the gas-oil twophase law.
     */
    OPM_HOST_DEVICE GasOilParams& gasOilParams()
    { EnsureFinalized::check(); return *gasOilParams_; }

    /*!
     * \brief Set the parameter object for the gas-oil twophase law.
     */
    OPM_HOST_DEVICE void setGasOilParams(GasOilParamsPtr val)
    { gasOilParams_ = val; }

    /*!
     * \brief The parameter object for the oil-water twophase law.
     */
    OPM_HOST_DEVICE const OilWaterParams& oilWaterParams() const
    { EnsureFinalized::check(); return *oilWaterParams_; }

    /*!
     * \brief The parameter object for the oil-water twophase law.
     */
    OPM_HOST_DEVICE OilWaterParams& oilWaterParams()
    { EnsureFinalized::check(); return *oilWaterParams_; }

    /*!
     * \brief Set the parameter object for the oil-water twophase law.
     */
    OPM_HOST_DEVICE void setOilWaterParams(OilWaterParamsPtr val)
    { oilWaterParams_ = val; }

  /*!
     * \brief The parameter object for the gas-water twophase law.
     */
    OPM_HOST_DEVICE const GasWaterParams& gasWaterParams() const
    { EnsureFinalized::check(); return *gasWaterParams_; }

    /*!
     * \brief The parameter object for the gas-water twophase law.
     */
    OPM_HOST_DEVICE GasWaterParams& gasWaterParams()
    { EnsureFinalized::check(); return *gasWaterParams_; }

    /*!
     * \brief Set the parameter object for the gas-water twophase law.
     */
    OPM_HOST_DEVICE void setGasWaterParams(GasWaterParamsPtr val)
    { gasWaterParams_ = val; }

    template<class Serializer>
    OPM_HOST_DEVICE void serializeOp(Serializer& serializer)
    {
        // This is for restart serialization.
        // Only dynamic state in the parameters need to be stored.
        serializer(*gasOilParams_);
        serializer(*oilWaterParams_);
        serializer(*gasWaterParams_);
    }

    OPM_HOST_DEVICE void setSwl(Scalar) {}

private:
    EclTwoPhaseApproach approach_{EclTwoPhaseApproach::GasOil};

    GasOilParamsPtr gasOilParams_;
    OilWaterParamsPtr oilWaterParams_;
    GasWaterParamsPtr gasWaterParams_;
};

} // namespace Opm

// namespace Opm::gpuistl{

// /// @brief this function is intented to make a GPU friendly view of the EclTwoPhaseMaterialParams
// /// @tparam TraitsT the same traits as in EclTwoPhaseMaterialParams
// /// @tparam ContainerType typically const gpuBuffer<scalarType>
// /// @tparam ViewType  typically gpuView<const scalarType>
// /// @param params the parameters object instansiated with gpuBuffers or similar
// /// @return the GPU view of the GPU EclTwoPhaseMaterialParams object
// template <class TraitsT, class GasOilParamsT, class OilWaterParamsT, class GasWwaterParamsT>
// EclTwoPhaseMaterialParams<TraitsT, GasOilParamsT, OilWaterParamsT, GasWwaterParamsT, false> make_view(const EclTwoPhaseMaterialParams<TraitsT, GasOilParamsT, OilWaterParamsT, GasWwaterParamsT, false>& params) {

//     // Since the view should mainly work for stuff already initialized on the GPU we know that we already have raw
//     // pointers that do not need to be further dealt with, this make_view should just pass it on to fit the API
//     // using twoPhaseParamsWithShared = EclTwoPhaseMaterialParams<TraitsT, GasOilParamsT, OilWaterParamsT, GasWwaterParamsT, true>;
//     // typename twoPhaseParamsWithShared::GasOilParams* gasOilParams = params.gasOilParams().get();
//     // typename twoPhaseParamsWithShared::OilWaterParams* oilWaterParams = params.oilWaterParams().get();
//     // typename twoPhaseParamsWithShared::GasWaterParams* gasWaterParams = params.gasWaterParams().get();

//     auto resultView = EclTwoPhaseMaterialParams<TraitsT, GasOilParamsT, OilWaterParamsT, GasWwaterParamsT, false>();
//     resultView.setGasOilParams(params.gasOilParams());
//     resultView.setOilWaterParams(params.oilWaterParams());
//     resultView.setGasWaterParams(params.gasWaterParams());
//     resultView.setApproach(params.approach());
//     resultView.finalize();

//     return resultView;
// }
// }

#endif
