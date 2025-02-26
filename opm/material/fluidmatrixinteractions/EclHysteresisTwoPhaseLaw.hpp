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
 * \copydoc Opm::EclHysteresisTwoPhaseLaw
 */
#ifndef OPM_ECL_HYSTERESIS_TWO_PHASE_LAW_HPP
#define OPM_ECL_HYSTERESIS_TWO_PHASE_LAW_HPP
#include <opm/common/TimingMacros.hpp>
#include "EclHysteresisTwoPhaseLawParams.hpp"
#include <stdexcept>

namespace Opm {

/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief This material law implements the hysteresis model of the ECL file format
 */
template <class EffectiveLawT,
          class ParamsT = EclHysteresisTwoPhaseLawParams<EffectiveLawT> >
class EclHysteresisTwoPhaseLaw : public EffectiveLawT::Traits
{
public:
    using EffectiveLaw = EffectiveLawT;
    using EffectiveLawParams = typename EffectiveLaw::Params;

    using Traits = typename EffectiveLaw::Traits;
    using Params = ParamsT;
    using Scalar = typename EffectiveLaw::Scalar;

    enum { wettingPhaseIdx = Traits::wettingPhaseIdx };
    enum { nonWettingPhaseIdx = Traits::nonWettingPhaseIdx };

    //! The number of fluid phases
    static constexpr int numPhases = EffectiveLaw::numPhases;
    static_assert(numPhases == 2,
                  "The endpoint scaling applies to the nested twophase laws, not to "
                  "the threephase one!");

    //! Specify whether this material law implements the two-phase
    //! convenience API
    static constexpr bool implementsTwoPhaseApi = true;

    static_assert(EffectiveLaw::implementsTwoPhaseApi,
                  "The material laws put into EclEpsTwoPhaseLaw must implement the "
                  "two-phase material law API!");

    //! Specify whether this material law implements the two-phase
    //! convenience API which only depends on the phase saturations
    static constexpr bool implementsTwoPhaseSatApi = true;

    static_assert(EffectiveLaw::implementsTwoPhaseSatApi,
                  "The material laws put into EclEpsTwoPhaseLaw must implement the "
                  "two-phase material law saturation API!");

    //! Specify whether the quantities defined by this material law
    //! are saturation dependent
    static constexpr bool isSaturationDependent = true;

    //! Specify whether the quantities defined by this material law
    //! are dependent on the absolute pressure
    static constexpr bool isPressureDependent = false;

    //! Specify whether the quantities defined by this material law
    //! are temperature dependent
    static constexpr bool isTemperatureDependent = false;

    //! Specify whether the quantities defined by this material law
    //! are dependent on the phase composition
    static constexpr bool isCompositionDependent = false;

    static constexpr bool isHysteresisDependent = true;

    /*!
     * \brief The capillary pressure-saturation curves depending on absolute saturations.
     *
     * \param values A random access container which stores the
     *               relative pressure of each fluid phase.
     * \param params The parameter object expressing the coefficients
     *               required by the van Genuchten law.
     * \param fs The fluid state for which the capillary pressure
     *           ought to be calculated
     */
    template <class Container, class FluidState>
    static void capillaryPressures(Container& /* values */,
                                   const Params& /* params */,
                                   const FluidState& /* fs */)
    {
        throw std::invalid_argument("The capillaryPressures(fs) method is not yet implemented");
    }

    /*!
     * \brief The relative permeability-saturation curves depending on absolute saturations.
     *
     * \param values A random access container which stores the
     *               relative permeability of each fluid phase.
     * \param params The parameter object expressing the coefficients
     *               required by the van Genuchten law.
     * \param fs The fluid state for which the relative permeabilities
     *           ought to be calculated
     */
    template <class Container, class FluidState>
    static void relativePermeabilities(Container& /* values */,
                                       const Params& /* params */,
                                       const FluidState& /* fs */)
    {
        throw std::invalid_argument("The pcnw(fs) method is not yet implemented");
    }

    /*!
     * \brief The capillary pressure-saturation curve.
     *
     *
     * \param params A object that stores the appropriate coefficients
     *                for the respective law.
     *
     * \return Capillary pressure [Pa] calculated by specific
     *         constitutive relation (e.g. Brooks & Corey, van
     *         Genuchten, linear...)
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation pcnw(const Params& /* params */,
                           const FluidState& /* fs */)
    {
        throw std::invalid_argument("The pcnw(fs) method is not yet implemented");
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatPcnw(const Params& params, const Evaluation& Sw)
    {
        OPM_TIMEFUNCTION_LOCAL();
        // if no pc hysteresis is enabled, use the drainage curve
        if (!params.config().enableHysteresis() || params.config().pcHysteresisModel() < 0)
            return EffectiveLaw::twoPhaseSatPcnw(params.drainageParams(), Sw);

        // Initial imbibition process
        if (params.initialImb()) {
            if (Sw >= params.pcSwMic()) {
                return EffectiveLaw::twoPhaseSatPcnw(params.imbibitionParams(), Sw);
            }
            else { // Reversal
                const Evaluation& F = (1.0/(params.pcSwMic()-Sw+params.curvatureCapPrs())-1.0/params.curvatureCapPrs())
                     / (1.0/(params.pcSwMic()-params.Swcrd()+params.curvatureCapPrs())-1.0/params.curvatureCapPrs());

                const Evaluation& Pcd = EffectiveLaw::twoPhaseSatPcnw(params.drainageParams(), Sw);
                const Evaluation& Pci = EffectiveLaw::twoPhaseSatPcnw(params.imbibitionParams(), Sw);
                const Evaluation& pc_Killough = Pci+F*(Pcd-Pci);

                return pc_Killough;
            }
        }

        // Initial drainage process
        if (Sw <= params.pcSwMdc())
            return EffectiveLaw::twoPhaseSatPcnw(params.drainageParams(), Sw);

        // Reversal
        Scalar Swma = 1.0-params.Sncrt();
        if (Sw >= Swma) {
            const Evaluation& Pci = EffectiveLaw::twoPhaseSatPcnw(params.imbibitionParams(), Sw);
            return Pci;
        }
        else {
            Scalar pciwght = params.pcWght(); // Align pci and pcd at Swir
            //const Evaluation& SwEff = params.Swcri()+(Sw-params.pcSwMdc())/(Swma-params.pcSwMdc())*(Swma-params.Swcri());
            const Evaluation& SwEff = Sw; // This is Killough 1976, Gives significantly better fit compared to benchmark then the above "scaling"
            const Evaluation& Pci = pciwght*EffectiveLaw::twoPhaseSatPcnw(params.imbibitionParams(), SwEff);

            const Evaluation& Pcd = EffectiveLaw::twoPhaseSatPcnw(params.drainageParams(), Sw);

            if (Pci == Pcd)
                return Pcd;

            const Evaluation& F = (1.0/(Sw-params.pcSwMdc()+params.curvatureCapPrs())-1.0/params.curvatureCapPrs())
                                / (1.0/(Swma-params.pcSwMdc()+params.curvatureCapPrs())-1.0/params.curvatureCapPrs());

            const Evaluation& pc_Killough = Pcd+F*(Pci-Pcd);

            return pc_Killough;
        }

        return 0.0;
    }

    /*!
     * \brief The saturation-capillary pressure curves.
     */
    template <class Container, class FluidState>
    static void saturations(Container& /* values */,
                            const Params& /* params */,
                            const FluidState& /* fs */)
    {
        throw std::invalid_argument("The saturations(fs) method is not yet implemented");
    }

    /*!
     * \brief Calculate wetting liquid phase saturation given that
     *        the rest of the fluid state has been initialized
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation Sw(const Params& /* params */,
                         const FluidState& /* fs */)
    {
        throw std::invalid_argument("The Sw(fs) method is not yet implemented");
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatSw(const Params& /* params */,
                                    const Evaluation& /* pc */)
    {
        throw std::invalid_argument("The twoPhaseSatSw(pc) method is not yet implemented");
    }

    /*!
     * \brief Calculate non-wetting liquid phase saturation given that
     *        the rest of the fluid state has been initialized
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation Sn(const Params& /* params */,
                         const FluidState& /* fs */)
    {
        throw std::invalid_argument("The Sn(pc) method is not yet implemented");
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatSn(const Params& /* params */,
                                    const Evaluation& /* pc */)
    {
        throw std::invalid_argument("The twoPhaseSatSn(pc) method is not yet implemented");
    }

    /*!
     * \brief The relative permeability for the wetting phase.
     *
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Relative permeability of the wetting phase calculated as implied by EffectiveLaw e.g. Brooks & Corey, van Genuchten, linear... .
     *
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation krw(const Params& /* params */,
                          const FluidState& /* fs */)
    {
        throw std::invalid_argument("The krw(fs) method is not yet implemented");
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatKrw(const Params& params, const Evaluation& Sw)
    {

        OPM_TIMEFUNCTION_LOCAL();
        // if no relperm hysteresis is enabled, use the drainage curve
        if (!params.config().enableHysteresis() || params.config().krHysteresisModel() < 0)
            return EffectiveLaw::twoPhaseSatKrw(params.drainageParams(), Sw);

        if (params.config().krHysteresisModel() == 0 || params.config().krHysteresisModel() == 2)
            // use drainage curve for wetting phase
            return EffectiveLaw::twoPhaseSatKrw(params.drainageParams(), Sw);

        // use imbibition curve for wetting phase
        if (params.config().krHysteresisModel() == 1 || params.config().krHysteresisModel() == 3)
            return EffectiveLaw::twoPhaseSatKrw(params.imbibitionParams(), Sw);

        if (Sw <= params.krnSwMdc()) {
            return EffectiveLaw::twoPhaseSatKrw(params.drainageParams(), Sw);
        }
        // Killough hysteresis for the wetting phase
        assert(params.config().krHysteresisModel() == 4);
        Evaluation Snorm = params.Sncri()+(1.0-Sw-params.Sncrt())*(params.Snmaxd()-params.Sncri())/(params.Snhy()-params.Sncrt());
        Evaluation Krwi_snorm = EffectiveLaw::twoPhaseSatKrw(params.imbibitionParams(), 1 - Snorm);
        return params.KrwdHy() +  params.krwWght() * (Krwi_snorm - params.Krwi_snmax());
    }

    /*!
     * \brief The relative permeability of the non-wetting phase.
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation krn(const Params& /* params */,
                          const FluidState& /* fs */)
    {
        throw std::invalid_argument("The krn(fs) method is not yet implemented");
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatKrn(const Params& params, const Evaluation& Sw)
    {
        OPM_TIMEFUNCTION_LOCAL();
        // If WAG hysteresis is enabled, the convential hysteresis model is ignored.
        // (Two-phase model, non-wetting: only gas in oil.)
        if (params.gasOilHysteresisWAG()) {

            // Primary drainage
            if (Sw <= params.krnSwMdc() + params.tolWAG() && params.nState() == 1) {
                return EffectiveLaw::twoPhaseSatKrn(params.drainageParams(), Sw);
            }

            // Imbibition or reversion to two-phase drainage retracing imb curve
            // (Shift along primary drainage curve.)
            if (params.nState() == 1) {
                Evaluation Swf = params.computeSwf(Sw);
                return EffectiveLaw::twoPhaseSatKrn(params.drainageParams(), Swf);
            }

            // Three-phase drainage along current secondary drainage curve
            if (Sw <= params.krnSwDrainRevert()+params.tolWAG() /*&& params.nState()>=1 */) {
                Evaluation Krg = EffectiveLaw::twoPhaseSatKrn(params.drainageParams(), Sw);
                Evaluation KrgDrain2 = (Krg-params.krnDrainStart())*params.reductionDrain() + params.krnImbStart();
                return KrgDrain2;
            }

            // Subsequent imbibition: Scanning curve derived from previous secondary drainage
            if (Sw >= params.krnSwWAG()-params.tolWAG() /*&& Sw > params.krnSwDrainRevert() && params.nState()>=1 */) {
                Evaluation KrgImb2 = params.computeKrImbWAG(Sw);
                return KrgImb2;
            }
            else {/* Sw < params.krnSwWAG() */  // Reversion along "next" drainage curve
                Evaluation Krg = EffectiveLaw::twoPhaseSatKrn(params.drainageParams(), Sw);
                Evaluation KrgDrainNxt = (Krg-params.krnDrainStartNxt())*params.reductionDrainNxt() + params.krnImbStartNxt();
                return KrgDrainNxt;
            }
        }
        // if no relperm hysteresis is enabled, use the drainage curve
        if (!params.config().enableHysteresis() || params.config().krHysteresisModel() < 0)
            return EffectiveLaw::twoPhaseSatKrn(params.drainageParams(), Sw);

        // if it is enabled, use either the drainage or the imbibition curve. if the
        // imbibition curve is used, the saturation must be shifted.
        if (Sw <= params.krnSwMdc()) {
            return EffectiveLaw::twoPhaseSatKrn(params.drainageParams(), Sw);
        }

        if (params.config().krHysteresisModel() <= 1) { //Carlson
            return EffectiveLaw::twoPhaseSatKrn(params.imbibitionParams(),
                                                Sw + params.deltaSwImbKrn());
        }

        // Killough
        assert(params.config().krHysteresisModel() == 2 || params.config().krHysteresisModel() == 3 || params.config().krHysteresisModel() == 4);
        Evaluation Snorm = params.Sncri()+(1.0-Sw-params.Sncrt())*(params.Snmaxd()-params.Sncri())/(params.Snhy()-params.Sncrt());
        return params.krnWght()*EffectiveLaw::twoPhaseSatKrn(params.imbibitionParams(),1.0-Snorm);
    }
};

} // namespace Opm

#endif
