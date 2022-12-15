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
 * \copydoc Opm::GasPvtThermal
 */
#ifndef OPM_GAS_PVT_THERMAL_HPP
#define OPM_GAS_PVT_THERMAL_HPP

#include <opm/common/utility/Visitor.hpp>

#include <opm/material/common/Tabulated1DFunction.hpp>

#include <opm/material/fluidsystems/blackoilpvt/Co2GasPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/DryGasPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/DryHumidGasPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/PvtEnums.hpp>
#include <opm/material/fluidsystems/blackoilpvt/WetGasPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/WetHumidGasPvt.hpp>

#if HAVE_ECL_INPUT
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SimpleTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>
#endif

namespace Opm {

/*!
 * \brief This class implements temperature dependence of the PVT properties of gas
 *
 * Note that this _only_ implements the temperature part, i.e., it requires the
 * isothermal properties as input.
 */
template <class Scalar>
class GasPvtThermal
{
public:
    using IsothermalPvt = std::variant<std::monostate,
                                       DryGasPvt<Scalar>,
                                       DryHumidGasPvt<Scalar>,
                                       WetHumidGasPvt<Scalar>,
                                       WetGasPvt<Scalar>,
                                       Co2GasPvt<Scalar>>;
    using TabulatedOneDFunction = Tabulated1DFunction<Scalar>;

    static IsothermalPvt initialize(GasPvtApproach gasPvtAppr)
    {
        switch (gasPvtAppr) {
        case GasPvtApproach::DryGas:
            return DryGasPvt<Scalar>{};

        case GasPvtApproach::DryHumidGas:
            return DryHumidGasPvt<Scalar>{};

        case GasPvtApproach::WetHumidGas:
            return WetHumidGasPvt<Scalar>{};

        case GasPvtApproach::WetGas:
            return WetGasPvt<Scalar>{};

        case GasPvtApproach::Co2Gas:
            return Co2GasPvt<Scalar>{};

        default:
            return std::monostate{};
        }
    }

#if HAVE_ECL_INPUT
    static GasPvtApproach chooseApproach(const EclipseState& eclState)
    {
        if (eclState.runspec().co2Storage())
            return GasPvtApproach::Co2Gas;
        else if (!eclState.getTableManager().getPvtgwTables().empty() && !eclState.getTableManager().getPvtgTables().empty())
            return GasPvtApproach::WetHumidGas;
        else if (!eclState.getTableManager().getPvtgTables().empty())
            return GasPvtApproach::WetGas;
        else if (eclState.getTableManager().hasTables("PVDG"))
            return GasPvtApproach::DryGas;
        else if (!eclState.getTableManager().getPvtgwTables().empty())
            return GasPvtApproach::DryHumidGas;

        return GasPvtApproach::NoGas;
    }

    /*!
     * \brief Implement the temperature part of the gas PVT properties.
     */
    void initFromState(const EclipseState& eclState, const Schedule& schedule)
    {
        //////
        // initialize the isothermal part
        //////
        isothermalPvt_ = initialize(chooseApproach(eclState));
        std::visit(VisitorOverloadSet{[&](auto& pvt)
                                      {
                                          pvt.initFromState(eclState, schedule);
                                      }, monoHandler_}, isothermalPvt_);

        //////
        // initialize the thermal part
        //////
        const auto& tables = eclState.getTableManager();

        enableThermalDensity_ = tables.GasDenT().size() > 0;
        enableJouleThomson_ = tables.GasJT().size() > 0;
        enableThermalViscosity_ = tables.hasTables("GASVISCT");
        enableInternalEnergy_ = tables.hasTables("SPECHEAT");

        unsigned numRegions;
        std::visit(VisitorOverloadSet{[&](const auto& pvt)
                                      {
                                          numRegions = pvt.numRegions();
                                      }, monoHandler_}, isothermalPvt_);
        setNumRegions(numRegions);

        // viscosity
        if (enableThermalViscosity_) {
            const auto& gasvisctTables = tables.getGasvisctTables();
            auto gasCompIdx = tables.gas_comp_index();
            std::string gasvisctColumnName = "Viscosity" + std::to_string(static_cast<long long>(gasCompIdx));

            for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
                const auto& T = gasvisctTables[regionIdx].getColumn("Temperature").vectorCopy();
                const auto& mu = gasvisctTables[regionIdx].getColumn(gasvisctColumnName).vectorCopy();
                gasvisctCurves_[regionIdx].setXYContainers(T, mu);
            }
        }

        // temperature dependence of gas density
        if (enableThermalDensity_) {
            const auto& gasDenT = tables.GasDenT();

            assert(gasDenT.size() == numRegions);
            for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
                const auto& record = gasDenT[regionIdx];

                gasdentRefTemp_[regionIdx] = record.T0;
                gasdentCT1_[regionIdx] = record.C1;
                gasdentCT2_[regionIdx] = record.C2;
            }
        }

        // Joule Thomson 
        if (enableJouleThomson_) {
            const auto& gasJT = tables.GasJT();

            assert(gasJT.size() == numRegions);
            for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
                const auto& record = gasJT[regionIdx];

                gasJTRefPres_[regionIdx] =  record.P0;
                gasJTC_[regionIdx] = record.C1;
            }

            const auto& densityTable = eclState.getTableManager().getDensityTable();

            assert(densityTable.size() == numRegions);
            for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
                 rhoRefO_[regionIdx] = densityTable[regionIdx].oil;
            }
        }

        if (enableInternalEnergy_) {
            // the specific internal energy of gas. be aware that ecl only specifies the heat capacity
            // (via the SPECHEAT keyword) and we need to integrate it ourselfs to get the
            // internal energy
            for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
                const auto& specHeatTable = tables.getSpecheatTables()[regionIdx];
                const auto& temperatureColumn = specHeatTable.getColumn("TEMPERATURE");
                const auto& cvGasColumn = specHeatTable.getColumn("CV_GAS");

                std::vector<double> uSamples(temperatureColumn.size());

                // the specific enthalpy of vaporization. since ECL does not seem to
                // feature a proper way to specify this quantity, we use the value for
                // methane. A proper model would also need to consider the enthalpy
                // change due to dissolution, i.e. the enthalpies of the gas and oil
                // phases should depend on the phase composition
                const Scalar hVap = 480.6e3; // [J / kg]

                Scalar u = temperatureColumn[0]*cvGasColumn[0] + hVap;
                for (size_t i = 0;; ++i) {
                    uSamples[i] = u;

                    if (i >= temperatureColumn.size() - 1)
                        break;

                    // integrate to the heat capacity from the current sampling point to the next
                    // one. this leads to a quadratic polynomial.
                    Scalar c_v0 = cvGasColumn[i];
                    Scalar c_v1 = cvGasColumn[i + 1];
                    Scalar T0 = temperatureColumn[i];
                    Scalar T1 = temperatureColumn[i + 1];
                    u += 0.5*(c_v0 + c_v1)*(T1 - T0);
                }

                internalEnergyCurves_[regionIdx].setXYContainers(temperatureColumn.vectorCopy(), uSamples);
            }
        }
    }
#endif // HAVE_ECL_INPUT

    /*!
     * \brief Set the number of PVT-regions considered by this object.
     */
    void setNumRegions(size_t numRegions)
    {
        gasvisctCurves_.resize(numRegions);
        internalEnergyCurves_.resize(numRegions);
        gasdentRefTemp_.resize(numRegions);
        gasdentCT1_.resize(numRegions);
        gasdentCT2_.resize(numRegions);
        gasJTRefPres_.resize(numRegions);
        gasJTC_.resize(numRegions);
        rhoRefO_.resize(numRegions);
    }

    /*!
     * \brief Finish initializing the thermal part of the gas phase PVT properties.
     */
    void initEnd()
    { }

    size_t numRegions() const
    { return gasvisctCurves_.size(); }

    /*!
     * \brief Returns true iff the density of the gas phase is temperature dependent.
     */
    bool enableThermalDensity() const
    { return enableThermalDensity_; }

     /*!
     * \brief Returns true iff Joule-Thomson effect for the gas phase is active.
     */
    bool enableJouleThomsony() const
    { return enableJouleThomson_; }

    /*!
     * \brief Returns true iff the viscosity of the gas phase is temperature dependent.
     */
    bool enableThermalViscosity() const
    { return enableThermalViscosity_; }

    /*!
     * \brief Returns the specific internal energy [J/kg] of gas given a set of parameters.
     */
    template <class Evaluation>
    Evaluation internalEnergy(unsigned regionIdx,
                              const Evaluation& temperature,
                              const Evaluation& pressure,
                              const Evaluation& Rv) const
    {
        if (!enableInternalEnergy_)
             throw std::runtime_error("Requested the internal energy of gas but it is disabled");

        if (!enableJouleThomson_) {
            // compute the specific internal energy for the specified tempature. We use linear
            // interpolation here despite the fact that the underlying heat capacities are
            // piecewise linear (which leads to a quadratic function)
            return internalEnergyCurves_[regionIdx].eval(temperature, /*extrapolate=*/true);
        }
        else {
            Evaluation Tref = gasdentRefTemp_[regionIdx];
            Evaluation Pref = gasJTRefPres_[regionIdx]; 
            Scalar JTC = gasJTC_[regionIdx]; // if JTC is default then JTC is calculated
            Evaluation Rvw = 0.0;

            Evaluation invB = inverseFormationVolumeFactor(regionIdx, temperature, pressure, Rv, Rvw);
            constexpr const Scalar hVap = 480.6e3; // [J / kg]
            Evaluation Cp = (internalEnergyCurves_[regionIdx].eval(temperature, /*extrapolate=*/true) - hVap)/temperature;
            Evaluation density = invB * (gasReferenceDensity(regionIdx) + Rv * rhoRefO_[regionIdx]);

            Evaluation enthalpyPres;
            if  (JTC != 0) {
                enthalpyPres = -Cp * JTC * (pressure -Pref);
            }
            else if(enableThermalDensity_) {
                Scalar c1T = gasdentCT1_[regionIdx];
                Scalar c2T = gasdentCT2_[regionIdx];
              
                Evaluation alpha = (c1T + 2 * c2T * (temperature - Tref)) /
                    (1 + c1T  *(temperature - Tref) + c2T * (temperature - Tref) * (temperature - Tref));

                constexpr const int N = 100; // value is experimental
                Evaluation deltaP = (pressure - Pref)/N;
                Evaluation enthalpyPresPrev = 0;
                for (size_t i = 0; i < N; ++i) {
                    Evaluation Pnew = Pref + i * deltaP;
                    Evaluation rho = inverseFormationVolumeFactor(regionIdx, temperature, Pnew, Rv, Rvw) *
                                     (gasReferenceDensity(regionIdx) + Rv * rhoRefO_[regionIdx]);
                    // see e.g.https://en.wikipedia.org/wiki/Joule-Thomson_effect for a derivation of the Joule-Thomson coeff.
                    Evaluation jouleThomsonCoefficient = -(1.0/Cp) * (1.0 - alpha * temperature)/rho;
                    Evaluation deltaEnthalpyPres = -Cp * jouleThomsonCoefficient * deltaP;
                    enthalpyPres = enthalpyPresPrev + deltaEnthalpyPres; 
                    enthalpyPresPrev = enthalpyPres;
                }
            }
            else {
                  throw std::runtime_error("Requested Joule-thomson calculation but thermal gas density (GASDENT) is not provided");
            }

            Evaluation enthalpy = Cp * (temperature - Tref) + enthalpyPres;

            return enthalpy - pressure/density;
        }
    }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation viscosity(unsigned regionIdx,
                         const Evaluation& temperature,
                         const Evaluation& pressure,
                         const Evaluation& Rv,
                         const Evaluation& Rvw) const
    {
        if (!enableThermalViscosity()) {
            Evaluation result;
            std::visit(VisitorOverloadSet{[&](const auto& pvt)
                                          {
                                              result = pvt.viscosity(regionIdx,
                                                                     temperature,
                                                                     pressure, Rv, Rvw);
                                          }, monoHandler_}, isothermalPvt_);
            return result;
        }

        // compute the viscosity deviation due to temperature
        const auto& muGasvisct = gasvisctCurves_[regionIdx].eval(temperature);
        return muGasvisct;
    }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the oil-saturated gas phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation saturatedViscosity(unsigned regionIdx,
                                  const Evaluation& temperature,
                                  const Evaluation& pressure) const
    {
        if (!enableThermalViscosity()) {
            Evaluation result;
            std::visit(VisitorOverloadSet{[&](const auto& pvt)
                                          {
                                              result = pvt.saturatedViscosity(regionIdx,
                                                                              temperature,
                                                                              pressure);
                                          }, monoHandler_}, isothermalPvt_);
            return result;
        }

        // compute the viscosity deviation due to temperature
        const auto& muGasvisct = gasvisctCurves_[regionIdx].eval(temperature, true);
        return muGasvisct;
    }

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation>
    Evaluation inverseFormationVolumeFactor(unsigned regionIdx,
                                            const Evaluation& temperature,
                                            const Evaluation& pressure,
                                            const Evaluation& Rv,
                                            const Evaluation& /*Rvw*/) const
    {
        Evaluation b;
        std::visit(VisitorOverloadSet{[&](const auto& pvt)
                                         {
                                             b = pvt.inverseFormationVolumeFactor(regionIdx,
                                                                                  temperature,
                                                                                  pressure,
                                                                                  Rv,
                                                                                  Evaluation{0.0});
                                          }, monoHandler_}, isothermalPvt_);

        if (!enableThermalDensity())
            return b;

        // we use the same approach as for the for water here, but with the OPM-specific
        // GASDENT keyword.
        //
        // TODO: Since gas is quite a bit more compressible than water, it might be
        //       necessary to make GASDENT to a table keyword. If the current temperature
        //       is relatively close to the reference temperature, the current approach
        //       should be good enough, though.
        Scalar TRef = gasdentRefTemp_[regionIdx];
        Scalar cT1 = gasdentCT1_[regionIdx];
        Scalar cT2 = gasdentCT2_[regionIdx];
        const Evaluation& Y = temperature - TRef;

        return b/(1 + (cT1 + cT2*Y)*Y);
    }

    /*!
     * \brief Returns the formation volume factor [-] of oil-saturated gas.
     */
    template <class Evaluation>
    Evaluation saturatedInverseFormationVolumeFactor(unsigned regionIdx,
                                                     const Evaluation& temperature,
                                                     const Evaluation& pressure) const
    {
        Evaluation b;
        std::visit(VisitorOverloadSet{[&](const auto& pvt)
                                         {
                                             b = pvt.saturatedInverseFormationVolumeFactor(regionIdx,
                                                                                          temperature,
                                                                                          pressure);
                                          }, monoHandler_}, isothermalPvt_);

        if (!enableThermalDensity())
            return b;

        // we use the same approach as for the for water here, but with the OPM-specific
        // GASDENT keyword.
        //
        // TODO: Since gas is quite a bit more compressible than water, it might be
        //       necessary to make GASDENT to a table keyword. If the current temperature
        //       is relatively close to the reference temperature, the current approach
        //       should be good enough, though.
        Scalar TRef = gasdentRefTemp_[regionIdx];
        Scalar cT1 = gasdentCT1_[regionIdx];
        Scalar cT2 = gasdentCT2_[regionIdx];
        const Evaluation& Y = temperature - TRef;

        return b/(1 + (cT1 + cT2*Y)*Y);
    }

    /*!
     * \brief Returns the water vaporization factor \f$R_v\f$ [m^3/m^3] of the water phase.
     */
    template <class Evaluation>
    Evaluation saturatedWaterVaporizationFactor(unsigned /*regionIdx*/,
                                              const Evaluation& /*temperature*/,
                                              const Evaluation& /*pressure*/) const
    { return 0.0; }

     /*!
     * \brief Returns the water vaporization factor \f$R_vw\f$ [m^3/m^3] of water saturated gas.
     */
    template <class Evaluation = Scalar>
    Evaluation saturatedWaterVaporizationFactor(unsigned /*regionIdx*/,
                                              const Evaluation& /*temperature*/,
                                              const Evaluation& /*pressure*/, 
                                              const Evaluation& /*saltConcentration*/) const
    { return 0.0; }

    

    /*!
     * \brief Returns the oil vaporization factor \f$R_v\f$ [m^3/m^3] of the gas phase.
     *
     * This method implements temperature dependence and requires the gas pressure,
     * temperature and the oil saturation as inputs. Currently it is just a dummy method
     * which passes through the isothermal oil vaporization factor.
     */
    template <class Evaluation>
    Evaluation saturatedOilVaporizationFactor(unsigned regionIdx,
                                              const Evaluation& temperature,
                                              const Evaluation& pressure) const
    {
        Evaluation result;
        std::visit(VisitorOverloadSet{[&](const auto& pvt)
                                         {
                                             result = pvt.saturatedOilVaporizationFactor(regionIdx,
                                                                                         temperature,
                                                                                         pressure);
                                          }, monoHandler_}, isothermalPvt_);

        return result;
    }

    /*!
     * \brief Returns the oil vaporization factor \f$R_v\f$ [m^3/m^3] of the gas phase.
     *
     * This method implements temperature dependence and requires the gas pressure,
     * temperature and the oil saturation as inputs. Currently it is just a dummy method
     * which passes through the isothermal oil vaporization factor.
     */
    template <class Evaluation>
    Evaluation saturatedOilVaporizationFactor(unsigned regionIdx,
                                              const Evaluation& temperature,
                                              const Evaluation& pressure,
                                              const Evaluation& oilSaturation,
                                              const Evaluation& maxOilSaturation) const
    {
        Evaluation result;
        std::visit(VisitorOverloadSet{[&](const auto& pvt)
                                         {
                                             result = pvt.saturatedOilVaporizationFactor(regionIdx,
                                                                                         temperature,
                                                                                         pressure,
                                                                                         oilSaturation,
                                                                                         maxOilSaturation);
                                          }, monoHandler_}, isothermalPvt_);

        return result;
    }

    /*!
     * \brief Returns the saturation pressure of the gas phase [Pa]
     *
     * This method implements temperature dependence and requires isothermal satuation
     * pressure and temperature as inputs. Currently it is just a dummy method which
     * passes through the isothermal saturation pressure.
     */
    template <class Evaluation>
    Evaluation saturationPressure(unsigned regionIdx,
                                  const Evaluation& temperature,
                                  const Evaluation& pressure) const
    {
        Evaluation result;
        std::visit(VisitorOverloadSet{[&](const auto& pvt)
                                         {
                                             result = pvt.saturationPressure(regionIdx,
                                                                             temperature,
                                                                             pressure);
                                          }, monoHandler_}, isothermalPvt_);

        return result;
    }

    template <class Evaluation>
    Evaluation diffusionCoefficient(const Evaluation& temperature,
                                    const Evaluation& pressure,
                                    unsigned compIdx) const
    {
        Evaluation result;
        std::visit(VisitorOverloadSet{[&](const auto& pvt)
                                         {
                                             result = pvt.diffusionCoefficient(temperature,
                                                                               pressure,
                                                                               compIdx);
                                          }, monoHandler_}, isothermalPvt_);

        return result;
    }

    Scalar gasReferenceDensity(unsigned regionIdx) const
    {
        Scalar result;
        std::visit(VisitorOverloadSet{[&](const auto& pvt)
                                         {
                                             result = pvt.gasReferenceDensity(regionIdx);
                                          }, monoHandler_}, isothermalPvt_);

        return result;
    }

    const std::vector<TabulatedOneDFunction>& gasvisctCurves() const
    { return gasvisctCurves_; }

    const std::vector<Scalar>& gasdentRefTemp() const
    { return gasdentRefTemp_; }

    const std::vector<Scalar>& gasdentCT1() const
    { return gasdentCT1_; }

    const std::vector<Scalar>& gasdentCT2() const
    { return gasdentCT2_; }

    const std::vector<TabulatedOneDFunction>& internalEnergyCurves() const
    { return internalEnergyCurves_; }

    bool enableInternalEnergy() const
    { return enableInternalEnergy_; }

    const std::vector<Scalar>& gasJTRefPres() const
    { return  gasJTRefPres_; }

     const std::vector<Scalar>&  gasJTC() const
    { return gasJTC_; }

private:
    IsothermalPvt isothermalPvt_;

    // The PVT properties needed for temperature dependence of the viscosity. We need
    // to store one value per PVT region.
    std::vector<TabulatedOneDFunction> gasvisctCurves_;

    std::vector<Scalar> gasdentRefTemp_;
    std::vector<Scalar> gasdentCT1_;
    std::vector<Scalar> gasdentCT2_;

    std::vector<Scalar> gasJTRefPres_;
    std::vector<Scalar> gasJTC_;

    std::vector<Scalar> rhoRefO_;

    // piecewise linear curve representing the internal energy of gas
    std::vector<TabulatedOneDFunction> internalEnergyCurves_;

    bool enableThermalDensity_ = false;
    bool enableJouleThomson_ = false;
    bool enableThermalViscosity_ = false;
    bool enableInternalEnergy_ = false;

    MonoThrowHandler<std::logic_error>
    monoHandler_{"Not implemented: Gas PVT of this deck!"}; // mono state handler
};

} // namespace Opm

#endif
