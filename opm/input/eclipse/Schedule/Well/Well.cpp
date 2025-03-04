/*
  Copyright 2019 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <opm/input/eclipse/Schedule/Well/Well.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/io/eclipse/rst/well.hpp>
#include <opm/output/eclipse/VectorItems/well.hpp>

#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/EclipseState/Phase.hpp>
#include <opm/input/eclipse/EclipseState/TracerConfig.hpp>

#include <opm/input/eclipse/Schedule/MSW/WellSegments.hpp>
#include <opm/input/eclipse/Schedule/ScheduleGrid.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQActive.hpp>
#include <opm/input/eclipse/Schedule/Well/WDFAC.hpp>
#include <opm/input/eclipse/Schedule/Well/WellBrineProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>
#include <opm/input/eclipse/Schedule/Well/WellEconProductionLimits.hpp>
#include <opm/input/eclipse/Schedule/Well/WellFoamProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/WellMICPProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/WellPolymerProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTracerProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/WVFPDP.hpp>
#include <opm/input/eclipse/Schedule/Well/WVFPEXP.hpp>

#include <opm/input/eclipse/Units/Units.hpp>

#include <opm/common/utility/OpmInputError.hpp>
#include <opm/common/utility/shmatch.hpp>

#include <opm/input/eclipse/Parser/ParserKeywords/C.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/S.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/W.hpp>

#include <opm/input/eclipse/Deck/DeckRecord.hpp>
#include <opm/input/eclipse/Deck/DeckKeyword.hpp>

#include "../MSW/Compsegs.hpp"
#include "../eval_uda.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <fmt/format.h>

namespace {

    bool defaulted(const Opm::DeckRecord& rec, const std::string& s)
    {
        const auto& item = rec.getItem(s);
        return item.defaultApplied(0) || (item.get<int>(0) == 0);
    }

    int limit(const Opm::DeckRecord& rec, const std::string& s, const int shift)
    {
        const auto& item = rec.getItem(s);
        return shift + item.get<int>(0);
    }

    bool match_le(const int              value,
                  const Opm::DeckRecord& rec,
                  const std::string&     s,
                  const int              shift = 0)
    {
        return defaulted(rec, s) || (value <= limit(rec, s, shift));
    }

    bool match_ge(const int              value,
                  const Opm::DeckRecord& rec,
                  const std::string&     s,
                  const int              shift = 0)
    {
        return defaulted(rec, s) || (value >= limit(rec, s, shift));
    }

    bool match_eq(const int              value,
                  const Opm::DeckRecord& rec,
                  const std::string&     s,
                  const int              shift = 0)
    {
        return defaulted(rec, s)|| (limit(rec, s, shift) == value);
    }

Opm::Connection::Order order_from_int(const int int_value)
{
    switch (int_value) {
    case 0: return Opm::Connection::Order::TRACK;
    case 1: return Opm::Connection::Order::DEPTH;
    case 2: return Opm::Connection::Order::INPUT;
    default:
        throw std::invalid_argument {
            fmt::format("Invalid integer value: {} encountered "
                        "when determining connection ordering", int_value)
        };
    }
}

Opm::Well::Status status_from_int(const int int_value)
{
    using Value = Opm::RestartIO::Helpers::VectorItems::IWell::Value::Status;

    switch (int_value) {
    case Value::Shut: return Opm::Well::Status::SHUT;
    case Value::Stop: return Opm::Well::Status::STOP;
    case Value::Open: return Opm::Well::Status::OPEN;
    case Value::Auto: return Opm::Well::Status::AUTO;
    default:
        throw std::logic_error {
            fmt::format("integer value: {} could not be "
                        "converted to a valid well status.", int_value)
        };
    }
}

Opm::Well::ProducerCMode producer_cmode_from_int(const int pmode)
{
    using CModeVal = ::Opm::RestartIO::Helpers::VectorItems::
        IWell::Value::WellCtrlMode;

    switch (pmode) {
    case CModeVal::Group:    return Opm::Well::ProducerCMode::GRUP;
    case CModeVal::OilRate:  return Opm::Well::ProducerCMode::ORAT;
    case CModeVal::WatRate:  return Opm::Well::ProducerCMode::WRAT;
    case CModeVal::GasRate:  return Opm::Well::ProducerCMode::GRAT;
    case CModeVal::LiqRate:  return Opm::Well::ProducerCMode::LRAT;
    case CModeVal::ResVRate: return Opm::Well::ProducerCMode::RESV;
    case CModeVal::THP:      return Opm::Well::ProducerCMode::THP;
    case CModeVal::BHP:      return Opm::Well::ProducerCMode::BHP;
    }

    throw std::invalid_argument {
        fmt::format("Cannot convert integer value {} to producer control mode", pmode)
    };
}

Opm::Well::InjectorCMode injector_cmode_from_int(const int imode)
{
    using CModeVal = ::Opm::RestartIO::Helpers::VectorItems::
        IWell::Value::WellCtrlMode;

    switch (imode) {
    case CModeVal::Group:    return Opm::Well::InjectorCMode::GRUP;
    case CModeVal::OilRate:  [[fallthrough]];
    case CModeVal::WatRate:  [[fallthrough]];
    case CModeVal::GasRate:  [[fallthrough]];
    case CModeVal::LiqRate:  return Opm::Well::InjectorCMode::RATE;
    case CModeVal::ResVRate: return Opm::Well::InjectorCMode::RESV;
    case CModeVal::THP:      return Opm::Well::InjectorCMode::THP;
    case CModeVal::BHP:      return Opm::Well::InjectorCMode::BHP;
    }

    throw std::invalid_argument {
        fmt::format("Cannot convert integer value {} to injector control mode", imode)
    };
}

bool haveEconomicLimits(const Opm::RestartIO::RstWell& rst_well)
{
    namespace Limits = Opm::RestartIO::Helpers::VectorItems::
        IWell::Value::EconLimit;

    auto is_finite = [](const float x)
    {
        const auto infty = 1.0e+20f;

        return std::abs(x) < infty;
    };

    auto is_nonzero = [](const float x)
    {
        return std::abs(x) > 0.0f;
    };

    return (rst_well.econ_workover_procedure   != Limits::WOProcedure::None)
        || (rst_well.econ_workover_procedure_2 != Limits::WOProcedure::None)
        || (rst_well.econ_limit_end_run        == Limits::EndRun::Yes)
        || (rst_well.econ_limit_quantity       != Limits::Rate)
        || is_nonzero(rst_well.econ_limit_min_oil)
        || is_nonzero(rst_well.econ_limit_min_gas)
        || is_nonzero(rst_well.econ_limit_min_liq)
        || is_finite (rst_well.econ_limit_max_wct)
        || is_finite (rst_well.econ_limit_max_gor)
        || is_finite (rst_well.econ_limit_max_wgr)
        || is_finite (rst_well.econ_limit_max_wct_2);
}

std::shared_ptr<Opm::WellEconProductionLimits>
economicLimits(const Opm::RestartIO::RstWell& rst_well)
{
    return haveEconomicLimits(rst_well)
        ? std::make_shared<Opm::WellEconProductionLimits>(rst_well)
        : std::make_shared<Opm::WellEconProductionLimits>();
}

Opm::Well::GuideRateTarget guideRatePhase(const int gr_phase)
{
    namespace WGrupCon = Opm::RestartIO::Helpers::VectorItems::
        IWell::Value::WGrupCon;

    using Target = Opm::Well::GuideRateTarget;

    switch (gr_phase) {
    case WGrupCon::GRPhase::Defaulted:            return Target::UNDEFINED;
    case WGrupCon::GRPhase::Oil:                  return Target::OIL;
    case WGrupCon::GRPhase::Water:                return Target::WAT;
    case WGrupCon::GRPhase::Gas:                  return Target::GAS;
    case WGrupCon::GRPhase::Liquid:               return Target::LIQ;
    case WGrupCon::GRPhase::SurfaceInjectionRate: return Target::RAT;
    case WGrupCon::GRPhase::ReservoirVolumeRate:  return Target::RES;
    }

    throw std::invalid_argument {
        fmt::format("Cannot convert integer value {} "
                    "to guiderate phase target", gr_phase)
    };
}

bool isGroupControllable(const int gr_controllable_flag)
{
    return gr_controllable_flag != Opm::RestartIO::Helpers::VectorItems::
        IWell::Value::WGrupCon::Controllable::No;
}

double guideRateValue(const float gr_value)
{
    return (std::abs(gr_value) < 1.0e+20f)
        ? gr_value
        : Opm::ParserKeywords::WGRUPCON::GUIDE_RATE::defaultValue;
}

Opm::Well::WellGuideRate guideRate(const Opm::RestartIO::RstWell& rst_well)
{
    return {
        isGroupControllable(rst_well.group_controllable_flag),
        guideRateValue(rst_well.grupcon_gr_value),
        guideRatePhase(rst_well.grupcon_gr_phase),
        rst_well.grupcon_gr_scaling
    };
}

std::shared_ptr<Opm::WVFPEXP>
explicitTHPOptions(const Opm::RestartIO::RstWell& rst_well)
{
    auto options = std::make_shared<Opm::WVFPEXP>();
    options->update(rst_well);

    return options;
}

const static bool def_automatic_shutin = true;
constexpr double def_solvent_fraction = 0;

}

namespace Opm {

Well::Well(const RestartIO::RstWell& rst_well,
           const int report_step,
           const int rst_whistctl_cmode,
           const TracerConfig& tracer_config,
           const UnitSystem& unit_system_arg,
           const double udq_undefined_arg,
           const std::optional<VFPProdTable::ALQ_TYPE>& alq_type) :
    wname(rst_well.name),
    group_name(rst_well.group),
    init_step(report_step),
    headI(rst_well.ij[0]),
    headJ(rst_well.ij[1]),
    ref_depth((std::abs(rst_well.datum_depth) < 1.0e+20) ? std::optional<double>{ rst_well.datum_depth } : std::nullopt),
    drainage_radius(rst_well.drainage_radius),
    allow_cross_flow(rst_well.allow_xflow == 1),
    automatic_shutin(def_automatic_shutin),
    pvt_table(rst_well.pvt_table),
    unit_system(&unit_system_arg),
    udq_undefined(udq_undefined_arg),
    wtype(rst_well.wtype),
    guide_rate(guideRate(rst_well)),
    efficiency_factor(rst_well.efficiency_factor),
    use_efficiency_in_network(true),   // @TODO@ Find and read the actual value from restart
    solvent_fraction(def_solvent_fraction),
    has_produced(rst_well.void_total != 0),
    has_injected(rst_well.void_inj_total != 0),    
    prediction_mode(rst_well.hist_requested_control == 0),
    econ_limits(economicLimits(rst_well)),
    foam_properties(std::make_shared<WellFoamProperties>()),
    polymer_properties(std::make_shared<WellPolymerProperties>()),
    micp_properties(std::make_shared<WellMICPProperties>()),
    brine_properties(std::make_shared<WellBrineProperties>()),
    tracer_properties(std::make_shared<WellTracerProperties>()),
    connections(std::make_shared<WellConnections>(order_from_int(rst_well.completion_ordering), headI, headJ)),
    production(std::make_shared<WellProductionProperties>(unit_system_arg, wname)),
    injection(std::make_shared<WellInjectionProperties>(unit_system_arg, wname)),
    wvfpdp(std::make_shared<WVFPDP>()),
    wvfpexp(explicitTHPOptions(rst_well)),
    wdfac(std::make_shared<WDFAC>(rst_well)),
    status(status_from_int(rst_well.well_status)),
    well_inj_temperature(std::nullopt),
    well_inj_mult(std::nullopt)
{

    auto is_defined = [](const double value) {
        return value < (Opm::RestartIO::RstWell::UNDEFINED_VALUE / 2);
    };

    auto update_if_defined = [&](Opm::UDAValue& target, const double value) {
        if (is_defined(value))
            target.update(value);
    };

    if (this->wtype.producer()) {
        auto p = std::make_shared<WellProductionProperties>(*this->unit_system, wname);

        p->whistctl_cmode = (rst_whistctl_cmode > 0)
            ? producer_cmode_from_int(rst_whistctl_cmode)
            : Well::ProducerCMode::CMODE_UNDEFINED;

        update_if_defined(p->BHPTarget, rst_well.bhp_target_float);
        update_if_defined(p->OilRate, rst_well.orat_target);
        update_if_defined(p->WaterRate, rst_well.wrat_target);
        update_if_defined(p->GasRate, rst_well.grat_target);
        update_if_defined(p->LiquidRate, rst_well.lrat_target);
        update_if_defined(p->ResVRate, rst_well.resv_target);
        update_if_defined(p->ALQValue, rst_well.alq_value);

        p->VFPTableNumber = rst_well.vfp_table;

        if (alq_type) {
            p->ALQValue.set_dim(VFPProdTable::ALQDimension(*alq_type, unit_system_arg));
        } else {
            if (rst_well.vfp_table != 0) {
                // This may be OK for restart wells that never re-open, so warning only.
                OpmLog::warning(fmt::format("Well {}: Trying to restore the ALQ dimension for VFP table {}, which has not been loaded.", wname, rst_well.vfp_table));
            }
        }

        p->predictionMode = this->prediction_mode;

        if (is_defined(rst_well.orat_target))
            p->addProductionControl( Well::ProducerCMode::ORAT );

        if (is_defined(rst_well.wrat_target))
            p->addProductionControl( Well::ProducerCMode::WRAT );

        if (is_defined(rst_well.grat_target))
            p->addProductionControl( Well::ProducerCMode::GRAT );

        if (is_defined(rst_well.lrat_target))
            p->addProductionControl( Well::ProducerCMode::LRAT );

        if (is_defined(rst_well.resv_target))
            p->addProductionControl( Well::ProducerCMode::RESV );

        if (rst_well.thp_target != 0.0f) {  // THP has defined default of 0.0
            p->THPTarget.update(rst_well.thp_target);
            p->addProductionControl( Well::ProducerCMode::THP );
        }

        if (! p->predictionMode)
            p->clearControls();

        p->controlMode = producer_cmode_from_int(rst_well.active_control);
        p->addProductionControl(p->controlMode);

        p->addProductionControl(Well::ProducerCMode::BHP);
        if (! p->predictionMode) {
            p->BHPTarget.update(0.0);
            p->setBHPLimit(unit_system->to_si(UnitSystem::measure::pressure, rst_well.bhp_target_float));
            p->controlMode = producer_cmode_from_int(rst_well.hist_requested_control);
        }
        else if (this->isAvailableForGroupControl())
            p->addProductionControl(Well::ProducerCMode::GRUP);

        this->updateProduction(std::move(p));
    }
    else {
        auto i = std::make_shared<WellInjectionProperties>(*this->unit_system, wname);
        i->VFPTableNumber = rst_well.vfp_table;
        i->predictionMode = this->prediction_mode;

        if (is_defined(rst_well.wrat_target) ||
            is_defined(rst_well.grat_target) )
            i->addInjectionControl(Well::InjectorCMode::RATE);

        if (is_defined(rst_well.resv_target)) {
            i->reservoirInjectionRate.update(rst_well.resv_target);
            i->addInjectionControl(Well::InjectorCMode::RESV);
        }

        i->injectorType = rst_well.wtype.injector_type();
        switch (i->injectorType) {
        case InjectorType::WATER:
            update_if_defined(i->surfaceInjectionRate, rst_well.wrat_target);
            break;
        case InjectorType::GAS:
            update_if_defined(i->surfaceInjectionRate, rst_well.grat_target);
            break;
        default:
            throw std::invalid_argument("What ...");
        }

        if (is_defined(rst_well.thp_target)) {
            i->THPTarget.update(rst_well.thp_target);
            i->addInjectionControl(Well::InjectorCMode::THP);
        }

        const auto active_control = i->predictionMode
            ? injector_cmode_from_int(rst_well.active_control)
            : injector_cmode_from_int(rst_well.hist_requested_control);

        if (! i->predictionMode) {
            i->clearControls();
            if ((active_control != Well::InjectorCMode::RATE) &&
                (active_control != Well::InjectorCMode::BHP))
            {
                throw std::invalid_argument {
                    fmt::format("Unsupported control mode '{}' for "
                                "history controlled injection well '{}'",
                                WellInjectorCMode2String(active_control), this->name())
                };
            }
        }

        i->controlMode = active_control;
        i->addInjectionControl(active_control);

        i->addInjectionControl(Well::InjectorCMode::BHP);
        i->BHPTarget.update(is_defined(rst_well.bhp_target_float) ? rst_well.bhp_target_float : 0.0);
        if (! i->predictionMode) {
            if (i->controlMode == Well::InjectorCMode::BHP)
                i->bhp_hist_limit = rst_well.hist_bhp_target;
            else  {
                // defaults (ie FBHPDEF) are not stored in the restart file.
                // we thus use an empirically obtained default value.
                i->resetDefaultHistoricalBHPLimit();
            }
        }
        else if (this->isAvailableForGroupControl())
            i->addInjectionControl(Well::InjectorCMode::GRUP);

        this->updateInjection(std::move(i));

        const auto isTemp = (rst_well.inj_temperature < RestartIO::RstWell::UNDEFINED_VALUE);
        std::size_t tracer_conc_index = 0;
        if (isTemp) {
            this->well_inj_temperature = rst_well.inj_temperature;
            ++tracer_conc_index;
        }

        if (!rst_well.tracer_concentration_injection.empty()) {
            auto tracer = std::make_shared<WellTracerProperties>(this->getTracerProperties());
            for (std::size_t tracer_index = 0; tracer_index < tracer_config.size(); ++tracer_index) {
                const auto& tname = tracer_config[tracer_index].name;
                const auto phase = tracer_config[tracer_index].phase;
                if (phase == Phase::WATER) {
                    const auto concentration = rst_well.tracer_concentration_injection[tracer_conc_index];
                    tracer->setConcentration(tname, concentration);
                } else {
                    const auto free_conc = rst_well.tracer_concentration_injection[tracer_conc_index];
                    const auto sol_conc = rst_well.tracer_concentration_injection[++tracer_conc_index];
                    if (WellType::gas_injector(this->wtype.ecl_wtype()) || WellType::oil_injector(this->wtype.ecl_wtype())) {
                        tracer->setConcentration(tname, free_conc);
                        if (sol_conc > 0.0) {
                            OpmLog::warning(fmt::format("Well {}: Restoring a non-zero solution concentration of tracer {} is not yet supported.", rst_well.name, tname));
                        }
                    }
                }
            }
            this->updateTracer(tracer);
        }
    }
}



Well::Well(const std::string& wname_arg,
           const std::string& gname,
           std::size_t init_step_arg,
           std::size_t insert_index_arg,
           int headI_arg,
           int headJ_arg,
           const std::optional<double>& ref_depth_arg,
           const WellType& wtype_arg,
           ProducerCMode whistctl_cmode,
           Connection::Order ordering_arg,
           const UnitSystem& unit_system_arg,
           double udq_undefined_arg,
           double dr,
           bool allow_xflow,
           bool auto_shutin,
           int pvt_table_,
           GasInflowEquation inflow_eq,
           bool temp_option):
    wname(wname_arg),
    group_name(gname),
    init_step(init_step_arg),
    insert_index(insert_index_arg),
    headI(headI_arg),
    headJ(headJ_arg),
    ref_depth(ref_depth_arg),
    drainage_radius(dr),
    allow_cross_flow(allow_xflow),
    automatic_shutin(auto_shutin),
    pvt_table(pvt_table_),
    gas_inflow(inflow_eq),
    unit_system(&unit_system_arg),
    udq_undefined(udq_undefined_arg),
    wtype(wtype_arg),
    guide_rate({true, -1, Well::GuideRateTarget::UNDEFINED,ParserKeywords::WGRUPCON::SCALING_FACTOR::defaultValue}),
    efficiency_factor(1.0),
    use_efficiency_in_network(true),
    solvent_fraction(0.0),
    derive_refdepth_from_conns_{ ! ref_depth_arg.has_value() || (*ref_depth_arg < 0.0) },
    econ_limits(std::make_shared<WellEconProductionLimits>()),
    foam_properties(std::make_shared<WellFoamProperties>()),
    polymer_properties(std::make_shared<WellPolymerProperties>()),
    micp_properties(std::make_shared<WellMICPProperties>()),
    brine_properties(std::make_shared<WellBrineProperties>()),
    tracer_properties(std::make_shared<WellTracerProperties>()),
    connections(std::make_shared<WellConnections>(ordering_arg, headI, headJ)),
    production(std::make_shared<WellProductionProperties>(*unit_system, wname)),
    injection(std::make_shared<WellInjectionProperties>(*unit_system, wname)),
    wvfpdp(std::make_shared<WVFPDP>()),
    wvfpexp(std::make_shared<WVFPEXP>()),
    wdfac(std::make_shared<WDFAC>()),
    status(Status::SHUT),
    well_inj_temperature(std::nullopt),
    well_inj_mult(std::nullopt)
{
    if (temp_option) {
        default_well_inj_temperature = Metric::TemperatureOffset + 0.0;
    }

    auto p = std::make_shared<WellProductionProperties>(*this->unit_system, this->wname);
    p->whistctl_cmode = whistctl_cmode;
    this->updateProduction(p);
}

Well Well::serializationTestObject()
{
    Well result;
    result.wname = "test1";
    result.group_name = "test2";
    result.init_step = 1;
    result.insert_index = 2;
    result.headI = 3;
    result.headJ = 4;
    result.ref_depth = 5;
    result.udq_undefined = 6.0;
    result.status = Status::AUTO;
    result.drainage_radius = 7.0;
    result.allow_cross_flow = true;
    result.automatic_shutin = false;
    result.pvt_table = 77;
    result.gas_inflow = GasInflowEquation::GPP;
    result.wtype = WellType(Phase::WATER);
    result.guide_rate = WellGuideRate::serializationTestObject();
    result.efficiency_factor = 8.0;
    result.use_efficiency_in_network = true;
    result.solvent_fraction = 9.0;
    result.prediction_mode = false;
    result.derive_refdepth_from_conns_ = false;
    result.econ_limits = std::make_shared<Opm::WellEconProductionLimits>(Opm::WellEconProductionLimits::serializationTestObject());
    result.foam_properties = std::make_shared<WellFoamProperties>(WellFoamProperties::serializationTestObject());
    result.polymer_properties =  std::make_shared<WellPolymerProperties>(WellPolymerProperties::serializationTestObject());
    result.micp_properties =  std::make_shared<WellMICPProperties>(WellMICPProperties::serializationTestObject());
    result.brine_properties = std::make_shared<WellBrineProperties>(WellBrineProperties::serializationTestObject());
    result.tracer_properties = std::make_shared<WellTracerProperties>(WellTracerProperties::serializationTestObject());
    result.connections = std::make_shared<WellConnections>(WellConnections::serializationTestObject());
    result.production = std::make_shared<Well::WellProductionProperties>(Well::WellProductionProperties::serializationTestObject());
    result.injection = std::make_shared<Well::WellInjectionProperties>(Well::WellInjectionProperties::serializationTestObject());
    result.segments = std::make_shared<WellSegments>(WellSegments::serializationTestObject());
    result.wvfpdp = std::make_shared<WVFPDP>(WVFPDP::serializationTestObject());
    result.wvfpexp = std::make_shared<WVFPEXP>(WVFPEXP::serializationTestObject());
    result.wdfac = std::make_shared<WDFAC>(WDFAC::serializationTestObject());
    result.m_pavg = PAvg();
    result.well_inj_temperature = 10.0;
    result.default_well_inj_temperature = 0.0;
    result.well_inj_mult = InjMult::serializationTestObject();
    result.m_filter_concentration = UDAValue::serializationTestObject();
    result.lgr_tag = "LGR-test";
    result.ref_type = WellRefinementType::LGR;
    return result;
}

bool Well::updateWPAVE(const PAvg& pavg)
{
    if (this->m_pavg == pavg) {
        return false;
    }

    this->m_pavg = pavg;
    return true;
}

void Well::flag_lgr_well(void)
{
    ref_type = WellRefinementType::LGR;
}

void Well::set_lgr_well_tag(const std::string& lgr_tag_name)
{
    lgr_tag = lgr_tag_name;
}

const std::string& Well::get_lgr_well_tag(void) const
{   
    if (this->ref_type == WellRefinementType::STANDARD)
        throw std::invalid_argument("Well is not an LGR well.");
     return lgr_tag;
}


bool Well::is_lgr_well(void) const
{
    return this->ref_type == WellRefinementType::LGR;
}


bool Well::updateEfficiencyFactor(double efficiency_factor_arg, bool use_efficiency_in_network_arg) {
    bool update = false;
    if (this->efficiency_factor != efficiency_factor_arg) {
        this->efficiency_factor = efficiency_factor_arg;
        update = true;
    }

    if (this->use_efficiency_in_network != use_efficiency_in_network_arg) {
        this->use_efficiency_in_network = use_efficiency_in_network_arg;
        update = true;
    }

    return update;
}

bool Well::updateWellGuideRate(double guide_rate_arg)
{
    if (this->guide_rate.guide_rate != guide_rate_arg) {
        this->guide_rate.guide_rate = guide_rate_arg;
        return true;
    }

    return false;
}


bool Well::updateFoamProperties(std::shared_ptr<WellFoamProperties> foam_properties_arg)
{
    if (this->wtype.producer()) {
        throw std::runtime_error {
            fmt::format("Assigning foam injection properties is "
                        "disallowed for production well {}.", this->name())
        };
    }

    if (*this->foam_properties != *foam_properties_arg) {
        this->foam_properties = std::move(foam_properties_arg);
        return true;
    }

    return false;
}

bool Well::updatePolymerProperties(std::shared_ptr<WellPolymerProperties> polymer_properties_arg)
{
    if (this->wtype.producer()) {
        throw std::runtime_error {
            fmt::format("Assigning polymer injection properties is "
                        "disallowed for production well {}.", this->name())
        };
    }

    if (*this->polymer_properties != *polymer_properties_arg) {
        this->polymer_properties = std::move(polymer_properties_arg);
        return true;
    }

    return false;
}

bool Well::updateMICPProperties(std::shared_ptr<WellMICPProperties> micp_properties_arg)
{
    if (this->wtype.producer()) {
        throw std::runtime_error {
            fmt::format("Assigning MICP injection properties is "
                        "disallowed for production well {}.", this->name())
        };
    }

    if (*this->micp_properties != *micp_properties_arg) {
        this->micp_properties = std::move(micp_properties_arg);
        return true;
    }

    return false;
}

bool Well::updateBrineProperties(std::shared_ptr<WellBrineProperties> brine_properties_arg)
{
    if (this->wtype.producer()) {
        throw std::runtime_error {
            fmt::format("Assigning brine injection properties is "
                        "disallowed for production well {}.", this->name())
        };
    }

    if (*this->brine_properties != *brine_properties_arg) {
        this->brine_properties = std::move(brine_properties_arg);
        return true;
    }

    return false;
}


bool Well::updateEconLimits(std::shared_ptr<WellEconProductionLimits> econ_limits_arg)
{
    if (*this->econ_limits != *econ_limits_arg) {
        this->econ_limits = std::move(econ_limits_arg);
        return true;
    }

    return false;
}

bool Well::updateWVFPDP(std::shared_ptr<WVFPDP> wvfpdp_arg) {
    if (*this->wvfpdp != *wvfpdp_arg) {
        this->wvfpdp = std::move(wvfpdp_arg);
        return true;
    }

    return false;
}

bool Well::updateWVFPEXP(std::shared_ptr<WVFPEXP> wvfpexp_arg)
{
    if (*this->wvfpexp != *wvfpexp_arg) {
        this->wvfpexp = std::move(wvfpexp_arg);
        return true;
    }

    return false;
}

bool Well::updateWDFAC(std::shared_ptr<WDFAC> wdfac_arg)
{
    if (*this->wdfac != *wdfac_arg) {
        this->wdfac = std::move(wdfac_arg);
        return true;
    }

    return false;
}

void Well::switchToProducer()
{
    auto p = std::make_shared<WellInjectionProperties>(this->getInjectionProperties());

    p->BHPTarget.update(0);
    p->dropInjectionControl(Well::InjectorCMode::BHP);

    this->updateInjection(std::move(p));

    this->wtype.update(true);
}


void Well::switchToInjector()
{
    auto p = std::make_shared<WellProductionProperties>(getProductionProperties());

    p->setBHPLimit(0);
    p->dropProductionControl(ProducerCMode::BHP);

    this->updateProduction(std::move(p));
}

bool Well::updateInjection(std::shared_ptr<WellInjectionProperties> injection_arg)
{
    auto update = this->wtype.update(injection_arg->injectorType);
    if (this->wtype.producer()) {
        this->switchToInjector();
        update = true;
    }

    if (*this->injection != *injection_arg) {
        this->injection = std::move(injection_arg);
        update = true;
    }

    return update;
}

bool Well::updateWellProductivityIndex()
{
    return this->connections->prepareWellPIScaling();
}

bool Well::updateHasProduced()
{
    if (this->wtype.producer() && (this->getStatus() == Status::OPEN)) {
        if (this->has_produced) {
            return false;
        }

        this->has_produced = true;
        return true;
    }

    return false;
}

bool Well::updateHasInjected()
{
    if (this->wtype.injector() && (this->getStatus() == Status::OPEN)) {
        if (this->has_injected) {
            return false;
        }

        this->has_injected= true;
        return true;
    }

    return false;
}

bool Well::updateProduction(std::shared_ptr<WellProductionProperties> production_arg)
{
    if (!this->wtype.producer()) {
        this->switchToProducer();
    }

    if (*this->production != *production_arg) {
        this->production = std::move(production_arg);
        return true;
    }

    return false;
}

bool Well::updateTracer(std::shared_ptr<WellTracerProperties> tracer_properties_arg)
{
    if (*this->tracer_properties != *tracer_properties_arg) {
        this->tracer_properties = std::move(tracer_properties_arg);
        return true;
    }

    return false;
}

bool Well::updateWellGuideRate(const bool available,
                               const double guide_rate_arg,
                               const GuideRateTarget guide_phase,
                               const double scale_factor)
{
    bool update = false;
    if (this->guide_rate.available != available) {
        this->guide_rate.available = available;
        update = true;
    }

    if (this->guide_rate.guide_rate != guide_rate_arg) {
        this->guide_rate.guide_rate = guide_rate_arg;
        update = true;
    }

    if (this->guide_rate.guide_phase != guide_phase) {
        this->guide_rate.guide_phase = guide_phase;
        update = true;
    }

    if (this->guide_rate.scale_factor != scale_factor) {
        this->guide_rate.scale_factor = scale_factor;
        update = true;
    }

    return update;
}

bool Well::updateAvailableForGroupControl(const bool available)
{
    bool update = false;

    if (this->guide_rate.available != available) {
        this->guide_rate.available = available;
        update = true;
    }

    return update;
}

bool Well::updateGroup(const std::string& group_arg)
{
    if (this->group_name != group_arg) {
        this->group_name = group_arg;
        return true;
    }

    return false;
}

bool Well::updateHead(std::optional<int> I, std::optional<int> J)
{
    bool update = false;

    if (I.has_value() && (this->headI != *I)) {
        this->headI = *I;
        update = true;
    }

    if (J.has_value() && (this->headJ != *J)) {
        this->headJ = *J;
        update = true;
    }

    return update;
}

bool Well::updateStatus(const Status well_state)
{
    this->status = well_state;
    return true;
}

bool Well::updateRefDepth(std::optional<double> ref_depth_arg)
{
    if (this->ref_depth == ref_depth_arg) {
        return false;
    }

    if (ref_depth_arg.has_value()) {
        this->derive_refdepth_from_conns_ = *ref_depth_arg < 0.0;

        this->ref_depth = ref_depth_arg;
        return true;
    }

    // If we get here, then this->ref_depth.has_value() &&
    // !ref_depth_arg.has_value().  Assign the argument's value--i.e., reset
    // this->ref_depth--if we're supposed to calculate the reference depth
    // from the connections.

    if (this->derive_refdepth_from_conns_) {
        this->ref_depth = ref_depth_arg;
        return true;
    }

    return false;
}

bool Well::updateDrainageRadius(std::optional<double> drainage_radius_arg)
{
    if (drainage_radius_arg.has_value() &&
        (this->drainage_radius != *drainage_radius_arg))
    {
        this->drainage_radius = *drainage_radius_arg;
        return true;
    }

    return false;
}


bool Well::updateCrossFlow(const bool allow_cross_flow_arg)
{
    if (this->allow_cross_flow != allow_cross_flow_arg) {
        this->allow_cross_flow = allow_cross_flow_arg;
        return true;
    }

    return false;
}

bool Well::updateAutoShutin(const bool auto_shutin)
{
    if (this->automatic_shutin != auto_shutin) {
        this->automatic_shutin = auto_shutin;
        return true;
    }

    return false;
}


bool Well::updateConnections(std::shared_ptr<WellConnections> connections_arg, bool force)
{
    connections_arg->order();

    if (force || (*this->connections != *connections_arg)) {
        this->connections = std::move(connections_arg);
        return true;
    }

    return false;
}

bool Well::updateConnections(std::shared_ptr<WellConnections> connections_arg, const ScheduleGrid& grid)
{
    bool update = this->updateConnections(std::move(connections_arg), false);

    if (this->pvt_table == 0 && !this->connections->empty()) {
        const auto& lowest = this->connections->lowest();
        const auto& props = grid.get_cell(lowest.getI(), lowest.getJ(), lowest.getK()).props;
        this->pvt_table = props->pvtnum;
        update = true;
    }

    return update;
}

bool Well::updateSolventFraction(const double solvent_fraction_arg)
{
    if (this->solvent_fraction != solvent_fraction_arg) {
        this->solvent_fraction = solvent_fraction_arg;
        return true;
    }

    return false;
}

bool Well::handleCOMPSEGS(const DeckKeyword& keyword,
                          const ScheduleGrid& grid,
                          const ParseContext& parseContext,
                          ErrorGuard& errors)
{
    if (this->segments == nullptr) {
        throw OpmInputError{
            fmt::format("The well segment structure for well '{}' must be "
                        "specified in the WELSEGS keyword before assigning "
                        "connections to segments in the COMPSEGS keyword.",
                        this->name()),
            keyword.location()
        };
    }

    auto [new_connections, new_segments] = Compsegs::processCOMPSEGS
        (keyword, *this->connections, *this->segments,
         grid, parseContext, errors);

    this->updateConnections(std::make_shared<WellConnections>(std::move(new_connections)), false);
    this->updateSegments(std::make_shared<WellSegments>(std::move(new_segments)));

    return true;
}

const std::string& Well::groupName() const
{
    return this->group_name;
}

bool Well::isMultiSegment() const
{
    return this->segments != nullptr;
}

bool Well::isProducer() const
{
    return this->wtype.producer();
}

bool Well::isInjector() const
{
    return this->wtype.injector();
}

const WellType& Well::wellType() const
{
    return this->wtype;
}

Well::InjectorCMode Well::injection_cmode() const
{
    if (this->isInjector()) {
        return this->injection->controlMode;
    }

    throw std::logic_error {
        fmt::format("Queried for INJECTION cmode for producer: {}", this->name())
    };
}

Well::ProducerCMode Well::production_cmode() const
{
    if (this->isProducer()) {
        return this->production->controlMode;
    }

    throw std::logic_error {
        fmt::format("Queried for PRODUCTION cmode for injector : {}", this->name())
    };
}

InjectorType Well::injectorType() const
{
    if (this->wtype.producer()) {
        throw std::runtime_error("Can not access injectorType attribute of a producer");
    }

    return this->injection->injectorType;
}

bool Well::isAvailableForGroupControl() const
{
    return this->guide_rate.available;
}

double Well::getGuideRate() const
{
    return this->guide_rate.guide_rate;
}

Well::GuideRateTarget Well::getGuideRatePhase() const
{
    const auto target = this->getRawGuideRatePhase();

    if (this->isInjector() && (target == GuideRateTarget::RAT)) {
        return this->preferredPhaseAsGuideRatePhase();
    }

    return target;
}

Well::GuideRateTarget Well::getRawGuideRatePhase() const
{
    return this->guide_rate.guide_phase;
}

Well::GuideRateTarget Well::preferredPhaseAsGuideRatePhase() const
{
    switch (this->getPreferredPhase()) {
    case Phase::OIL:   return GuideRateTarget::OIL;
    case Phase::GAS:   return GuideRateTarget::GAS;
    case Phase::WATER: return GuideRateTarget::WAT;

    default:
        throw std::logic_error {
            fmt::format("Unable to convert well preferred "
                        "phase {} to GuideRate target phase",
                        static_cast<int>(this->getPreferredPhase()))
        };
    }
}

double Well::getGuideRateScalingFactor() const
{
    return this->guide_rate.scale_factor;
}

double Well::getEfficiencyFactor(bool network) const {
    if (network && !this->use_efficiency_in_network) {
        return 1.0;
    }
    return this->efficiency_factor;
}

double Well::getSolventFraction() const {
    return this->solvent_fraction;
}

std::size_t Well::seqIndex() const
{
    return this->insert_index;
}

int Well::getHeadI() const
{
    return this->headI;
}

int Well::getHeadJ() const
{
    return this->headJ;
}

bool Well::getAutomaticShutIn() const
{
    return this->automatic_shutin;
}

bool Well::getAllowCrossFlow() const
{
    return this->allow_cross_flow;
}

bool Well::hasRefDepth() const
{
    return this->ref_depth.has_value();
}

double Well::getRefDepth() const
{
    if (!this->hasRefDepth()) {
        throw std::logic_error {
            fmt::format("Well: {} - tried to access not "
                        "initialized well reference depth", this->name())
        };
    }

    return *this->ref_depth;
}

double Well::getWPaveRefDepth() const
{
    return this->wpave_ref_depth.value_or(this->getRefDepth());
}

void Well::updateRefDepth()
{
    if ((! this->ref_depth.has_value() || (*this->ref_depth < 0.0)) &&
        this->derive_refdepth_from_conns_)
    {
        // Reference depth was defaulted and we get the depth of the first
        // connection.

        if (this->connections->empty()) {
            throw std::invalid_argument {
                fmt::format("No reservoir connection defined for well {}. "
                            "Cannot infer reference depth.", this->name())
            };
        }

        this->ref_depth = this->connections->get(0).depth();
    }
}

void Well::updateWPaveRefDepth(const double depth)
{
    this->wpave_ref_depth = depth;
}

double Well::getDrainageRadius() const
{
    return this->drainage_radius;
}

//const std::vector<std::string>& Well::wListNames() const {
//    return this->w_list_names;
//}

const std::string& Well::name() const
{
    return this->wname;
}

bool Well::hasSameConnectionsPointers(const Well& other) const
{
    // Note: This is *supposed* to be a pointer comparison.  We need to know
    // if the two connection structures represent the exact same object, not
    // just if they have the same value.
    return this->connections == other.connections;
}

void Well::setInsertIndex(const std::size_t index)
{
    this->insert_index = index;
}

double Well::convertDeckPI(double deckPI) const
{
    using M = UnitSystem::measure;

    // XXX: Should really have LIQUID here too, but the 'Phase' type does
    //      not provide that enumerator.
    switch (this->getPreferredPhase()) {
    case Phase::GAS:
        return this->unit_system->to_si(M::gas_productivity_index, deckPI);

    case Phase::OIL:
    case Phase::WATER:
        return this->unit_system->to_si(M::liquid_productivity_index, deckPI);

    default:
        throw std::invalid_argument {
            "Preferred phase " + std::to_string(static_cast<int>(this->getPreferredPhase())) +
            " is not supported. Must be one of 'OIL', 'GAS', or 'WATER'"
        };
    }
}

void Well::applyWellProdIndexScaling(const double scalingFactor,
                                     std::vector<bool>& scalingApplicable)
{
    this->connections->applyWellPIScaling(scalingFactor, scalingApplicable);
}

bool Well::hasConnections() const
{
    return !this->connections->empty();
}

const WellConnections& Well::getConnections() const
{
    return *this->connections;
}

WellConnections& Well::getConnections()
{
    return *this->connections;
}

std::vector<const Connection *> Well::getConnections(int completion) const
{
    std::vector<const Connection *> connvector;

    for (const auto& conn : this->getConnections()) {
        if (conn.complnum() == completion) {
            connvector.push_back(&conn);
        }
    }

    return connvector;
}

const WellFoamProperties& Well::getFoamProperties() const
{
    return *this->foam_properties;
}

const WellPolymerProperties& Well::getPolymerProperties() const
{
    return *this->polymer_properties;
}

const WellMICPProperties& Well::getMICPProperties() const
{
    return *this->micp_properties;
}

const WellBrineProperties& Well::getBrineProperties() const
{
    return *this->brine_properties;
}

const WellTracerProperties& Well::getTracerProperties() const
{
    return *this->tracer_properties;
}

const WVFPDP& Well::getWVFPDP() const
{
    return *this->wvfpdp;
}

const WVFPEXP& Well::getWVFPEXP() const
{
    return *this->wvfpexp;
}

const WDFAC& Well::getWDFAC() const
{
    return *this->wdfac;
}

const WellEconProductionLimits& Well::getEconLimits() const
{
    return *this->econ_limits;
}

const Well::WellProductionProperties& Well::getProductionProperties() const
{
    return *this->production;
}

const WellSegments& Well::getSegments() const
{
    if (this->segments != nullptr) {
        return *this->segments;
    }

    throw std::logic_error {
        fmt::format("Multi-segment information requested "
                    "for regular well '{}'", this->name())
    };
}

int Well::maxSegmentID() const
{
    return (this->segments == nullptr)
        ? 0
        : this->segments->maxSegmentID();
}

int Well::maxBranchID() const
{
    return (this->segments == nullptr)
        ? 0
        : this->segments->maxBranchID();
}

const Well::WellInjectionProperties& Well::getInjectionProperties() const
{
    return *this->injection;
}

Well::Status Well::getStatus() const
{
    return this->status;
}

const PAvg& Well::pavg() const
{
    return this->m_pavg;
}


std::map<int, std::vector<Connection>> Well::getCompletions() const
{
    std::map<int, std::vector<Connection>> completions;

    for (const auto& conn : *this->connections) {
        auto pair = completions.find( conn.complnum() );
        if (pair == completions.end()) {
            completions[conn.complnum()] = {};
        }

        pair = completions.find(conn.complnum());
        pair->second.push_back(conn);
    }

    return completions;
}

bool Well::hasCompletion(int completion) const
{
    return std::any_of(this->connections->begin(), this->connections->end(),
                       [completion](const auto& conn)
                       {
                           return conn.complnum() == completion;
                       });
}

Phase Well::getPreferredPhase() const
{
    return this->wtype.preferred_phase();
}

int Well::pvt_table_number() const
{
    return this->pvt_table;
}

int Well::fip_region_number() const
{
    return ParserKeywords::WELSPECS::FIP_REGION::defaultValue;
}

// When all connections of a well are closed with the WELOPEN keywords, the
// well itself should also be SHUT.  In the main parsing code this is
// handled by the function checkIfAllConnectionsIsShut() which is called at
// the end of every report step in Schedule::iterateScheduleSection().  This
// is done in this way because there is some twisted logic aggregating
// connection changes over a complete report step.
//
// However, when WELOPEN is called at runtime, typically in an ACTIONX
// block, the full Schedule::iterateScheduleSection() is not run and the
// check for all connections closed is not performed.  Therefore, we have a
// runtime flag here which makes sure to close the well in this case.

bool Well::handleWELOPENConnections(const DeckRecord& record,
                                    const Connection::State state_arg)
{
    auto match = [&record](const Connection &c) -> bool
    {
        if (!match_eq(c.getI(), record, "I" , -1)) { return false; }
        if (!match_eq(c.getJ(), record, "J" , -1)) { return false; }
        if (!match_eq(c.getK(), record, "K" , -1)) { return false; }
        if (!match_ge(c.complnum(), record, "C1")) { return false; }
        if (!match_le(c.complnum(), record, "C2")) { return false; }

        return true;
    };

    auto new_connections = std::make_shared<WellConnections>
        (this->connections->ordering(), this->headI, this->headJ);

    for (const auto& connection : *this->connections) {
        if (! match(connection)) {
            // No state change needed here.  Include connection as-is into
            // new connection set.
            new_connections->add(connection);
            continue;
        }

        auto connection_copy = connection;
        connection_copy.setState(state_arg);

        new_connections->add(connection_copy);
    }

    return this->updateConnections(std::move(new_connections), false);
}

bool Well::handleCSKIN(const DeckRecord&      record,
                       const KeywordLocation& location)
{
    using Kw = ParserKeywords::CSKIN;

    auto need_skin_adjustment = [&record](const Connection &c) {
        const auto value_shift = -1;

        return match_eq(c.getI(), record, Kw::I::itemName,       value_shift)
            && match_eq(c.getJ(), record, Kw::J::itemName,       value_shift)
            && match_ge(c.getK(), record, Kw::K_UPPER::itemName, value_shift)
            && match_le(c.getK(), record, Kw::K_LOWER::itemName, value_shift);
    };

    // New connection set which will be updated with new connection level
    // skin factors.
    auto new_connections = std::make_shared<WellConnections>
        (this->connections->ordering(), this->headI, this->headJ);

    const auto skin_factor = record.getItem<Kw::CONNECTION_SKIN_FACTOR>().getSIDouble(0);
    for (const auto& connection : *this->connections) {
        if (! need_skin_adjustment(connection)) {
            // No CSKIN adjustment needed here.  Include connection as-is
            // into new connection set.
            new_connections->add(connection);
            continue;
        }

        // If we get here, we must make the connection's skin factor be
        // 'skin_factor'.
        //
        // First guard against this adjustment making the CTF go negative,
        // typically because the 'skin_factor' value is large and negative
        // itself.
        if (const auto& ctf_props = connection.ctfProperties();
            ctf_props.peaceman_denom + skin_factor - ctf_props.skin_factor < 0.0)
        {
            throw OpmInputError {
                fmt::format("Negative connection transmissibility "
                            "factor generated by skin factor {} "
                            "in connection ({},{},{}) for well {}.",
                            skin_factor,
                            connection.getI() + 1,
                            connection.getJ() + 1,
                            connection.getK() + 1,
                            this->name()),
                location
            };
        }

        auto connection_copy = connection;
        connection_copy.setSkinFactor(skin_factor);

        new_connections->add(connection_copy);
    }

    return this->updateConnections(std::move(new_connections), false);
}

bool Well::handleCOMPLUMP(const DeckRecord& record)
{
    auto match = [&record](const Connection &c) -> bool {
        if (!match_eq(c.getI(), record, "I" , -1)) { return false; }
        if (!match_eq(c.getJ(), record, "J" , -1)) { return false; }
        if (!match_ge(c.getK(), record, "K1", -1)) { return false; }
        if (!match_le(c.getK(), record, "K2", -1)) { return false; }

        return true;
    };

    auto new_connections = std::make_shared<WellConnections>
        (this->connections->ordering(), this->headI, this->headJ);

    const int complnum = record.getItem("N").get<int>(0);
    if (complnum <= 0) {
        throw std::invalid_argument {
            fmt::format("Completion number must be >= 1. COMPLNUM={} is invalid", complnum)
        };
    }

    for (const auto& connection : *this->connections) {
        if (! match(connection)) {
            new_connections->add(connection);
            continue;
        }

        auto connection_copy = connection;
        connection_copy.setComplnum(complnum);

        new_connections->add(connection_copy);
    }

    return this->updateConnections(std::move(new_connections), false);
}

bool Well::handleWPIMULT(const DeckRecord& record)
{
    auto match = [&record](const Connection &c) -> bool {
        if (!match_ge(c.complnum(), record, "FIRST")) { return false; }
        if (!match_le(c.complnum(), record, "LAST"))  { return false; }
        if (!match_eq(c.getI()    , record, "I", -1)) { return false; }
        if (!match_eq(c.getJ()    , record, "J", -1)) { return false; }
        if (!match_eq(c.getK()    , record, "K", -1)) { return false; }

        return true;
    };

    auto new_connections = std::make_shared<WellConnections>
        (this->connections->ordering(), this->headI, this->headJ);

    const auto wellPi = record.getItem("WELLPI").get<double>(0);

    for (const auto& connection : *this->connections) {
        if (! match(connection)) {
            new_connections->add(connection);
            continue;
        }

        auto connection_copy = connection;
        connection_copy.scaleWellPi(wellPi);

        new_connections->add(connection_copy);
    }

    return this->updateConnections(std::move(new_connections), false);
}

bool Well::handleWINJCLN(const DeckRecord& record, const KeywordLocation& location)
{
    const double fraction_removal = record
        .getItem<ParserKeywords::WINJCLN::FRAC_REMOVE>()
        .getSIDouble(0);

    if ((fraction_removal < 0.0) || (fraction_removal > 1.0)) {
        const auto reason = fmt::format("Item 2 in keyword WINJCLN "
                                        "must be between 0 and 1, while a "
                                        "value {} is given.", fraction_removal);

        throw OpmInputError(reason, location);
    }

    auto match = [&record](const Connection& c) -> bool {
        if (!match_eq(c.getI(), record, "I", -1)) { return false; }
        if (!match_eq(c.getJ(), record, "J", -1)) { return false; }
        if (!match_eq(c.getK(), record, "K", -1)) { return false; }

        return true;
    };

    auto new_connections = std::make_shared<WellConnections>
        (this->connections->ordering(), this->headI, this->headJ);

    for (const auto& connection : *this->connections) {
        if (! match(connection)) {
            new_connections->add(connection);
            continue;
        }

        auto filter_cake = connection.getFilterCake();
        filter_cake.applyCleanMultiplier(1.0 - fraction_removal);

        auto connection_copy = connection;
        connection_copy.setFilterCake(filter_cake);

        new_connections->add(connection_copy);
    }

    return this->updateConnections(std::move(new_connections), false);
}

bool Well::handleWINJDAM(const DeckRecord& record, const KeywordLocation& location)
{
    auto match = [&record](const Connection& c) -> bool {
        if (!match_eq(c.getI(), record, "I", -1)) { return false; }
        if (!match_eq(c.getJ(), record, "J", -1)) { return false; }
        if (!match_eq(c.getK(), record, "K", -1)) { return false; }

        return true;
    };

    const FilterCake filter_cake {record, location};

    auto new_connections = std::make_shared<WellConnections>
        (this->connections->ordering(), this->headI, this->headJ);

    for (const auto& connection : *this->connections) {
        if (! match(connection)) {
            new_connections->add(connection);
            continue;
        }

        auto connection_copy = connection;
        connection_copy.setFilterCake(filter_cake);

        new_connections->add(connection_copy);
    }

    return this->updateConnections(std::move(new_connections), false);
}

bool Well::handleWINJMULT(const Opm::DeckRecord& record, const KeywordLocation& location)
{
    // For this keyword, the default for I, J, K will be negative
    //
    // It is not totally clear how specifying 0 or a negative value will
    // work current match_eq function only treats 0 and default values for
    // all connections, we might need to revisit this part later when
    // complication regarding this occurs.  it is possible that changing
    // (item.get<int>(0) == 0); to (item.get<int>(0) <= 0) is solution to go
    // while it remains to be discussed.
    auto match = [&record](const Connection& c) -> bool {
        if (!match_eq(c.getI(), record, "I", -1)) { return false; }
        if (!match_eq(c.getJ(), record, "J", -1)) { return false; }
        if (!match_eq(c.getK(), record, "K", -1)) { return false; }

        return true;
    };

    using Kw = ParserKeywords::WINJMULT;

    const InjMultMode mode = InjMult::injMultModeFromString
        (record.getItem<Kw::MODE>().getTrimmedString(0), location);

    const bool mode_change = (this->inj_mult_mode != mode);
    if (mode_change) {
        this->inj_mult_mode = mode;
    }

    const InjMult inj_mult {
        record.getItem<Kw::FRACTURING_PRESSURE>().getSIDouble(0),
        record.getItem<Kw::MULTIPLIER_GRADIENT>().getSIDouble(0)
    };

    bool connections_update = false;
    bool well_inj_update = false;

    if (mode == InjMultMode::WREV) {
        // All the connections will share the same INJMULT setup when under
        // WREV it is stored in the Well object
        this->well_inj_mult = inj_mult;
        well_inj_update = true;
    }
    else if ((mode == InjMultMode::CREV) ||
             (mode == InjMultMode::CIRR))
    {
        auto new_connections = std::make_shared<WellConnections>
            (this->connections->ordering(), this->headI, this->headJ);

        for (const auto& connection : *this->connections) {
            if (! match(connection)) {
                new_connections->add(connection);
                continue;
            }

            auto connection_copy = connection;
            connection_copy.setInjMult(inj_mult);

            new_connections->add(connection_copy);
        }

        connections_update = this->updateConnections(std::move(new_connections), false);
    }

    return mode_change || connections_update || well_inj_update;
}


bool Opm::Well::applyGlobalWPIMULT(const double scaling_factor)
{
    auto new_connections = std::make_shared<WellConnections>
        (this->connections->ordering(), this->headI, this->headJ);

    for (const auto& connection : *this->connections) {
        auto connection_copy = connection;

        connection_copy.scaleWellPi(scaling_factor);
        new_connections->add(connection_copy);
    }

    return this->updateConnections(std::move(new_connections), false);
}

void Well::updateSegments(std::shared_ptr<WellSegments> segments_arg)
{
    this->segments = std::move(segments_arg);
    this->updateRefDepth(this->segments->depthTopSegment());
    this->derive_refdepth_from_conns_ = false;
}

bool Well::handleWELSEGS(const DeckKeyword& keyword)
{
    if (this->segments != nullptr) {
        auto new_segments = std::make_shared<WellSegments>(*this->segments);
        new_segments->loadWELSEGS(keyword);

        this->updateSegments(std::move(new_segments));
    }
    else {
        this->updateSegments(std::make_shared<WellSegments>(keyword));
    }

    return true;
}

bool Well::updatePVTTable(std::optional<int> pvt_table_)
{
    if (pvt_table_.has_value() && (this->pvt_table != *pvt_table_)) {
        this->pvt_table = *pvt_table_;
        return true;
    }

    return false;
}

bool Well::updateWSEGSICD(const std::vector<std::pair<int, SICD>>& sicd_pairs)
{
    auto new_segments = std::make_shared<WellSegments>(*this->segments);

    if (new_segments->updateWSEGSICD(sicd_pairs)) {
        this->segments = std::move(new_segments);
        return true;
    }

    return false;
}

bool Well::updateWSEGAICD(const std::vector<std::pair<int, AutoICD>>& aicd_pairs,
                          const KeywordLocation& location)
{
    auto new_segments = std::make_shared<WellSegments>(*this->segments);

    if (new_segments->updateWSEGAICD(aicd_pairs, location)) {
        this->segments = std::move(new_segments);
        return true;
    }

    return false;
}

bool Well::updateWSEGVALV(const std::vector<std::pair<int, Valve>>& valve_pairs)
{
    auto new_segments = std::make_shared<WellSegments>(*this->segments);

    if (new_segments->updateWSEGVALV(valve_pairs)) {
        this->segments = std::move(new_segments);
        return true;
    }

    return false;
}

void Well::filterConnections(const ActiveGridCells& grid)
{
    this->connections->filter(grid);
}

std::size_t Well::firstTimeStep() const
{
    return this->init_step;
}

bool Well::hasBeenDefined(const std::size_t timeStep) const
{
    return timeStep >= this->init_step;
}

Well::GasInflowEquation Well::gas_inflow_equation() const
{
    return this->gas_inflow;
}

bool Well::predictionMode() const
{
    return this->prediction_mode;
}

bool Well::hasProduced() const
{
    return this->has_produced;
}

bool Well::hasInjected() const
{
    return this->has_injected;
}

bool Well::updatePrediction(const bool prediction_mode_arg)
{
    if (this->prediction_mode != prediction_mode_arg) {
        this->prediction_mode = prediction_mode_arg;
        return true;
    }

    return false;
}

double Well::production_rate(const SummaryState& st, const Phase prod_phase) const
{
    if (! this->isProducer()) {
        return 0.0;
    }

    const auto controls = this->productionControls(st);

    auto zero_if_undefined = [&st](const double value)
    {
        return st.is_undefined_value(value)
            ? 0.0 : value;
    };

    switch (prod_phase) {
    case Phase::WATER: return zero_if_undefined(controls.water_rate);
    case Phase::OIL:   return zero_if_undefined(controls.oil_rate);
    case Phase::GAS:   return zero_if_undefined(controls.gas_rate);

    case Phase::SOLVENT:
        throw std::invalid_argument( "Production of 'SOLVENT' requested." );

    case Phase::POLYMER:
        throw std::invalid_argument( "Production of 'POLYMER' requested." );

    case Phase::ENERGY:
        throw std::invalid_argument( "Production of 'ENERGY' requested." );

    case Phase::POLYMW:
        throw std::invalid_argument( "Production of 'POLYMW' requested.");

    case Phase::FOAM:
        throw std::invalid_argument( "Production of 'FOAM' requested.");

    case Phase::BRINE:
        throw std::invalid_argument( "Production of 'BRINE' requested.");

    case Phase::ZFRACTION:
        throw std::invalid_argument( "Production of 'ZFRACTION' requested.");
    }

    throw std::logic_error( "Unreachable state. Invalid Phase value. "
                            "This is likely a programming error." );
}

double Well::injection_rate(const SummaryState& st, Phase phase_arg) const
{
    if (!this->isInjector()) {
        return 0.0;
    }

    const auto controls = this->injectionControls(st);

    const auto type = controls.injector_type;

    if ((phase_arg == Phase::WATER) && (type != InjectorType::WATER)) { return 0.0; }
    if ((phase_arg == Phase::OIL)   && (type != InjectorType::OIL)  ) { return 0.0; }
    if ((phase_arg == Phase::GAS)   && (type != InjectorType::GAS)  ) { return 0.0; }

    if (st.is_undefined_value(controls.surface_rate)) {
        return 0.0;
    }

    return controls.surface_rate;
}


bool Well::wellNameInWellNamePattern(const std::string& wellName,
                                     const std::string& wellNamePattern)
{
    return shmatch(wellNamePattern, wellName);
}

Well::ProductionControls Well::productionControls(const SummaryState& st) const
{
    if (this->isProducer()) {
        return this->production->controls(st, this->udq_undefined);
    }

    throw std::logic_error("Trying to get production data from an injector");
}

Well::InjectionControls Well::injectionControls(const SummaryState& st) const
{
    if (!this->isProducer()) {
        return this->injection->controls(*this->unit_system, st, this->udq_undefined);
    }

    throw std::logic_error("Trying to get injection data from a producer");
}

double Well::alq_value(const SummaryState& st) const
{
    if (this->wtype.producer()) {
        auto controls = this->production->controls(st, this->udq_undefined);
        return controls.alq_value;
    }

    throw std::runtime_error("Cannot request ALQ value in an injector");
}

// These accessor functions are at the "wrong" level of abstraction; the
// same properties are part of the InjectionControls and ProductionControls
// structs.  They are made available here to avoid passing a SummaryState
// instance in situations where it is not really needed.

int Well::vfp_table_number() const
{
    return this->wtype.producer()
        ? this->production->VFPTableNumber
        : this->injection->VFPTableNumber;
}

double Well::inj_temperature() const
{
    if (!this->wtype.injector()) {
        throw std::logic_error {
            fmt::format("Well {}: Cannot request injection "
                        "temperature for a non-injector", this->name())
        };
    }

    if (!this->well_inj_temperature) {
        if (this->default_well_inj_temperature) {
            OpmLog::warning(fmt::format("Well {}: Injection temperature not specified, "
                                        "using default value of {}", this->name(),
                                        *this->default_well_inj_temperature));

            return *this->default_well_inj_temperature;
        }
        else {
            throw std::logic_error {
                fmt::format("Well {}: Unable to obtain injection temperature - "
                            "not specified in deck and no default defined.", this->name())
            };
        }
    }

    return *this->well_inj_temperature;
}

bool Well::hasInjTemperature() const
{
    return this->well_inj_temperature.has_value();
}

void Well::setWellInjTemperature(const double temp)
{
    this->well_inj_temperature = temp;
}

bool Well::cmp_structure(const Well& other) const
{
    if ((this->segments == nullptr) != (other.segments == nullptr)) {
        return false;
    }

    if ((this->segments != nullptr) && (this->getSegments() != other.getSegments())) {
        return false;
    }

    if ((this->unit_system != nullptr) && (*this->unit_system != *other.unit_system)) {
        return false;
    }

    return (this->name() == other.name())
        && (this->groupName() == other.groupName())
        && (this->firstTimeStep() == other.firstTimeStep())
        && (this->seqIndex() == other.seqIndex())
        && (this->getHeadI() == other.getHeadI())
        && (this->getHeadJ() == other.getHeadJ())
        && (this->hasRefDepth() == other.hasRefDepth())
        && (!this->hasRefDepth() || (this->getRefDepth() == other.getRefDepth()))
        && (this->getDrainageRadius() == other.getDrainageRadius())
        && (this->getAllowCrossFlow() == other.getAllowCrossFlow())
        && (this->getAutomaticShutIn() == other.getAutomaticShutIn())
        && (this->udq_undefined == other.udq_undefined)
        && (this->getPreferredPhase() == other.getPreferredPhase()) // wellType()
        && (this->efficiency_factor == other.efficiency_factor)
        && (this->use_efficiency_in_network == other.use_efficiency_in_network)
        && (this->getConnections() == other.getConnections())
        ;
}


bool Well::operator==(const Well& data) const
{
    return this->cmp_structure(data)
        && (this->wpave_ref_depth == data.wpave_ref_depth)
        && (this->gas_inflow == data.gas_inflow)
        && (this->pvt_table == data.pvt_table)
        && (this->isProducer() == data.isProducer()) // wellType()
        && (this->guide_rate == data.guide_rate)
        && (this->hasProduced() == data.hasProduced())
        && (this->hasInjected() == data.hasInjected())
        && (this->predictionMode() == data.predictionMode())
        && (this->getSolventFraction() == data.getSolventFraction())
        && (this->derive_refdepth_from_conns_ == data.derive_refdepth_from_conns_)
        && (this->getEconLimits() == data.getEconLimits())
        && (this->getFoamProperties() == data.getFoamProperties())
        && (this->getPolymerProperties() == data.getPolymerProperties())
        && (this->getMICPProperties() == data.getMICPProperties())
        && (this->getBrineProperties() == data.getBrineProperties())
        && (this->getTracerProperties() == data.getTracerProperties())
        && (this->getProductionProperties() == data.getProductionProperties())
        && (this->getInjectionProperties() == data.getInjectionProperties())
        && (this->getWVFPDP() == data.getWVFPDP())
        && (this->getWVFPEXP() == data.getWVFPEXP())
        && (this->getWDFAC() == data.getWDFAC())
        && (this->getStatus() == data.getStatus())
        && (this->m_pavg == data.m_pavg)
        && (this->well_inj_temperature == data.well_inj_temperature)
        && (this->inj_mult_mode == data.inj_mult_mode)
        && (this->default_well_inj_temperature == data.default_well_inj_temperature)
        && (this->well_inj_mult == data.well_inj_mult)
        && (this->m_filter_concentration == data.m_filter_concentration)
        && (this->lgr_tag == data.lgr_tag)
        && (this->ref_type == data.ref_type)

        ;
}

} // namespace Opm

int Opm::Well::eclipseControlMode(const Well::InjectorCMode imode,
                                  const InjectorType        itype)
{
    using IMode = ::Opm::Well::InjectorCMode;
    using Val   = ::Opm::RestartIO::Helpers::VectorItems::IWell::Value::WellCtrlMode;
    using IType = ::Opm::InjectorType;

    switch (imode) {
    case IMode::RATE: {
        switch (itype) {
        case IType::OIL:   return Val::OilRate;
        case IType::WATER: return Val::WatRate;
        case IType::GAS:   return Val::GasRate;
        case IType::MULTI: return Val::WMCtlUnk;
        }}
        break;

    case IMode::RESV: return Val::ResVRate;
    case IMode::THP:  return Val::THP;
    case IMode::BHP:  return Val::BHP;
    case IMode::GRUP: return Val::Group;
    default:
        return Val::WMCtlUnk;
    }

    return Val::WMCtlUnk;
}

int Opm::Well::eclipseControlMode(const Opm::Well::ProducerCMode pmode)
{
    using PMode = ::Opm::Well::ProducerCMode;
    using Val   = ::Opm::RestartIO::Helpers::VectorItems::IWell::Value::WellCtrlMode;

    switch (pmode) {
    case PMode::ORAT: return Val::OilRate;
    case PMode::WRAT: return Val::WatRate;
    case PMode::GRAT: return Val::GasRate;
    case PMode::LRAT: return Val::LiqRate;
    case PMode::RESV: return Val::ResVRate;
    case PMode::THP:  return Val::THP;
    case PMode::BHP:  return Val::BHP;
    case PMode::CRAT: return Val::CombRate;
    case PMode::GRUP: return Val::Group;
    default:
        return Val::WMCtlUnk;
    }

    return Val::WMCtlUnk;
}

// This function converts OPM well status values to an integer value
// suitable for output to the restart file.  OPM tracks the status and the
// active control of a well separately, but when this is written to a
// restart file they are combined to a single integer.  Moreover, OPM
// permits a well to have an active control while still being shut, but when
// this is converted to an integer value suitable for the restart file, the
// value 0 will be used to signal a SHUT well and the active control will be
// lost.
//
// In the case of a well which is in state 'STOP' or 'AUTO' an integer
// corresponding to the currently active control is written to the restart
// file.

int Opm::Well::eclipseControlMode(const Well&         well,
                                  const SummaryState& st)
{
    if (well.isProducer()) {
        const auto& ctrl = well.productionControls(st);

        return eclipseControlMode(ctrl.cmode);
    }
    else { // Injector
        const auto& ctrl = well.injectionControls(st);

        return eclipseControlMode(ctrl.cmode, well.injectorType());
    }
}

Opm::Well::InjMultMode Opm::Well::getInjMultMode() const
{
    return this->inj_mult_mode;
}

const Opm::InjMult& Opm::Well::getWellInjMult() const
{
    assert(this->aciveWellInjMult());
    return this->well_inj_mult.value();
}

bool Opm::Well::aciveWellInjMult() const
{
   return this->well_inj_mult.has_value();
}

void Opm::Well::setFilterConc(const UDAValue& conc)
{
    this->m_filter_concentration = conc;
}

double Opm::Well::evalFilterConc(const SummaryState& summary_sate) const
{
    return UDA::eval_well_uda(this->m_filter_concentration,
                              this->name(), summary_sate, 0.0);
}
