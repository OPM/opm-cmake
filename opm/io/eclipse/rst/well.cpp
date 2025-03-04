/*
  Copyright 2020 Equinor ASA.

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

#include <opm/io/eclipse/rst/well.hpp>

#include <opm/io/eclipse/rst/header.hpp>
#include <opm/io/eclipse/rst/connection.hpp>

#include <opm/output/eclipse/VectorItems/connection.hpp>
#include <opm/output/eclipse/VectorItems/msw.hpp>
#include <opm/output/eclipse/VectorItems/well.hpp>

#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <opm/common/utility/String.hpp>

#include <opm/input/eclipse/Parser/ParserItem.hpp>
#include <opm/input/eclipse/Parser/ParserKeyword.hpp>
#include <opm/input/eclipse/Parser/ParserRecord.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/W.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {
    bool is_sentinel(const float raw_value)
    {
        return ! (std::abs(raw_value) < Opm::RestartIO::RstWell::UNDEFINED_VALUE);
    }

    double swel_value(const float raw_value)
    {
        return is_sentinel(raw_value) ? 0.0 : raw_value;
    }

    template <typename Convert>
    double keep_sentinel(const float raw_value, Convert&& convert)
    {
        return is_sentinel(raw_value) ? raw_value : convert(raw_value);
    }

    float dfactor_correlation_coefficient_a(const Opm::UnitSystem& unit_system,
                                            const float            coeff_a)
    {
        const auto dimension = Opm::ParserKeywords::WDFACCOR{}
            .getRecord(0).get(Opm::ParserKeywords::WDFACCOR::A::itemName)
            .dimensions().front();

        return static_cast<float>(unit_system.to_si(dimension, coeff_a));
    }

    constexpr int def_ecl_phase = 1;
    constexpr int def_pvt_table = 0;
} // Anonymous namespace

namespace VI = ::Opm::RestartIO::Helpers::VectorItems;
using M = ::Opm::UnitSystem::measure;

Opm::RestartIO::RstWell::RstWell(const UnitSystem&  unit_system,
                                 const RstHeader&   header,
                                 const std::string& group_arg,
                                 const std::string* zwel,
                                 const int*         iwel,
                                 const float*       swel,
                                 const double*      xwel,
                                 const int*         icon,
                                 const float*       scon,
                                 const double*      xcon) :
    name(rtrim_copy(zwel[0])),
    group(group_arg),
    ij(                                                              {iwel[VI::IWell::IHead] - 1, iwel[VI::IWell::JHead] - 1}),
    k1k2(                                                            std::make_pair(iwel[VI::IWell::FirstK] - 1, iwel[VI::IWell::LastK] - 1)),
    wtype(                                                           iwel[VI::IWell::WType], def_ecl_phase),
    well_status(                                                     iwel[VI::IWell::Status]),
    active_control(                                                  iwel[VI::IWell::ActWCtrl]),
    vfp_table(                                                       iwel[VI::IWell::VFPTab]),
    econ_workover_procedure(                                         iwel[VI::IWell::EconWorkoverProcedure]),
    preferred_phase(                                                 iwel[VI::IWell::PreferredPhase]),
    allow_xflow(                                                     iwel[VI::IWell::XFlow] == 1),
    group_controllable_flag(                                         iwel[VI::IWell::WGrupConControllable]),
    econ_limit_end_run(                                              iwel[VI::IWell::EconLimitEndRun]),
    grupcon_gr_phase(                                                iwel[VI::IWell::WGrupConGRPhase]),
    hist_requested_control(                                          iwel[VI::IWell::HistReqWCtrl]),
    msw_index(                                                       iwel[VI::IWell::MsWID]),
    completion_ordering(                                             iwel[VI::IWell::CompOrd]),
    pvt_table(                                                       def_pvt_table),
    msw_pressure_drop_model(                                         iwel[VI::IWell::MSW_PlossMod]),
    wtest_config_reasons(                                            iwel[VI::IWell::WTestConfigReason]),
    wtest_close_reason(                                              iwel[VI::IWell::WTestCloseReason]),
    wtest_remaining(                                                 iwel[VI::IWell::WTestRemaining] -1),
    econ_limit_quantity(                                             iwel[VI::IWell::EconLimitQuantity]),
    econ_workover_procedure_2(                                       iwel[VI::IWell::EconWorkoverProcedure_2]),
    thp_lookup_procedure_vfptable(                                   iwel[VI::IWell::THPLookupVFPTable]),
    close_if_thp_stabilised(                                         iwel[VI::IWell::CloseWellIfTHPStabilised]),
    prevent_thpctrl_if_unstable(                                     iwel[VI::IWell::PreventTHPIfUnstable]),
    glift_active(                                                    iwel[VI::IWell::LiftOpt] == 1),
    glift_alloc_extra_gas(                                           iwel[VI::IWell::LiftOptAllocExtra] == 1),
    // The values orat_target -> bhp_target_float will be used in UDA values. The
    // UDA values are responsible for unit conversion and raw values are
    // internalized here.
    orat_target(                                                     swel[VI::SWell::OilRateTarget]),
    wrat_target(                                                     swel[VI::SWell::WatRateTarget]),
    grat_target(                                                     swel[VI::SWell::GasRateTarget]),
    lrat_target(                                                     swel[VI::SWell::LiqRateTarget]),
    resv_target(                                                     swel[VI::SWell::ResVRateTarget]),
    thp_target(                                                      swel[VI::SWell::THPTarget]),
    bhp_target_float(                                                swel[VI::SWell::BHPTarget]),
    vfp_bhp_adjustment(  unit_system.to_si(M::pressure,              swel[VI::SWell::VfpBhpAdjustment])),
    vfp_bhp_scaling_factor(                                          swel[VI::SWell::VfpBhpScalingFact]),
    hist_lrat_target(    unit_system.to_si(M::liquid_surface_rate,   swel[VI::SWell::HistLiqRateTarget])),
    hist_grat_target(    unit_system.to_si(M::gas_surface_rate,      swel[VI::SWell::HistGasRateTarget])),
    hist_bhp_target(     unit_system.to_si(M::pressure,              swel[VI::SWell::HistBHPTarget])),
    datum_depth(         keep_sentinel(swel[VI::SWell::DatumDepth], [&unit_system](const double depth) { return unit_system.to_si(M::length, depth); })),
    drainage_radius(     unit_system.to_si(M::length,                swel_value(swel[VI::SWell::DrainageRadius]))),
    grupcon_gr_value(                                                swel[VI::SWell::WGrupConGuideRate]), // No unit conversion
    efficiency_factor(   unit_system.to_si(M::identity,              swel[VI::SWell::EfficiencyFactor1])),
    alq_value(                                                       swel[VI::SWell::Alq_value]),
    econ_limit_min_oil(  unit_system.to_si(M::liquid_surface_rate,   swel[VI::SWell::EconLimitMinOil])),
    econ_limit_min_gas(  unit_system.to_si(M::gas_surface_rate,      swel[VI::SWell::EconLimitMinGas])),
    econ_limit_max_wct(  keep_sentinel(swel[VI::SWell::EconLimitMaxWct], [](const double wct_) { return wct_; })),
    econ_limit_max_gor(  keep_sentinel(swel[VI::SWell::EconLimitMaxGor], [&unit_system](const double gor_) { return unit_system.to_si(M::gas_oil_ratio, gor_); })),
    econ_limit_max_wgr(  keep_sentinel(swel[VI::SWell::EconLimitMaxWgr], [&unit_system](const double wgr_) { return unit_system.to_si(M::oil_gas_ratio, wgr_); })),
    econ_limit_max_wct_2(keep_sentinel(swel[VI::SWell::EconLimitMaxWct_2], [](const double wct_) { return wct_; })),
    econ_limit_min_liq(  unit_system.to_si(M::liquid_surface_rate,   swel[VI::SWell::EconLimitMinLiq])),
    wtest_interval(      unit_system.to_si(M::time,                  swel[VI::SWell::WTestInterval])),
    wtest_startup(       unit_system.to_si(M::time,                  swel[VI::SWell::WTestStartupTime])),
    grupcon_gr_scaling(                                              swel[VI::SWell::WGrupConGRScaling]),
    glift_max_rate(      unit_system.to_si(M::gas_surface_rate,      swel[VI::SWell::LOmaxRate])),
    glift_min_rate(      unit_system.to_si(M::gas_surface_rate,      swel[VI::SWell::LOminRate])),
    glift_weight_factor(                                             swel[VI::SWell::LOweightFac]),
    glift_inc_weight_factor(                                         swel[VI::SWell::LOincFac]),
    dfac_corr_coeff_a(dfactor_correlation_coefficient_a(unit_system, swel[VI::SWell::DFacCorrCoeffA])),
    dfac_corr_exponent_b(                                            swel[VI::SWell::DFacCorrExpB]),
    dfac_corr_exponent_c(                                            swel[VI::SWell::DFacCorrExpC]),
    //
    oil_rate(            unit_system.to_si(M::liquid_surface_rate,   xwel[VI::XWell::OilPrRate])),
    water_rate(          unit_system.to_si(M::liquid_surface_rate,   xwel[VI::XWell::WatPrRate])),
    gas_rate(            unit_system.to_si(M::gas_surface_rate,      xwel[VI::XWell::GasPrRate])),
    liquid_rate(         unit_system.to_si(M::rate,                  xwel[VI::XWell::LiqPrRate])),
    void_rate(           unit_system.to_si(M::rate,                  xwel[VI::XWell::VoidPrRate])),
    thp(                 unit_system.to_si(M::pressure,              xwel[VI::XWell::TubHeadPr])),
    flow_bhp(            unit_system.to_si(M::pressure,              xwel[VI::XWell::FlowBHP])),
    wct(                 unit_system.to_si(M::water_cut,             xwel[VI::XWell::WatCut])),
    gor(                 unit_system.to_si(M::gas_oil_ratio,         xwel[VI::XWell::GORatio])),
    oil_total(           unit_system.to_si(M::liquid_surface_volume, xwel[VI::XWell::OilPrTotal])),
    water_total(         unit_system.to_si(M::liquid_surface_volume, xwel[VI::XWell::WatPrTotal])),
    gas_total(           unit_system.to_si(M::gas_surface_volume,    xwel[VI::XWell::GasPrTotal])),
    void_total(          unit_system.to_si(M::volume,                xwel[VI::XWell::VoidPrTotal])),
    water_inj_total(     unit_system.to_si(M::liquid_surface_volume, xwel[VI::XWell::WatInjTotal])),
    gas_inj_total(       unit_system.to_si(M::gas_surface_volume,    xwel[VI::XWell::GasInjTotal])),
    void_inj_total(      unit_system.to_si(M::volume,                xwel[VI::XWell::VoidInjTotal])),
    gas_fvf(                                                         xwel[VI::XWell::GasFVF]),
    bhp_target_double(   unit_system.to_si(M::pressure,              xwel[VI::XWell::BHPTarget])),
    hist_oil_total(      unit_system.to_si(M::liquid_surface_volume, xwel[VI::XWell::HistOilPrTotal])),
    hist_wat_total(      unit_system.to_si(M::liquid_surface_volume, xwel[VI::XWell::HistWatPrTotal])),
    hist_gas_total(      unit_system.to_si(M::gas_surface_volume,    xwel[VI::XWell::HistGasPrTotal])),
    hist_water_inj_total(unit_system.to_si(M::liquid_surface_volume, xwel[VI::XWell::HistWatInjTotal])),
    hist_gas_inj_total(  unit_system.to_si(M::gas_surface_volume,    xwel[VI::XWell::HistGasInjTotal])),
    water_void_rate(     unit_system.to_si(M::liquid_surface_volume, xwel[VI::XWell::WatVoidPrRate])),
    gas_void_rate(       unit_system.to_si(M::gas_surface_volume,    xwel[VI::XWell::GasVoidPrRate]))
{
    // For E100 it appears that +1 instead of -1 is written for group_controllable_flag when the
    // group control is active, so using this to correct active_control (where ind.ctrl. is written)
    if (group_controllable_flag > 0) {
        active_control = VI::IWell::Value::WellCtrlMode::Group;
    }

    // If the TEMP option is active, this is the first tracer
    std::size_t tracer_index = 0;
    const bool isTemp = header.runspec.temp();
    if (isTemp) {
        this->inj_temperature = unit_system.to_si(M::temperature, swel[VI::SWell::TracerOffset + tracer_index++]);
    }
    const auto& tracers = header.runspec.tracers();
    for (const auto num_tracer_injconcs = tracers.water_tracers() + 2*(tracers.gas_tracers() + tracers.oil_tracers()) + tracer_index;
         tracer_index < static_cast<std::size_t>(num_tracer_injconcs);
         ++tracer_index)
    {
        this->tracer_concentration_injection.push_back(swel[VI::SWell::TracerOffset + tracer_index]);
    }

    for (int ic = 0; ic < iwel[VI::IWell::NConn]; ++ic) {
        const std::size_t icon_offset = ic * header.niconz;
        const std::size_t scon_offset = ic * header.nsconz;
        const std::size_t xcon_offset = ic * header.nxconz;

        this->connections.emplace_back(unit_system, ic, header.nsconz,
                                       icon + icon_offset,
                                       scon + scon_offset,
                                       xcon + xcon_offset);
    }
}


Opm::RestartIO::RstWell::RstWell(const UnitSystem&          unit_system,
                                 const RstHeader&           header,
                                 const std::string&         group_arg,
                                 const std::string*         zwel,
                                 const int*                 iwel,
                                 const float*               swel,
                                 const double*              xwel,
                                 const int*                 icon,
                                 const float*               scon,
                                 const double*              xcon,
                                 const std::vector<int>&    iseg,
                                 const std::vector<double>& rseg)
    : RstWell { unit_system, header, group_arg,
                zwel, iwel, swel, xwel,
                icon, scon, xcon }
{
    if (this->msw_index == 0) {
        // Not a multi-segmented well.  Don't create RstSegment objects.
        return;
    }

    // Recall: There are 'nsegmx' segments per MS well in [ir]seg.
    const auto skippedSegments = (this->msw_index - 1) * header.nsegmx;
    const auto* const isegWell = &iseg[skippedSegments * header.nisegz];
    const auto* const rsegWell = &rseg[skippedSegments * header.nrsegz];

    // ------------------------------------------------------------------------

    // 1: Create RstSegment objects for all active segments attatched to this well.
    auto segNumToIx = std::unordered_map<int, std::vector<RstSegment>::size_type>{};
    for (auto is = 0*header.nsegmx; is < header.nsegmx; ++is) {
        const auto* const isegSeg = &isegWell[is * header.nisegz];
        const auto* const rsegSeg = &rsegWell[is * header.nrsegz];

        if (isegSeg[VI::ISeg::BranchNo] > 0) {
            // Segment is on a branch and therefore active.  Create an
            // RstSegment object to represent this segment.
            const auto segNum = is + 1;
            segNumToIx.emplace(segNum, this->segments.size());

            this->segments.emplace_back(unit_system, segNum, isegSeg, rsegSeg);
        }
    }

    // ------------------------------------------------------------------------

    // 2: Compute inlet segments for each segment in this well.
    for (const auto& segment : this->segments) {
        if (const auto outlet = segment.outlet_segment; outlet != 0) {
            auto segIxPos = segNumToIx.find(outlet);
            if (segIxPos == segNumToIx.end()) {
                continue;
            }

            this->segments[segIxPos->second]
                .inflow_segments.push_back(segment.segment);
        }
    }
}


const Opm::RestartIO::RstSegment&
Opm::RestartIO::RstWell::segment(int segment_number) const
{
    auto iter = std::find_if(this->segments.begin(), this->segments.end(),
                             [segment_number](const RstSegment& segment)
                             { return segment.segment == segment_number; });

    if (iter == this->segments.end()) {
        throw std::invalid_argument("No such segment");
    }

    return *iter;
}
