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

#ifndef RST_WELL
#define RST_WELL

#include <opm/io/eclipse/rst/connection.hpp>
#include <opm/io/eclipse/rst/segment.hpp>

#include <opm/input/eclipse/Schedule/ScheduleTypes.hpp>

#include <array>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace Opm {
class UnitSystem;
} // namespace Opm

namespace Opm { namespace RestartIO {

struct RstHeader;

struct RstWell
{
    RstWell(const ::Opm::UnitSystem& unit_system,
            const RstHeader& header,
            const std::string& group_arg,
            const std::string* zwel,
            const int * iwel,
            const float * swel,
            const double * xwel,
            const int * icon,
            const float * scon,
            const double * xcon);

    RstWell(const ::Opm::UnitSystem& unit_system,
            const RstHeader& header,
            const std::string& group_arg,
            const std::string* zwel,
            const int * iwel,
            const float * swel,
            const double * xwel,
            const int * icon,
            const float * scon,
            const double * xcon,
            const std::vector<int>& iseg,
            const std::vector<double>& rseg);

    std::string name;
    std::string group;
    std::array<int, 2> ij;
    std::pair<int,int> k1k2;
    WellType wtype;
    int well_status;
    int active_control;
    int vfp_table;
    int econ_workover_procedure;
    int preferred_phase;
    bool allow_xflow;
    int group_controllable_flag;
    int econ_limit_end_run;
    int grupcon_gr_phase;
    int hist_requested_control;
    int msw_index;
    int completion_ordering;
    int pvt_table;
    int msw_pressure_drop_model;
    int wtest_config_reasons;
    int wtest_close_reason;
    int wtest_remaining;
    int econ_limit_quantity;
    int econ_workover_procedure_2;
    int thp_lookup_procedure_vfptable;
    int close_if_thp_stabilised;
    int prevent_thpctrl_if_unstable;
    bool glift_active;
    bool glift_alloc_extra_gas;

    float orat_target;
    float wrat_target;
    float grat_target;
    float lrat_target;
    float resv_target;
    float thp_target;
    float bhp_target_float;
    float vfp_bhp_adjustment;
    float vfp_bhp_scaling_factor;
    float hist_lrat_target;
    float hist_grat_target;
    float hist_bhp_target;
    float datum_depth;
    float drainage_radius;
    float grupcon_gr_value;
    float efficiency_factor;
    float alq_value;
    float econ_limit_min_oil;
    float econ_limit_min_gas;
    float econ_limit_max_wct;
    float econ_limit_max_gor;
    float econ_limit_max_wgr;
    float econ_limit_max_wct_2;
    float econ_limit_min_liq;
    float wtest_interval;
    float wtest_startup;
    float grupcon_gr_scaling;
    float glift_max_rate;
    float glift_min_rate;
    float glift_weight_factor;
    float glift_inc_weight_factor;
    float dfac_corr_coeff_a{};
    float dfac_corr_exponent_b{};
    float dfac_corr_exponent_c{};
    float inj_temperature{UNDEFINED_VALUE};
    std::vector<float> tracer_concentration_injection;

    double oil_rate;
    double water_rate;
    double gas_rate;
    double liquid_rate;
    double void_rate;
    double thp;
    double flow_bhp;
    double wct;
    double gor;
    double oil_total;
    double water_total;
    double gas_total;
    double void_total;
    double water_inj_total;
    double gas_inj_total;
    double void_inj_total;
    double gas_fvf;
    double bhp_target_double;
    double hist_oil_total;
    double hist_wat_total;
    double hist_gas_total;
    double hist_water_inj_total;
    double hist_gas_inj_total;
    double water_void_rate;
    double gas_void_rate;

    static constexpr auto UNDEFINED_VALUE = 1.0e20f;

    const RstSegment& segment(int segment_number) const;
    std::vector<RstConnection> connections;
    std::vector<RstSegment> segments;
};

}} // namespace Opm::RestartIO

#endif // RST_WELL
