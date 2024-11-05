/*
  Copyright 2021 Equinor ASA.

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

#include <opm/input/eclipse/Schedule/ScheduleState.hpp>

#include <opm/input/eclipse/Parser/ParserKeywords/F.hpp>

#include <opm/input/eclipse/Schedule/Action/Actions.hpp>
#include <opm/input/eclipse/Schedule/GasLiftOpt.hpp>
#include <opm/input/eclipse/Schedule/Group/GConSale.hpp>
#include <opm/input/eclipse/Schedule/Group/GConSump.hpp>
#include <opm/input/eclipse/Schedule/Group/GSatProd.hpp>
#include <opm/input/eclipse/Schedule/Group/GroupEconProductionLimits.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRateConfig.hpp>
#include <opm/input/eclipse/Schedule/Network/Balance.hpp>
#include <opm/input/eclipse/Schedule/Network/ExtNetwork.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/ReservoirCouplingInfo.hpp>
#include <opm/input/eclipse/Schedule/RFTConfig.hpp>
#include <opm/input/eclipse/Schedule/RPTConfig.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQActive.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQConfig.hpp>
#include <opm/input/eclipse/Schedule/VFPInjTable.hpp>
#include <opm/input/eclipse/Schedule/VFPProdTable.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellMatcher.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTestConfig.hpp>

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <ctime>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <fmt/format.h>

namespace {

    std::pair<std::size_t, std::size_t>
    date_diff(const Opm::time_point& t2, const Opm::time_point& t1)
    {
        const auto ts1 = Opm::TimeStampUTC { Opm::TimeService::to_time_t(t1) };
        const auto ts2 = Opm::TimeStampUTC { Opm::TimeService::to_time_t(t2) };

        const auto year_diff  = ts2.year() - ts1.year();
        const auto month_diff = year_diff*12 + ts2.month() - ts1.month();

        return { year_diff, month_diff };
    }

} // Anonymous namespace

namespace Opm {

void ScheduleState::updateSAVE(bool save) {
    this->m_save_step = save;
}

bool ScheduleState::save() const {
    return this->m_save_step;
}

ScheduleState::ScheduleState(const time_point& t1)
    : m_start_time(t1)
    , m_first_in_month(true)
    , m_first_in_year(true)
{
    auto ts1 = TimeStampUTC(TimeService::to_time_t(this->m_start_time));
    this->m_month_num = ts1.month() - 1;
}

ScheduleState::ScheduleState(const time_point& start_time, const time_point& end_time)
    : ScheduleState(start_time)
{
    this->m_end_time = end_time;
}

void ScheduleState::update_date(const time_point& prev_time)
{
    auto [year_diff, month_diff] = date_diff(this->m_start_time, prev_time);
    this->m_year_num += year_diff;
    this->m_first_in_month = (month_diff > 0);
    this->m_first_in_year = (year_diff > 0);

    auto ts1 = TimeStampUTC(TimeService::to_time_t(this->m_start_time));
    this->m_month_num = ts1.month() - 1;
}

ScheduleState::ScheduleState(const ScheduleState& src, const time_point& start_time)
    : ScheduleState { src }     // Copy constructor
{
    this->m_start_time = start_time;
    this->m_end_time = std::nullopt;
    this->m_sim_step = src.sim_step() + 1;
    this->m_events.reset();
    this->m_wellgroup_events.reset();
    this->m_geo_keywords.clear();
    this->target_wellpi.clear();
    this->m_save_step = false;

    {
        auto next_rft = this->rft_config().next();

        if (next_rft.has_value()) {
            this->rft_config.update(std::move(*next_rft));
        }
    }

    {
        this->update_date(src.m_start_time);

        if (this->rst_config().save) {
            auto new_rst = this->rst_config();
            new_rst.save = false;
            this->rst_config.update(std::move(new_rst));
        }
    }

    if (this->next_tstep.has_value()) {
        if (!this->next_tstep->every_report()) {
            this->next_tstep = std::nullopt;
        }
        // Need to signal an event also for the persistance to take effect
        this->events().addEvent(ScheduleEvents::TUNING_CHANGE);
    }

    // TSINIT from TUNING should only apply to one report step, but TUNING
    // was copied from last ScheduleState. If that has TSINIT set then
    // the first time step would be limited if a TUNING_CHANGE event happens
    // e.g. because of above or because of NEXTSTEP in ACTIONX
    this->m_tuning.TSINIT.reset();

    {
        auto new_udq = this->udq();

        if (new_udq.clear_pending_assignments()) {
            // New report step.  All ASSIGNments from previous report steps
            // have been performed.
            this->udq.update(std::move(new_udq));
        }
    }
}

ScheduleState::ScheduleState(const ScheduleState& src,
                             const time_point& start_time,
                             const time_point& end_time)
    : ScheduleState { src, start_time }
{
    this->m_end_time = end_time;
}

time_point ScheduleState::start_time() const {
    return this->m_start_time;
}

time_point ScheduleState::end_time() const {
    return this->m_end_time.value();
}

std::size_t ScheduleState::sim_step() const {
    return this->m_sim_step;
}

std::size_t ScheduleState::month_num() const {
    return this->m_month_num;
}

std::size_t ScheduleState::year_num() const {
    return this->m_year_num;
}

bool ScheduleState::first_in_month() const {
    return this->m_first_in_month;
}

bool ScheduleState::first_in_year() const {
    return this->m_first_in_year;
}

void ScheduleState::init_nupcol(Nupcol nupcol) {
    this->m_nupcol = std::move(nupcol);
}

void ScheduleState::update_nupcol(int nupcol) {
    this->m_nupcol.update(nupcol);
}

int ScheduleState::nupcol() const {
    return this->m_nupcol.value();
}

void ScheduleState::update_oilvap(OilVaporizationProperties oilvap) {
    this->m_oilvap = std::move(oilvap);
}

const OilVaporizationProperties& ScheduleState::oilvap() const {
    return this->m_oilvap;
}

OilVaporizationProperties& ScheduleState::oilvap() {
    return this->m_oilvap;
}

void ScheduleState::update_geo_keywords(std::vector<DeckKeyword> geo_keywords) {
    this->m_geo_keywords = std::move(geo_keywords);
}

std::vector<DeckKeyword>& ScheduleState::geo_keywords() {
    return this->m_geo_keywords;
}

const std::vector<DeckKeyword>& ScheduleState::geo_keywords() const {
    return this->m_geo_keywords;
}

void ScheduleState::update_message_limits(MessageLimits message_limits) {
    this->m_message_limits = std::move(message_limits);
}

const MessageLimits& ScheduleState::message_limits() const {
    return this->m_message_limits;
}

MessageLimits& ScheduleState::message_limits() {
    return this->m_message_limits;
}

Well::ProducerCMode ScheduleState::whistctl() const {
    return this->m_whistctl_mode;
}

void ScheduleState::update_whistctl(Well::ProducerCMode whistctl) {
    this->m_whistctl_mode = whistctl;
}

const std::optional<double>& ScheduleState::sumthin() const {
    return this->m_sumthin;
}

void ScheduleState::update_sumthin(double sumthin) {
    if (! (sumthin > 0.0))
        this->m_sumthin.reset();
    else
        this->m_sumthin = sumthin;
}

bool ScheduleState::rptonly() const
{
    return this->m_rptonly;
}

void ScheduleState::rptonly(const bool only)
{
    this->m_rptonly = only;
}

bool ScheduleState::operator==(const ScheduleState& other) const {

    return this->m_start_time == other.m_start_time
        && this->m_oilvap == other.m_oilvap
        && this->m_sim_step == other.m_sim_step
        && this->m_month_num == other.m_month_num
        && this->m_save_step == other.m_save_step
        && this->m_first_in_month == other.m_first_in_month
        && this->m_first_in_year == other.m_first_in_year
        && this->m_year_num == other.m_year_num
        && this->target_wellpi == other.target_wellpi
        && this->m_tuning == other.m_tuning
        && this->m_end_time == other.m_end_time
        && this->m_events == other.m_events
        && this->m_wellgroup_events == other.m_wellgroup_events
        && this->m_geo_keywords == other.m_geo_keywords
        && this->m_message_limits == other.m_message_limits
        && this->m_whistctl_mode == other.m_whistctl_mode
        && this->m_nupcol == other.m_nupcol
        && this->network.get() == other.network.get()
        && this->network_balance.get() == other.network_balance.get()
        && this->wtest_config.get() == other.wtest_config.get()
        && this->well_order.get() == other.well_order.get()
        && this->group_order.get() == other.group_order.get()
        && this->gconsale.get() == other.gconsale.get()
        && this->gconsump.get() == other.gconsump.get()
        && this->gsatprod.get() == other.gsatprod.get()
        && this->wlist_manager.get() == other.wlist_manager.get()
        && this->rpt_config.get() == other.rpt_config.get()
        && this->actions.get() == other.actions.get()
        && this->udq_active.get() == other.udq_active.get()
        && this->glo.get() == other.glo.get()
        && this->guide_rate.get() == other.guide_rate.get()
        && this->rft_config.get() == other.rft_config.get()
        && this->udq.get() == other.udq.get()
        && this->bhp_defaults.get() == other.bhp_defaults.get()
        && this->source.get() == other.source.get()
        && this->wells == other.wells
        && this->inj_streams == other.inj_streams
        && this->groups == other.groups
        && this->vfpprod == other.vfpprod
        && this->vfpinj == other.vfpinj
        && this->m_sumthin == other.m_sumthin
        && this->next_tstep == other.next_tstep
        && this->m_rptonly == other.m_rptonly
        ;
}



ScheduleState ScheduleState::serializationTestObject() {
    auto t1 = TimeService::now();
    auto t2 = t1 + std::chrono::hours(48);
    ScheduleState ts(t1, t2);
    ts.m_sim_step = 123;
    ts.m_month_num = 12;
    ts.m_year_num = 66;
    ts.vfpprod = map_member<int, VFPProdTable>::serializationTestObject();
    ts.vfpinj = map_member<int, VFPInjTable>::serializationTestObject();
    ts.groups = map_member<std::string, Group>::serializationTestObject();
    ts.m_events = Events::serializationTestObject();
    ts.m_nupcol = Nupcol::serializationTestObject();
    ts.update_oilvap( Opm::OilVaporizationProperties::serializationTestObject() );
    ts.m_message_limits = MessageLimits::serializationTestObject();
    ts.m_whistctl_mode = Well::ProducerCMode::THP;
    ts.target_wellpi = {{"WELL1", 1000}, {"WELL2", 2000}};

    ts.m_sumthin = 12.345;
    ts.m_rptonly = true;

    ts.bhp_defaults.update( BHPDefaults::serializationTestObject() );
    ts.pavg.update( PAvg::serializationTestObject() );
    ts.wtest_config.update( WellTestConfig::serializationTestObject() );
    ts.gconsump.update( GConSump::serializationTestObject() );
    ts.gconsale.update( GConSale::serializationTestObject() );
    ts.gsatprod.update( GSatProd::serializationTestObject() );
    ts.gecon.update( GroupEconProductionLimits::serializationTestObject() );
    ts.rescoup.update( ReservoirCoupling::CouplingInfo::serializationTestObject() );
    ts.wlist_manager.update( WListManager::serializationTestObject() );
    ts.rpt_config.update( RPTConfig::serializationTestObject() );
    ts.actions.update( Action::Actions::serializationTestObject() );
    ts.udq_active.update( UDQActive::serializationTestObject() );
    ts.network.update( Network::ExtNetwork::serializationTestObject() );
    ts.network_balance.update( Network::Balance::serializationTestObject() );
    ts.well_order.update( NameOrder::serializationTestObject() );
    ts.group_order.update( GroupOrder::serializationTestObject() );
    ts.udq.update( UDQConfig::serializationTestObject() );
    ts.guide_rate.update( GuideRateConfig::serializationTestObject() );
    ts.glo.update( GasLiftOpt::serializationTestObject() );
    ts.rft_config.update( RFTConfig::serializationTestObject() );
    ts.rst_config.update( RSTConfig::serializationTestObject() );
    ts.source.update( Source::serializationTestObject() );

    return ts;
}

void ScheduleState::update_tuning(Tuning tuning) {
    this->m_tuning = std::move(tuning);
}

const Tuning& ScheduleState::tuning() const {
    return this->m_tuning;
}

Tuning& ScheduleState::tuning() {
    return this->m_tuning;
}

// Returns -1 if there is no active limit on next step (from TUNING or NEXT[STEP])
double ScheduleState::max_next_tstep(const bool enableTUNING) const {
    double tuning_value = (enableTUNING && this->m_tuning.TSINIT.has_value())  ? this->m_tuning.TSINIT.value() : -1.0;
    double next_value = this->next_tstep.has_value() ? this->next_tstep->value() : -1.0;
    return std::max(next_value, tuning_value);
}

void ScheduleState::update_events(Events events) {
    this->m_events = events;
}

Events& ScheduleState::events() {
    return this->m_events;
}


const Events& ScheduleState::events() const {
    return this->m_events;
}

void ScheduleState::update_wellgroup_events(WellGroupEvents wgevents) {
    this->m_wellgroup_events = std::move(wgevents);
}

WellGroupEvents& ScheduleState::wellgroup_events() {
    return this->m_wellgroup_events;
}

const WellGroupEvents& ScheduleState::wellgroup_events() const {
    return this->m_wellgroup_events;
}


/*
  Observe that the decision to write a restart file will typically be a
  combination of the RST configuration from the previous report step, and the
  first_in_year++ attributes of this report step. That is the reason the
  function takes a RSTConfig argument - instead of using the rst_config member.

*/

bool ScheduleState::rst_file(const RSTConfig&  rst,
                             const time_point& previous_restart_output_time) const
{
    if (rst.save)
        return true;

    if (rst.write_rst_file.has_value())
        return rst.write_rst_file.value();

    const auto freq = rst.freq.value_or(1);
    const auto basic = rst.basic.value_or(0);

    if (basic == 0)
        return false;

    if (basic == 3)
        return (this->sim_step() % freq) == 0;

    const auto [year_diff, month_diff] =
        date_diff(this->m_start_time, previous_restart_output_time);

    if (basic == 4) {
        return this->first_in_year()
            && (year_diff >= static_cast<std::size_t>(freq));
    }

    if (basic == 5) {
        return this->first_in_month()
            && (month_diff >= static_cast<std::size_t>(freq));
    }

    throw std::logic_error(fmt::format("Unsupported BASIC={} value", basic));
}


bool ScheduleState::has_gpmaint() const
{
    return std::any_of(this->groups.begin(), this->groups.end(), [](const auto& name_group) {
            return name_group.second->gpmaint().has_value();
    });
}


} // namespace Opm
