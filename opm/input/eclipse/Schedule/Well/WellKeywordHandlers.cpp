/*
  Copyright 2020 Statoil ASA.

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

#include "WellKeywordHandlers.hpp"
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/utility/OpmInputError.hpp>
#include <opm/common/utility/String.hpp>

#include <opm/input/eclipse/Deck/DeckKeyword.hpp>

#include <opm/input/eclipse/Parser/ParseContext.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/F.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/W.hpp>

#include <opm/input/eclipse/Schedule/Network/ExtNetwork.hpp>
#include <opm/input/eclipse/Schedule/ScheduleState.hpp>
#include <opm/input/eclipse/Schedule/ScheduleStatic.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQActive.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQConfig.hpp>
#include <opm/input/eclipse/Schedule/Well/NameOrder.hpp>
#include <opm/input/eclipse/Schedule/Well/WDFAC.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>
#include <opm/input/eclipse/Schedule/Well/WellEconProductionLimits.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTestConfig.hpp>
#include <opm/input/eclipse/Schedule/Well/WListManager.hpp>
#include <opm/input/eclipse/Schedule/Well/WVFPDP.hpp>
#include <opm/input/eclipse/Schedule/Well/WVFPEXP.hpp>

#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include "../HandlerContext.hpp"

#include <fmt/format.h>
#include <fmt/ranges.h>

#include <algorithm>
#include <memory>

namespace Opm {

namespace {

/*
  The function trim_wgname() is used to trim the leading and trailing spaces
  away from the group and well arguments given in the WELSPECS
  keyword. If the deck argument contains a leading or trailing space that is
  treated as an input error, and the action taken is regulated by the setting
  ParseContext::PARSE_WGNAME_SPACE.

  Observe that the spaces are trimmed *unconditionally* - i.e. if the
  ParseContext::PARSE_WGNAME_SPACE setting is set to InputError::IGNORE that
  means that we do not inform the user about "our fix", but it is *not* possible
  to configure the parser to leave the spaces intact.
*/
std::string trim_wgname(const DeckKeyword& keyword,
                        const std::string& wgname_arg,
                        const ParseContext& parseContext,
                        ErrorGuard& errors)
{
    std::string wgname = trim_copy(wgname_arg);
    if (wgname != wgname_arg)  {
        const auto& location = keyword.location();
        std::string msg_fmt = fmt::format("Problem with keyword {{keyword}}\n"
                                          "In {{file}} line {{line}}\n"
                                          "Illegal space in {} when defining WELL/GROUP.", wgname_arg);
        parseContext.handleError(ParseContext::PARSE_WGNAME_SPACE, msg_fmt, location, errors);
    }
    return wgname;
}

void handleWCONHIST(HandlerContext& handlerContext)
{
    for (const auto& record : handlerContext.keyword) {
        const std::string& wellNamePattern = record.getItem("WELL").getTrimmedString(0);
        const auto well_names = handlerContext.wellNames(wellNamePattern, false);

        const Well::Status status = WellStatusFromString(record.getItem("STATUS").getTrimmedString(0));

        for (const auto& well_name : well_names) {
            handlerContext.updateWellStatus(well_name, status,
                                            handlerContext.keyword.location());

            std::optional<VFPProdTable::ALQ_TYPE> alq_type;
            auto well2 = handlerContext.state().wells.get( well_name );
            const bool switching_from_injector = !well2.isProducer();
            auto properties = std::make_shared<Well::WellProductionProperties>(well2.getProductionProperties());
            bool update_well = false;

            auto table_nr = record.getItem("VFP_TABLE").get< int >(0);
            if (record.getItem("VFP_TABLE").defaultApplied(0)) { // Default 1* use the privious set vfp table
                table_nr = properties->VFPTableNumber;
            }

            if (table_nr != 0) {
                const auto& vfpprod = handlerContext.state().vfpprod;
                if (vfpprod.has(table_nr)) {
                    alq_type = handlerContext.state().vfpprod(table_nr).getALQType();
                } else {
                    std::string reason = fmt::format("Problem with well:{} VFP table: {} not defined", well_name, table_nr);
                    throw OpmInputError(reason, handlerContext.keyword.location());
                }
            }
            double default_bhp;
            if (handlerContext.state().bhp_defaults.get().prod_target) {
                default_bhp = *handlerContext.state().bhp_defaults.get().prod_target;
            } else {
                default_bhp = UnitSystem::newMETRIC().to_si(UnitSystem::measure::pressure,
                                                            ParserKeywords::FBHPDEF::TARGET_BHP::defaultValue);
            }

            // Injectors at a restart time will not have any WellProductionProperties with the
            // proper whistctl_cmode, so this needs to be set before the call to handleWCONHIST
            if (switching_from_injector) {
                properties->whistctl_cmode = handlerContext.state().whistctl();
            }

            properties->handleWCONHIST(alq_type,
                                       table_nr,
                                       default_bhp,
                                       handlerContext.static_schedule().m_unit_system, record);

            if (switching_from_injector) {
                if (properties->bhp_hist_limit_defaulted) {
                    properties->setBHPLimit(default_bhp);
                }

                auto inj_props = std::make_shared<Well::WellInjectionProperties>(well2.getInjectionProperties());
                inj_props->resetBHPLimit();
                well2.updateInjection(inj_props);
                update_well = true;
                handlerContext.state().wellgroup_events().addEvent( well2.name(), ScheduleEvents::WELL_SWITCHED_INJECTOR_PRODUCER);
            }

            if (well2.updateProduction(properties)) {
                update_well = true;
            }

            if (well2.updatePrediction(false)) {
                update_well = true;
            }

            if (well2.updateHasProduced()) {
                update_well = true;
            }

            if (update_well) {
                handlerContext.state().events().addEvent( ScheduleEvents::PRODUCTION_UPDATE );
                handlerContext.state().wellgroup_events().addEvent( well2.name(), ScheduleEvents::PRODUCTION_UPDATE);
                handlerContext.state().wells.update( well2 );
            }

            // Always check if well can be opened (it could have been closed for numerical reasons and possible to operate with new params)
            if (handlerContext.getWellStatus(well_name) == WellStatus::OPEN) {
                handlerContext.state().wellgroup_events().addEvent( well_name, ScheduleEvents::REQUEST_OPEN_WELL);
            }

        }
    }
}

void handleWCONINJE(HandlerContext& handlerContext)
{
    for (const auto& record : handlerContext.keyword) {
        const std::string& wellNamePattern = record.getItem("WELL").getTrimmedString(0);
        const auto well_names = handlerContext.wellNames(wellNamePattern);

        const Well::Status status = WellStatusFromString(record.getItem("STATUS").getTrimmedString(0));

        for (const auto& well_name : well_names) {
            handlerContext.updateWellStatus(well_name, status,
                                            handlerContext.keyword.location());

            bool update_well = false;
            auto well2 = handlerContext.state().wells.get( well_name );

            auto injection = std::make_shared<Well::WellInjectionProperties>(well2.getInjectionProperties());
            auto previousInjectorType = injection->injectorType;

            double default_bhp_limit;
            const auto& usys = handlerContext.static_schedule().m_unit_system;
            if (handlerContext.state().bhp_defaults.get().inj_limit) {
                default_bhp_limit = usys.from_si(UnitSystem::measure::pressure,
                                                 *handlerContext.state().bhp_defaults.get().inj_limit);
            } else {
                default_bhp_limit = UnitSystem::newMETRIC().to_si(UnitSystem::measure::pressure,
                                                                  ParserKeywords::WCONINJE::BHP::defaultValue.get<double>());
                default_bhp_limit = usys.from_si(UnitSystem::measure::pressure,
                                                 default_bhp_limit);
            }

            auto table_nr = record.getItem("VFP_TABLE").get< int >(0);

            if (table_nr != 0) {
                const auto& vfpinj = handlerContext.state().vfpinj;
                if (!vfpinj.has(table_nr)) {
                    std::string reason = fmt::format("Problem with well:{} VFP table: {} not defined", well_name, table_nr);
                    throw OpmInputError(reason, handlerContext.keyword.location());
                }
            }

            injection->handleWCONINJE(record, default_bhp_limit,
                                      well2.isAvailableForGroupControl(), well_name, handlerContext.keyword.location());

            const bool switching_from_producer = well2.isProducer();
            if (well2.updateInjection(injection)) {
                update_well = true;
            }

            if (switching_from_producer) {
                handlerContext.state().wellgroup_events().addEvent( well2.name(), ScheduleEvents::WELL_SWITCHED_INJECTOR_PRODUCER);
            }

            if (well2.updatePrediction(true)) {
                update_well = true;
            }

            if (well2.updateHasInjected()) {
                update_well = true;
            }

            if (update_well) {
                handlerContext.state().events().addEvent(ScheduleEvents::INJECTION_UPDATE);
                handlerContext.state().wellgroup_events().addEvent( well_name, ScheduleEvents::INJECTION_UPDATE);
                if (previousInjectorType != injection->injectorType) {
                    handlerContext.state().wellgroup_events().addEvent( well_name, ScheduleEvents::INJECTION_TYPE_CHANGED);
                }
                handlerContext.state().wells.update( std::move(well2) );
            }

            if (handlerContext.state().wells.get( well_name ).getStatus() == Well::Status::OPEN) {
                handlerContext.state().wellgroup_events().addEvent(well_name, ScheduleEvents::REQUEST_OPEN_WELL);
            }

            auto udq_active = handlerContext.state().udq_active.get();
            if (injection->updateUDQActive(handlerContext.state().udq.get(), udq_active)) {
                handlerContext.state().udq_active.update( std::move(udq_active) );
            }

            handlerContext.affected_well(well_name);
        }
    }
}

void handleWCONINJH(HandlerContext& handlerContext)
{
    for (const auto& record : handlerContext.keyword) {
        const std::string& wellNamePattern = record.getItem("WELL").getTrimmedString(0);
        const auto well_names = handlerContext.wellNames(wellNamePattern, false);
        const Well::Status status = WellStatusFromString( record.getItem("STATUS").getTrimmedString(0));

        for (const auto& well_name : well_names) {
            handlerContext.updateWellStatus(well_name, status,
                                            handlerContext.keyword.location());
            bool update_well = false;
            auto well2 = handlerContext.state().wells.get( well_name );
            auto injection = std::make_shared<Well::WellInjectionProperties>(well2.getInjectionProperties());
            auto previousInjectorType = injection->injectorType;

            double default_bhp_limit;
            if (handlerContext.state().bhp_defaults.get().inj_limit) {
                default_bhp_limit = *handlerContext.state().bhp_defaults.get().inj_limit;
            } else {
                default_bhp_limit = UnitSystem::newMETRIC().to_si(UnitSystem::measure::pressure,
                                                                  6891.2);
            }

            auto table_nr = record.getItem("VFP_TABLE").get< int >(0);
            if (record.getItem("VFP_TABLE").defaultApplied(0)) { // Default 1* use the privious set vfp table
                table_nr = injection->VFPTableNumber;
            }
            if (table_nr != 0) {
                const auto& vfpinj = handlerContext.state().vfpinj;
                if (!vfpinj.has(table_nr)) {
                    std::string reason = fmt::format("Problem with well:{} VFP table: {} not defined", well_name, table_nr);
                    throw OpmInputError(reason, handlerContext.keyword.location());
                }
            }

            injection->handleWCONINJH(record, table_nr, default_bhp_limit,
                                      well2.isProducer(), well_name,
                                      handlerContext.keyword.location());

            const bool switching_from_producer = well2.isProducer();
            if (well2.updateInjection(injection)) {
                update_well = true;
            }

            if (switching_from_producer) {
                handlerContext.state().wellgroup_events().addEvent( well2.name(), ScheduleEvents::WELL_SWITCHED_INJECTOR_PRODUCER);
            }

            if (well2.updatePrediction(false)) {
                update_well = true;
            }

            if (well2.updateHasInjected()) {
                update_well = true;
            }

            if (update_well) {
                handlerContext.state().events().addEvent( ScheduleEvents::INJECTION_UPDATE );
                handlerContext.state().wellgroup_events().addEvent( well_name, ScheduleEvents::INJECTION_UPDATE);
                if (previousInjectorType != injection->injectorType) {
                    handlerContext.state().wellgroup_events().addEvent( well_name, ScheduleEvents::INJECTION_TYPE_CHANGED);
                }
                handlerContext.state().wells.update( std::move(well2) );
            }

            // Always check if well can be opened (it could have been closed for numerical reasons and possible to operate with new params)
            if (handlerContext.getWellStatus(well_name) == WellStatus::OPEN) {
                handlerContext.state().wellgroup_events().addEvent( well_name, ScheduleEvents::REQUEST_OPEN_WELL);
            }

        }
    }
}

bool belongsToAutoChokeGroup(const Well& well, const ScheduleState& state) {
    const auto& network = state.network.get();
    if (!network.active())
        return false;
    auto group_name = well.groupName();
    while (group_name != "FIELD") {
        if (network.has_node(group_name)) {
            auto node_name = group_name;
            if (network.node(node_name).as_choke())
                return true;
            while (network.uptree_branch(node_name)) {
                node_name = network.uptree_branch(node_name)->uptree_node();
                if (network.node(node_name).as_choke())
                    return true;
            }
        }
        group_name = state.groups.get(group_name).parent();
    }
    return false;
}

void handleWCONPROD(HandlerContext& handlerContext)
{
    for (const auto& record : handlerContext.keyword) {
        const std::string& wellNamePattern = record.getItem("WELL").getTrimmedString(0);
        const auto well_names = handlerContext.wellNames(wellNamePattern, false);

        const Well::Status status = WellStatusFromString(record.getItem("STATUS").getTrimmedString(0));

        for (const auto& well_name : well_names) {
            bool update_well = handlerContext.updateWellStatus(well_name, status,
                                                               handlerContext.keyword.location());
            std::optional<VFPProdTable::ALQ_TYPE> alq_type;
            auto well2 = handlerContext.state().wells.get( well_name );
            const bool switching_from_injector = !well2.isProducer();
            auto properties = std::make_shared<Well::WellProductionProperties>(well2.getProductionProperties());
            properties->clearControls();
            if (well2.isAvailableForGroupControl() || belongsToAutoChokeGroup(well2, handlerContext.state())) {
                properties->addProductionControl(Well::ProducerCMode::GRUP);
            }

            auto table_nr = record.getItem("VFP_TABLE").get< int >(0);
            if (table_nr != 0) {
                const auto& vfpprod = handlerContext.state().vfpprod;
                if (vfpprod.has(table_nr)) {
                    alq_type = handlerContext.state().vfpprod(table_nr).getALQType();
                } else {
                    std::string reason = fmt::format("Problem with well:{} VFP table: {} not defined", well_name, table_nr);
                    throw OpmInputError(reason, handlerContext.keyword.location());
                }
            }

            double default_bhp_target;
            if (handlerContext.state().bhp_defaults.get().prod_target) {
                default_bhp_target = *handlerContext.state().bhp_defaults.get().prod_target;
            } else {
                default_bhp_target = UnitSystem::newMETRIC().to_si(UnitSystem::measure::pressure,
                                                                   ParserKeywords::WCONPROD::BHP::defaultValue.get<double>());
            }

            properties->handleWCONPROD(alq_type, table_nr, default_bhp_target,
                                       handlerContext.static_schedule().m_unit_system,
                                       well_name, record, handlerContext.keyword.location());

            if (switching_from_injector) {
                if (properties->bhp_hist_limit_defaulted) {
                    properties->setBHPLimit(default_bhp_target);
                }
                update_well = true;
                handlerContext.state().wellgroup_events().addEvent( well2.name(), ScheduleEvents::WELL_SWITCHED_INJECTOR_PRODUCER);
            }

            if (well2.updateProduction(properties)) {
                update_well = true;
            }

            if (well2.updatePrediction(true)) {
                update_well = true;
            }

            if (well2.updateHasProduced()) {
                update_well = true;
            }

            if (well2.getStatus() == WellStatus::OPEN) {
                handlerContext.state().wellgroup_events().addEvent(well2.name(), ScheduleEvents::REQUEST_OPEN_WELL);
            }

            if (update_well) {
                handlerContext.state().events().addEvent( ScheduleEvents::PRODUCTION_UPDATE );
                handlerContext.state().wellgroup_events().addEvent( well2.name(), ScheduleEvents::PRODUCTION_UPDATE);
                handlerContext.state().wells.update( std::move(well2) );
            }

            auto udq_active = handlerContext.state().udq_active.get();
            if (properties->updateUDQActive(handlerContext.state().udq.get(), udq_active)) {
                handlerContext.state().udq_active.update( std::move(udq_active));
            }

            handlerContext.affected_well(well_name);
        }
    }
}

void handleWCYCLE(HandlerContext& handlerContext)
{
    auto new_config = handlerContext.state().wcycle();
    for (const auto& record : handlerContext.keyword) {
        new_config.addRecord(record);
    }
    handlerContext.state().wcycle.update(std::move(new_config));
}

void handleWELLSTRE(HandlerContext& handlerContext)
{
    auto& inj_streams = handlerContext.state().inj_streams;
    for (const auto& record : handlerContext.keyword) {
        const std::string& stream_name = record.getItem("STREAM").getTrimmedString(0);
        const auto& compostion = record.getItem("COMPOSITIONS").getData<double>();
        // TODO: we need to check the number of the values in the composition is the same as the number of components

        const double sum = std::accumulate(compostion.begin(), compostion.end(), 0.0);
        if (std::abs(sum - 1.0) > std::numeric_limits<double>::epsilon()) {
            std::string msg = fmt::format("The sum of the composition values for stream '{}' is not 1.0, but {}.", stream_name, sum);
            throw OpmInputError(msg, handlerContext.keyword.location());
        }
        auto composition_ptr = std::make_shared<std::vector<double>>(std::move(compostion));
        inj_streams.update(stream_name, composition_ptr);
    }

}

void handleWELOPEN(HandlerContext& handlerContext)
{
    const auto& keyword = handlerContext.keyword;

    auto conn_defaulted = []( const DeckRecord& rec ) {
        auto defaulted = []( const DeckItem& item ) {
            return item.defaultApplied( 0 );
        };

        return std::all_of( rec.begin() + 2, rec.end(), defaulted );
    };

    constexpr auto open = Well::Status::OPEN;

    for (const auto& record : keyword) {
        const auto& wellNamePattern = record.getItem( "WELL" ).getTrimmedString(0);
        const auto& status_str = record.getItem( "STATUS" ).getTrimmedString( 0 );
        const auto well_names = handlerContext.wellNames(wellNamePattern);

        /* if all records are defaulted or just the status is set, only
         * well status is updated
         */
        if (conn_defaulted(record)) {
            const auto new_well_status = WellStatusFromString(status_str);
            for (const auto& wname : well_names) {
                const auto did_update_well_status =
                    handlerContext.updateWellStatus(wname, new_well_status);

                handlerContext.affected_well(wname);

                if (did_update_well_status) {
                    handlerContext.record_well_structure_change();
                }

                if (did_update_well_status && (new_well_status == open)) {
                    // Record possible well injection/production status change
                    auto well2 = handlerContext.state().wells.get(wname);

                    const auto did_flow_update =
                        (well2.isProducer() && well2.updateHasProduced())
                        ||
                        (well2.isInjector() && well2.updateHasInjected());

                    if (did_flow_update) {
                        handlerContext.state().wells.update(std::move(well2));
                    }
                }

                if (new_well_status == open) {
                    handlerContext.state().wellgroup_events().addEvent( wname, ScheduleEvents::REQUEST_OPEN_WELL);
                }
            }
            continue;
        }

        /*
          Some of the connection information has been entered, in this case
          we *only* update the status of the connections, and not the well
          itself. Unless all connections are shut - then the well is also
          shut.
         */
        for (const auto& wname : well_names) {
            {
                auto well = handlerContext.state().wells.get(wname);
                handlerContext.state().wells.update( std::move(well) );
            }

            const auto connection_status = Connection::StateFromString( status_str );
            {
                auto well = handlerContext.state().wells.get(wname);
                well.handleWELOPENConnections(record, connection_status);
                handlerContext.state().wells.update( std::move(well) );
            }

            handlerContext.affected_well(wname);
            handlerContext.record_well_structure_change();

            handlerContext.state().events().addEvent(ScheduleEvents::COMPLETION_CHANGE);
        }
    }
}

void handleWINJGAS(HandlerContext& HandlerContext)
{
    for (const auto& record : HandlerContext.keyword) {
        const std::string& wellNamePattern = record.getItem("WELL").getTrimmedString(0);
        const auto well_names = HandlerContext.wellNames(wellNamePattern, false);
        const std::string& fluid_nature = record.getItem("FLUID").getTrimmedString(0);
        // TODO: technically, only the firt two characters are significant
        // TODO: we need to test it out whether only the first two characters matter
        if (fluid_nature != "STREAM") {
            std::string msg = fmt::format("The fluid nature '{}' is not supported in WINJGAS keyword.", fluid_nature);
            throw OpmInputError(msg, HandlerContext.keyword.location());
        }
        const std::string& stream_name = record.getItem("STREAM").getTrimmedString(0);
        // TODO: we do not handle other records for now

        // we make sure the stream is defined in WELLSTRE keyword
        const auto& inj_streams = HandlerContext.state().inj_streams;
        if (!inj_streams.has(stream_name)) {
            std::string msg = fmt::format("The stream '{}' is not defined in WELLSTRE keyword.", stream_name);
            throw OpmInputError(msg, HandlerContext.keyword.location());
        }

        auto well2 = HandlerContext.state().wells.get(well_names[0]);
        // if (well2.isProducer()) {
        //     std::string msg = fmt::format("The well '{}' is a producer, not an injector.", well_names[0]);
        //     throw OpmInputError(msg, HandlerContext.keyword.location());
        // }
        auto injection = std::make_shared<Well::WellInjectionProperties>(well2.getInjectionProperties());

        // TODO: should we make it a injection event?
        const auto& inj_stream = inj_streams.get(stream_name);
        injection->setGasInjComposition(inj_stream);
    }
}

void handleWELSPECS(HandlerContext& handlerContext)
{
    using Kw = ParserKeywords::WELSPECS;

    auto getTrimmedName = [&handlerContext](const auto& item)
    {
        return trim_wgname(handlerContext.keyword,
                           item.template get<std::string>(0),
                           handlerContext.parseContext,
                           handlerContext.errors);
    };

    auto fieldWells = std::vector<std::string>{};
    for (const auto& record : handlerContext.keyword) {
        if (const auto fip_region_number = record.getItem<Kw::FIP_REGION>().get<int>(0);
            fip_region_number != Kw::FIP_REGION::defaultValue)
        {
            const auto& location = handlerContext.keyword.location();
            const auto msg = fmt::format("Non-defaulted FIP region {} in WELSPECS keyword "
                                         "in file {} line {} is not supported. "
                                         "Reset to default value {}.",
                                         fip_region_number,
                                         location.filename,
                                         location.lineno,
                                         Kw::FIP_REGION::defaultValue);
            OpmLog::warning(msg);
        }

        if (const auto& density_calc_type = record.getItem<Kw::DENSITY_CALC>().get<std::string>(0);
            density_calc_type != Kw::DENSITY_CALC::defaultValue)
        {
            const auto& location = handlerContext.keyword.location();
            const auto msg = fmt::format("Non-defaulted density calculation method '{}' "
                                         "in WELSPECS keyword in file {} line {} is "
                                         "not supported. Reset to default value {}.",
                                         density_calc_type,
                                         location.filename,
                                         location.lineno,
                                         Kw::DENSITY_CALC::defaultValue);
            OpmLog::warning(msg);
        }

        const auto wellName = getTrimmedName(record.getItem<Kw::WELL>());
        const auto groupName = getTrimmedName(record.getItem<Kw::GROUP>());

        // We might get here from an ACTIONX context, or we might get
        // called on a well (list) template, to reassign certain well
        // properties--e.g, the well's controlling group--so check if
        // 'wellName' matches any existing well names through pattern
        // matching before treating the wellName as a simple well name.
        //
        // An empty list of well names is okay since that means we're
        // creating a new well in this case.
        const auto allowEmptyWellList = true;
        const auto existingWells = handlerContext.wellNames(wellName, allowEmptyWellList);

        if (groupName == "FIELD") {
            if (existingWells.empty()) {
                fieldWells.push_back(wellName);
            }
            else {
                fieldWells.insert(fieldWells.end(), existingWells.begin(), existingWells.end());
            }
        }

        if (! handlerContext.state().groups.has(groupName)) {
            handlerContext.addGroup(groupName);
        }

        if (existingWells.empty()) {
            // 'wellName' does not match any existing wells.  Create a
            // new Well object for this well.
            handlerContext.welspecsCreateNewWell(record,
                                                 wellName,
                                                 groupName);
        }
        else {
            // 'wellName' matches one or more existing wells.  Assign
            // new properties for those wells.
            handlerContext.welspecsUpdateExistingWells(record,
                                                       existingWells,
                                                       groupName);
        }
    }

    if (! fieldWells.empty()) {
        std::sort(fieldWells.begin(), fieldWells.end());
        fieldWells.erase(std::unique(fieldWells.begin(), fieldWells.end()),
                         fieldWells.end());

        const auto* plural = (fieldWells.size() == 1) ? "" : "s";

        const auto msg_fmt = fmt::format(R"(Well{0} parented directly to 'FIELD'; this is allowed but discouraged.
Well{0} entered with 'FIELD' parent group:
* {1})", plural, fmt::join(fieldWells, "\n * "));

        handlerContext.parseContext.handleError(ParseContext::SCHEDULE_WELL_IN_FIELD_GROUP,
                                                msg_fmt,
                                                handlerContext.keyword.location(),
                                                handlerContext.errors);
    }

    if (! handlerContext.keyword.empty()) {
        handlerContext.record_well_structure_change();
    }
}


void handleWELSPECL(HandlerContext& handlerContext)
{
    using Kw = ParserKeywords::WELSPECL;
    auto getTrimmedName = [&handlerContext](const auto& item)
    {
        return trim_wgname(handlerContext.keyword,
                           item.template get<std::string>(0),
                           handlerContext.parseContext,
                           handlerContext.errors);
    };
    handleWELSPECS(handlerContext);
    for (const auto& record : handlerContext.keyword) {
        const auto wellName = getTrimmedName(record.getItem<Kw::WELL>());
        const auto lgrTag = getTrimmedName(record.getItem<Kw::LGR>());
        handlerContext.state().wells.get(wellName).flag_lgr_well();
        handlerContext.state().wells.get(wellName).set_lgr_well_tag(lgrTag);
    }
}

/*
  The documentation for the WELTARG keyword says that the well
  must have been fully specified and initialized using one of the
  WCONxxxx keywords prior to modifying the well using the WELTARG
  keyword.

  The following implementation of handling the WELTARG keyword
  does not check or enforce in any way that this is done (i.e. it
  is not checked or verified that the well is initialized with any
  WCONxxxx keyword).

  Update: See the discussion following the definitions of the SI factors, due
  to a bad design we currently need the well to be specified with
  WCONPROD / WCONHIST before WELTARG is applied, if not the units for the
  rates will be wrong.
*/

void handleWELTARG(HandlerContext& handlerContext)
{
    const double SiFactorP = handlerContext.static_schedule().m_unit_system.parse("Pressure").getSIScaling();
    for (const auto& record : handlerContext.keyword) {
        const std::string& wellNamePattern = record.getItem("WELL").getTrimmedString(0);
        const auto well_names = handlerContext.wellNames(wellNamePattern);
        if (well_names.empty()) {
            handlerContext.invalidNamePattern( wellNamePattern);
        }

        const auto cmode = WellWELTARGCModeFromString(record.getItem("CMODE").getTrimmedString(0));
        const auto new_arg = record.getItem("NEW_VALUE").get<UDAValue>(0);

        for (const auto& well_name : well_names) {
            auto well2 = handlerContext.state().wells.get(well_name);
            bool update = false;
            if (well2.isProducer()) {
                auto prop = std::make_shared<Well::WellProductionProperties>(well2.getProductionProperties());
                prop->handleWELTARG(cmode, new_arg, SiFactorP);
                update = well2.updateProduction(prop);
                if (cmode == Well::WELTARGCMode::GUID) {
                    update |= well2.updateWellGuideRate(new_arg.get<double>());
                }

                auto udq_active = handlerContext.state().udq_active.get();
                if (prop->updateUDQActive(handlerContext.state().udq.get(), cmode, udq_active)) {
                    handlerContext.state().udq_active.update( std::move(udq_active));
                }
            }
            else {
                auto inj = std::make_shared<Well::WellInjectionProperties>(well2.getInjectionProperties());
                inj->handleWELTARG(cmode, new_arg, SiFactorP);
                update = well2.updateInjection(inj);
                if (cmode == Well::WELTARGCMode::GUID) {
                    update |= well2.updateWellGuideRate(new_arg.get<double>());
                }

                auto udq_active = handlerContext.state().udq_active.get();
                if (inj->updateUDQActive(handlerContext.state().udq.get(), cmode, udq_active)) {
                    handlerContext.state().udq_active.update(std::move(udq_active));
                }
            }

            if (update) {
                if (well2.isProducer()) {
                    handlerContext.state().wellgroup_events().addEvent( well_name, ScheduleEvents::PRODUCTION_UPDATE);
                    handlerContext.state().events().addEvent( ScheduleEvents::PRODUCTION_UPDATE );
                } else {
                    handlerContext.state().wellgroup_events().addEvent( well_name, ScheduleEvents::INJECTION_UPDATE);
                    handlerContext.state().events().addEvent( ScheduleEvents::INJECTION_UPDATE );
                }
                handlerContext.state().wells.update( std::move(well2) );
            }

            handlerContext.affected_well(well_name);
        }
    }
}

void handleWELTRAJ(HandlerContext& handlerContext)
{
    for (const auto& record : handlerContext.keyword) {
        const std::string& wellNamePattern = record.getItem("WELL").getTrimmedString(0);
        const auto wellnames = handlerContext.wellNames(wellNamePattern, false);

        for (const auto& name : wellnames) {
            auto well2 = handlerContext.state().wells.get(name);
            auto connections = std::make_shared<WellConnections>(WellConnections(well2.getConnections()));
            connections->loadWELTRAJ(record, handlerContext.grid, name, handlerContext.keyword.location());
            if (well2.updateConnections(connections, handlerContext.grid)) {
                handlerContext.state().wells.update( well2 );
                handlerContext.record_well_structure_change();
            }
            handlerContext.state().wellgroup_events().addEvent( name, ScheduleEvents::COMPLETION_CHANGE);
            const auto& md = connections->getMD();
            if (!std::is_sorted(std::begin(md), std::end(md))) {
                auto msg = fmt::format("Well {} measured depth column is not strictly increasing", name);
                throw OpmInputError(msg, handlerContext.keyword.location());
            }
        }
    }
    handlerContext.state().events().addEvent(ScheduleEvents::COMPLETION_CHANGE);
}

void handleWHISTCTL(HandlerContext& handlerContext)
{
    const auto& record = handlerContext.keyword.getRecord(0);
    const std::string& cmodeString = record.getItem("CMODE").getTrimmedString(0);
    const auto controlMode = WellProducerCModeFromString(cmodeString);

    if (controlMode != Well::ProducerCMode::NONE && !Well::WellProductionProperties::effectiveHistoryProductionControl(controlMode) ) {
        std::string msg = "The WHISTCTL keyword specifies an un-supported control mode " + cmodeString
            + ", which makes WHISTCTL keyword not affect the simulation at all";
        OpmLog::warning(msg);
    }
    handlerContext.state().update_whistctl( controlMode );

    const std::string bhp_terminate = record.getItem("BPH_TERMINATE").getTrimmedString(0);
    if (bhp_terminate == "YES") {
        std::string msg_fmt = "Problem with {keyword}\n"
                              "In {file} line {line}\n"
                              "Setting item 2 in {keyword} to 'YES' to stop the run is not supported";
        handlerContext.parseContext.handleError( ParseContext::UNSUPPORTED_TERMINATE_IF_BHP , msg_fmt, handlerContext.keyword.location(), handlerContext.errors );
    }

    for (const auto& well_ref : handlerContext.state().wells()) {
        auto well2 = well_ref.get();
        auto prop = std::make_shared<Well::WellProductionProperties>(well2.getProductionProperties());

        if (prop->whistctl_cmode != controlMode) {
            prop->whistctl_cmode = controlMode;
            well2.updateProduction(prop);
            handlerContext.state().wells.update( std::move(well2) );
        }
    }
}

void handleWLIST(HandlerContext& handlerContext)
{
    const std::string legal_actions = "NEW:ADD:DEL:MOV";
    for (const auto& record : handlerContext.keyword) {
        const std::string& name = record.getItem("NAME").getTrimmedString(0);
        const std::string& action = record.getItem("ACTION").getTrimmedString(0);
        const std::vector<std::string>& well_args = record.getItem("WELLS").getData<std::string>();
        std::vector<std::string> wells;
        auto new_wlm = handlerContext.state().wlist_manager.get();

        if (legal_actions.find(action) == std::string::npos)
            throw std::invalid_argument("The action:" + action + " is not recognized.");

        for (const auto& well_arg : well_args) {
            // does not use overload for context to avoid throw
            const auto names = handlerContext.wellNames(well_arg, true);
            if (names.empty() && well_arg.find("*") == std::string::npos) {
                std::string msg_fmt = "Problem with {keyword}\n"
                                      "In {file} line {line}\n"
                                      "The well '" + well_arg + "' has not been defined with WELSPECS and will not be added to the list.";
                const auto& parseContext = handlerContext.parseContext;
                parseContext.handleError(ParseContext::SCHEDULE_INVALID_NAME, msg_fmt, handlerContext.keyword.location(),handlerContext.errors);
                continue;
            }

            std::move(names.begin(), names.end(), std::back_inserter(wells));
        }

        if (name[0] != '*')
            throw std::invalid_argument("The list name in WLIST must start with a '*'");

        if (action == "NEW") {
            new_wlm.newList(name, wells);
        }

        if (!new_wlm.hasList(name))
            throw std::invalid_argument("Invalid well list: " + name);

        if (action == "MOV") {
            for (const auto& well : wells) {
                new_wlm.delWell(well);
            }
        }

        if (action == "DEL") {
            for (const auto& well : wells) {
                new_wlm.delWListWell(well, name);
            }
        } else if (action != "NEW"){
            for (const auto& well : wells) {
                new_wlm.addWListWell(well, name);
            }
        }
        handlerContext.state().wlist_manager.update( std::move(new_wlm) );
    }
}

void handleWPAVE(HandlerContext& handlerContext)
{
    auto wpave = PAvg( handlerContext.keyword.getRecord(0) );

    if (wpave.inner_weight() > 1.0) {
        const auto reason =
            fmt::format("Inner block weighting F1 "
                        "must not exceed 1.0. Got {}",
                        wpave.inner_weight());

        throw OpmInputError {
            reason, handlerContext.keyword.location()
        };
    }

    if ((wpave.conn_weight() < 0.0) ||
        (wpave.conn_weight() > 1.0))
    {
        const auto reason =
            fmt::format("Connection weighting factor F2 "
                        "must be between zero and one "
                        "inclusive. Got {} instead.",
                        wpave.conn_weight());

        throw OpmInputError {
            reason, handlerContext.keyword.location()
        };
    }

    for (const auto& wname : handlerContext.state().well_order()) {
        const auto& well = handlerContext.state().wells.get(wname);
        if (well.pavg() != wpave) {
            auto new_well = handlerContext.state().wells.get(wname);
            new_well.updateWPAVE(wpave);
            handlerContext.state().wells.update(std::move(new_well));
        }
    }

    handlerContext.state().pavg.update(std::move(wpave));
}

void handleWPAVEDEP(HandlerContext& handlerContext)
{
    for (const auto& record : handlerContext.keyword) {
        const std::string& wellNamePattern = record.getItem<ParserKeywords::WPAVEDEP::WELL>().getTrimmedString(0);
        const auto well_names = handlerContext.wellNames(wellNamePattern, false);

        if (well_names.empty()) {
            handlerContext.invalidNamePattern(wellNamePattern);
        }

        const auto& item = record.getItem<ParserKeywords::WPAVEDEP::REFDEPTH>();
        if (item.hasValue(0)) {
            auto ref_depth = item.getSIDouble(0);
            for (const auto& well_name : well_names) {
                auto well = handlerContext.state().wells.get(well_name);
                well.updateWPaveRefDepth( ref_depth );
                handlerContext.state().wells.update( std::move(well) );
            }
        }
    }
}

void handleWTEST(HandlerContext& handlerContext)
{
    auto new_config = handlerContext.state().wtest_config.get();
    for (const auto& record : handlerContext.keyword) {
        const std::string& wellNamePattern = record.getItem("WELL").getTrimmedString(0);
        const auto well_names = handlerContext.wellNames(wellNamePattern);
        if (well_names.empty()) {
            handlerContext.invalidNamePattern(wellNamePattern);
        }

        const double test_interval = record.getItem("INTERVAL").getSIDouble(0);
        const std::string& reasons = record.getItem("REASON").get<std::string>(0);
        const int num_test = record.getItem("TEST_NUM").get<int>(0);
        const double startup_time = record.getItem("START_TIME").getSIDouble(0);

        for (const auto& well_name : well_names) {
            if (reasons.empty())
                new_config.drop_well(well_name);
            else
                new_config.add_well(well_name, reasons, test_interval, num_test, startup_time, handlerContext.currentStep);
        }
    }
    handlerContext.state().wtest_config.update( std::move(new_config) );
}

}

std::vector<std::pair<std::string,KeywordHandlers::handler_function>>
getWellHandlers()
{
    return {
        { "WCONHIST", &handleWCONHIST },
        { "WCONINJE", &handleWCONINJE },
        { "WCONINJH", &handleWCONINJH },
        { "WCONPROD", &handleWCONPROD },

        { "WCYCLE",   &handleWCYCLE   },

        { "WELOPEN" , &handleWELOPEN  },
        { "WELLSTRE", &handleWELLSTRE },
        { "WELSPECS", &handleWELSPECS },
        { "WELSPECL", &handleWELSPECL },
        { "WELTARG" , &handleWELTARG  },
        { "WELTRAJ" , &handleWELTRAJ  },
        { "WHISTCTL", &handleWHISTCTL },
        { "WINJGAS",  &handleWINJGAS  },
        { "WLIST"   , &handleWLIST    },
        { "WPAVE"   , &handleWPAVE    },
        { "WPAVEDEP", &handleWPAVEDEP },
        { "WTEST"   , &handleWTEST    },
    };
}

}
