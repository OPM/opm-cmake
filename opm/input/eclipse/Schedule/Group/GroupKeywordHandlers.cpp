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

#include "GroupKeywordHandlers.hpp"

#include <opm/common/utility/String.hpp>

#include <opm/input/eclipse/Deck/DeckKeyword.hpp>

#include <opm/input/eclipse/Parser/ParseContext.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/G.hpp>

#include <opm/input/eclipse/Schedule/GasLiftOpt.hpp>
#include <opm/input/eclipse/Schedule/Group/GConSale.hpp>
#include <opm/input/eclipse/Schedule/Group/GConSump.hpp>
#include <opm/input/eclipse/Schedule/Group/GSatProd.hpp>
#include <opm/input/eclipse/Schedule/Group/Group.hpp>
#include <opm/input/eclipse/Schedule/Group/GroupEconProductionLimits.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRateConfig.hpp>
#include <opm/input/eclipse/Schedule/ScheduleState.hpp>
#include <opm/input/eclipse/Schedule/ScheduleStatic.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQActive.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQConfig.hpp>

#include "../HandlerContext.hpp"

#include <fmt/format.h>

namespace Opm {

namespace {

/*
  The function trim_wgname() is used to trim the leading and trailing spaces
  away from the group and well arguments given in the GRUPTREE
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

void handleGCONINJE(HandlerContext& handlerContext)
{
    using GI = ParserKeywords::GCONINJE;
    const auto& keyword = handlerContext.keyword;
    for (const auto& record : keyword) {
        const std::string& groupNamePattern = record.getItem<GI::GROUP>().getTrimmedString(0);
        const auto group_names = handlerContext.groupNames(groupNamePattern);
        if (group_names.empty()) {
            handlerContext.invalidNamePattern(groupNamePattern);
        }

        const Group::InjectionCMode controlMode = Group::InjectionCModeFromString(record.getItem<GI::CONTROL_MODE>().getTrimmedString(0));
        const Phase phase = get_phase( record.getItem<GI::PHASE>().getTrimmedString(0));
        const auto surfaceInjectionRate = record.getItem<GI::SURFACE_TARGET>().get<UDAValue>(0);
        const auto reservoirInjectionRate = record.getItem<GI::RESV_TARGET>().get<UDAValue>(0);
        const auto reinj_target = record.getItem<GI::REINJ_TARGET>().get<UDAValue>(0);
        const auto voidage_target = record.getItem<GI::VOIDAGE_TARGET>().get<UDAValue>(0);
        const bool is_free = DeckItem::to_bool(record.getItem<GI::RESPOND_TO_PARENT>().getTrimmedString(0));

        std::optional<std::string> guide_rate_str;
        {
            const auto& item = record.getItem("GUIDE_RATE_DEF");
            if (item.hasValue(0)) {
                const auto& string_value = record.getItem("GUIDE_RATE_DEF").getTrimmedString(0);
                if (string_value.size() > 0)
                    guide_rate_str = string_value;
            }

        }

        for (const auto& group_name : group_names) {
            const bool is_field { group_name == "FIELD" } ;

            auto guide_rate_def = Group::GuideRateInjTarget::NO_GUIDE_RATE;
            double guide_rate = 0;
            if (!is_field) {
                if (guide_rate_str) {
                    guide_rate_def = Group::GuideRateInjTargetFromString(guide_rate_str.value());
                    guide_rate = record.getItem("GUIDE_RATE").get<double>(0);
                }
            }

            {
                // FLD overrides item 8 (is_free i.e if FLD the group is available for higher up groups)
                const bool availableForGroupControl = (is_free || controlMode == Group::InjectionCMode::FLD)&& !is_field;
                auto new_group = handlerContext.state().groups.get(group_name);
                Group::GroupInjectionProperties injection{group_name};
                injection.phase = phase;
                injection.cmode = controlMode;
                injection.surface_max_rate = surfaceInjectionRate;
                injection.resv_max_rate = reservoirInjectionRate;
                injection.target_reinj_fraction = reinj_target;
                injection.target_void_fraction = voidage_target;
                injection.injection_controls = 0;
                injection.guide_rate = guide_rate;
                injection.guide_rate_def = guide_rate_def;
                injection.available_group_control = availableForGroupControl;

                if (!record.getItem("SURFACE_TARGET").defaultApplied(0))
                    injection.injection_controls += static_cast<int>(Group::InjectionCMode::RATE);

                if (!record.getItem("RESV_TARGET").defaultApplied(0))
                    injection.injection_controls += static_cast<int>(Group::InjectionCMode::RESV);

                if (!record.getItem("REINJ_TARGET").defaultApplied(0))
                    injection.injection_controls += static_cast<int>(Group::InjectionCMode::REIN);

                if (!record.getItem("VOIDAGE_TARGET").defaultApplied(0))
                    injection.injection_controls += static_cast<int>(Group::InjectionCMode::VREP);

                if (record.getItem("REINJECT_GROUP").hasValue(0))
                    injection.reinj_group = record.getItem("REINJECT_GROUP").getTrimmedString(0);

                if (record.getItem("VOIDAGE_GROUP").hasValue(0))
                    injection.voidage_group = record.getItem("VOIDAGE_GROUP").getTrimmedString(0);

                if (new_group.updateInjection(injection)) {
                    auto new_config = handlerContext.state().guide_rate();
                    new_config.update_injection_group(group_name, injection);
                    handlerContext.state().guide_rate.update( std::move(new_config));

                    handlerContext.state().groups.update( std::move(new_group));
                    handlerContext.state().events().addEvent( ScheduleEvents::GROUP_INJECTION_UPDATE );
                    handlerContext.state().wellgroup_events().addEvent( group_name, ScheduleEvents::GROUP_INJECTION_UPDATE);

                    auto udq_active = handlerContext.state().udq_active.get();
                    if (injection.updateUDQActive(handlerContext.state().udq.get(), udq_active))
                        handlerContext.state().udq_active.update( std::move(udq_active));
                }
            }
        }
    }
}

void handleGCONPROD(HandlerContext& handlerContext)
{
    const auto& keyword = handlerContext.keyword;
    for (const auto& record : keyword) {
        const std::string& groupNamePattern = record.getItem("GROUP").getTrimmedString(0);
        const auto group_names = handlerContext.groupNames(groupNamePattern);
        if (group_names.empty()) {
            handlerContext.invalidNamePattern(groupNamePattern);
        }

        const Group::ProductionCMode controlMode = Group::ProductionCModeFromString(record.getItem("CONTROL_MODE").getTrimmedString(0));
        Group::GroupLimitAction groupLimitAction;
        groupLimitAction.allRates = Group::ExceedActionFromString(record.getItem("EXCEED_PROC").getTrimmedString(0));

        groupLimitAction.water = record.getItem("WATER_EXCEED_PROCEDURE").defaultApplied(0)
        ? groupLimitAction.allRates
        : Group::ExceedActionFromString(record.getItem("WATER_EXCEED_PROCEDURE").getTrimmedString(0));
        
        groupLimitAction.gas = record.getItem("GAS_EXCEED_PROCEDURE").defaultApplied(0)
        ? groupLimitAction.allRates
        : Group::ExceedActionFromString(record.getItem("GAS_EXCEED_PROCEDURE").getTrimmedString(0));

        groupLimitAction.liquid = record.getItem("LIQUID_EXCEED_PROCEDURE").defaultApplied(0)
        ? groupLimitAction.allRates
        : Group::ExceedActionFromString(record.getItem("LIQUID_EXCEED_PROCEDURE").getTrimmedString(0));

        const bool respond_to_parent = DeckItem::to_bool(record.getItem("RESPOND_TO_PARENT").getTrimmedString(0));

        const auto oil_target = record.getItem("OIL_TARGET").get<UDAValue>(0);
        const auto gas_target = record.getItem("GAS_TARGET").get<UDAValue>(0);
        const auto water_target = record.getItem("WATER_TARGET").get<UDAValue>(0);
        const auto liquid_target = record.getItem("LIQUID_TARGET").get<UDAValue>(0);
        const auto resv_target = record.getItem("RESERVOIR_FLUID_TARGET").get<UDAValue>(0);

        const bool apply_default_oil_target = record.getItem("OIL_TARGET").defaultApplied(0);
        const bool apply_default_gas_target = record.getItem("GAS_TARGET").defaultApplied(0);
        const bool apply_default_water_target = record.getItem("WATER_TARGET").defaultApplied(0);
        const bool apply_default_liquid_target = record.getItem("LIQUID_TARGET").defaultApplied(0);
        const bool apply_default_resv_target = record.getItem("RESERVOIR_FLUID_TARGET").defaultApplied(0);

        std::optional<std::string> guide_rate_str;
        {
            const auto& item = record.getItem("GUIDE_RATE_DEF");
            if (item.hasValue(0)) {
                const auto& string_value = record.getItem("GUIDE_RATE_DEF").getTrimmedString(0);
                if (string_value.size() > 0)
                    guide_rate_str = string_value;
            }

        }

        for (const auto& group_name : group_names) {
            const bool is_field { group_name == "FIELD" } ;

            auto guide_rate_def = Group::GuideRateProdTarget::NO_GUIDE_RATE;
            double guide_rate = 0;
            if (!is_field) {
                if (guide_rate_str) {
                    guide_rate_def = Group::GuideRateProdTargetFromString(guide_rate_str.value());

                    if ((guide_rate_def == Group::GuideRateProdTarget::INJV ||
                         guide_rate_def == Group::GuideRateProdTarget::POTN ||
                         guide_rate_def == Group::GuideRateProdTarget::FORM)) {
                        std::string msg_fmt = "Problem with {keyword}\n"
                            "In {file} line {line}\n"
                            "The supplied guide rate will be ignored";
                        const auto& parseContext = handlerContext.parseContext;
                        auto& errors = handlerContext.errors;
                        parseContext.handleError(ParseContext::SCHEDULE_IGNORED_GUIDE_RATE, msg_fmt, keyword.location(), errors);
                    } else {
                        guide_rate = record.getItem("GUIDE_RATE").get<double>(0);
                        if (guide_rate == 0)
                            guide_rate_def = Group::GuideRateProdTarget::POTN;
                    }
                }
            }

            {
                // FLD overrides item 8 (respond_to_parent i.e if FLD the group is available for higher up groups)
                const bool availableForGroupControl { (respond_to_parent || controlMode == Group::ProductionCMode::FLD) && !is_field } ;
                auto new_group = handlerContext.state().groups.get(group_name);
                Group::GroupProductionProperties production(handlerContext.static_schedule().m_unit_system, group_name);
                production.cmode = controlMode;
                production.oil_target = oil_target;
                production.gas_target = gas_target;
                production.water_target = water_target;
                production.liquid_target = liquid_target;
                production.guide_rate = guide_rate;
                production.guide_rate_def = guide_rate_def;
                production.resv_target = resv_target;
                production.available_group_control = availableForGroupControl;
                production.group_limit_action = groupLimitAction;

                production.production_controls = 0;
                // GCONPROD
                // 'G1' 'ORAT' 1000 100 200 300 NONE =>  constraints 100,200,300 should be ignored
                //
                // GCONPROD
                // 'G1' 'ORAT' 1000 100 200 300 RATE =>  constraints 100,200,300 should be honored
                if (production.cmode == Group::ProductionCMode::ORAT ||
                    (groupLimitAction.allRates != Group::ExceedAction::NONE &&
                    !apply_default_oil_target)) {
                    production.production_controls |= static_cast<int>(Group::ProductionCMode::ORAT);
                }
                if (production.cmode == Group::ProductionCMode::WRAT ||
                    ((groupLimitAction.water != Group::ExceedAction::NONE) &&
                    !apply_default_water_target)) {
                    production.production_controls |= static_cast<int>(Group::ProductionCMode::WRAT);
                }
                if (production.cmode == Group::ProductionCMode::GRAT ||
                    ((groupLimitAction.gas  != Group::ExceedAction::NONE) &&
                    !apply_default_gas_target)) {
                    production.production_controls |= static_cast<int>(Group::ProductionCMode::GRAT);
                }
                if (production.cmode == Group::ProductionCMode::LRAT ||
                    ((groupLimitAction.liquid != Group::ExceedAction::NONE) &&
                    !apply_default_liquid_target)) {
                    production.production_controls |= static_cast<int>(Group::ProductionCMode::LRAT);
                }

                if (!apply_default_resv_target)
                    production.production_controls |= static_cast<int>(Group::ProductionCMode::RESV);

                if (new_group.updateProduction(production)) {
                    auto new_config = handlerContext.state().guide_rate();
                    new_config.update_production_group(new_group);
                    handlerContext.state().guide_rate.update( std::move(new_config));

                    handlerContext.state().groups.update( std::move(new_group));
                    handlerContext.state().events().addEvent(ScheduleEvents::GROUP_PRODUCTION_UPDATE);
                    handlerContext.state().wellgroup_events().addEvent( group_name, ScheduleEvents::GROUP_PRODUCTION_UPDATE);

                    auto udq_active = handlerContext.state().udq_active.get();
                    if (production.updateUDQActive(handlerContext.state().udq.get(), udq_active))
                        handlerContext.state().udq_active.update( std::move(udq_active));
                }
            }
        }
    }
}

void handleGSATPROD(HandlerContext& handlerContext)
{
    const auto& keyword = handlerContext.keyword;
    auto new_gsatprod = handlerContext.state().gsatprod.get();

    for (const auto& record : keyword) {
        const std::string& groupNamePattern = record.getItem("SATELLITE_GROUP_NAME_OR_GROUP_NAME_ROOT").getTrimmedString(0);
        const auto group_names = handlerContext.groupNames(groupNamePattern);
        if (group_names.empty()) {
            handlerContext.invalidNamePattern(groupNamePattern);
        }
        const auto oil_rate = record.getItem("OIL_PRODUCTION_RATE").getSIDouble(0);
        const auto gas_rate = record.getItem("GAS_PRODUCTION_RATE").getSIDouble(0);
        const auto water_rate = record.getItem("WATER_PRODUCTION_RATE").getSIDouble(0);
        const auto resv_rate = record.getItem("RES_FLUID_VOL_PRODUCTION_RATE").getSIDouble(0);
        const auto glift_rate = record.getItem("LIFT_GAS_SUPPLY_RATE").getSIDouble(0);

        for (const auto& group_name : group_names) {
            const bool is_field { group_name == "FIELD" } ;
            if (is_field) {
                std::string msg_fmt = "Problem with {keyword}\n"
                    "In {file} line {line}\n"
                    "GSATPROD group cannot be named FIELD ";
                const auto& parseContext = handlerContext.parseContext;
                auto& errors = handlerContext.errors;
                parseContext.handleError(ParseContext::SCHEDULE_GROUP_ERROR, msg_fmt, keyword.location(), errors);
            }
            new_gsatprod.assign(group_name, oil_rate, gas_rate, water_rate, resv_rate, glift_rate);
        }
    }
    handlerContext.state().gsatprod.update( std::move(new_gsatprod) );
}

void handleGCONSALE(HandlerContext& handlerContext)
{
    auto new_gconsale = handlerContext.state().gconsale.get();
    for (const auto& record : handlerContext.keyword) {
        const std::string& groupName = record.getItem("GROUP").getTrimmedString(0);
        auto sales_target = record.getItem("SALES_TARGET").get<UDAValue>(0);
        auto max_rate = record.getItem("MAX_SALES_RATE").get<UDAValue>(0);
        auto min_rate = record.getItem("MIN_SALES_RATE").get<UDAValue>(0);
        std::string procedure = record.getItem("MAX_PROC").getTrimmedString(0);
        auto udqconfig = handlerContext.state().udq.get().params().undefinedValue();

        new_gconsale.add(groupName, sales_target, max_rate, min_rate, procedure,
                         udqconfig, handlerContext.static_schedule().m_unit_system);

        auto new_group = handlerContext.state().groups.get( groupName );
        Group::GroupInjectionProperties injection{groupName};
        injection.phase = Phase::GAS;
        if (new_group.updateInjection(injection))
            handlerContext.state().groups.update(new_group);
    }
    handlerContext.state().gconsale.update( std::move(new_gconsale) );
}

void handleGCONSUMP(HandlerContext& handlerContext)
{
    auto new_gconsump = handlerContext.state().gconsump.get();
    for (const auto& record : handlerContext.keyword) {
        const std::string& groupName = record.getItem("GROUP").getTrimmedString(0);
        auto consumption_rate = record.getItem("GAS_CONSUMP_RATE").get<UDAValue>(0);
        auto import_rate = record.getItem("GAS_IMPORT_RATE").get<UDAValue>(0);

        std::string network_node_name;
        auto network_node = record.getItem("NETWORK_NODE");
        if (!network_node.defaultApplied(0))
            network_node_name = network_node.getTrimmedString(0);

        auto udqconfig = handlerContext.state().udq.get().params().undefinedValue();

        new_gconsump.add(groupName, consumption_rate, import_rate, network_node_name,
                         udqconfig, handlerContext.static_schedule().m_unit_system);
    }
    handlerContext.state().gconsump.update( std::move(new_gconsump) );
}

void handleGECON(HandlerContext& handlerContext)
{
    auto gecon = handlerContext.state().gecon();
    const auto& keyword = handlerContext.keyword;
    auto report_step = handlerContext.currentStep;
    for (const auto& record : keyword) {
        const std::string& groupNamePattern
            = record.getItem<ParserKeywords::GECON::GROUP>().getTrimmedString(0);
        const auto group_names = handlerContext.groupNames(groupNamePattern);
        if (group_names.empty()) {
            handlerContext.invalidNamePattern(groupNamePattern);
        }
        for (const auto& gname : group_names) {
            gecon.add_group(report_step, gname, record);
        }
    }
    handlerContext.state().gecon.update(std::move(gecon));
}

void handleGEFAC(HandlerContext& handlerContext)
{
    for (const auto& record : handlerContext.keyword) {
        const std::string& groupNamePattern = record.getItem("GROUP").getTrimmedString(0);
        const auto group_names = handlerContext.groupNames(groupNamePattern);
        if (group_names.empty()) {
            handlerContext.invalidNamePattern(groupNamePattern);
        }

        const bool use_efficiency_in_network = DeckItem::to_bool(record.getItem("USE_GEFAC_IN_NETWORK").getTrimmedString(0));
        const auto gefac = record.getItem("EFFICIENCY_FACTOR").get<double>(0);

        for (const auto& group_name : group_names) {
            auto new_group = handlerContext.state().groups.get(group_name);
            if (new_group.update_gefac(gefac, use_efficiency_in_network)) {
                handlerContext.state().wellgroup_events().addEvent( group_name, ScheduleEvents::WELLGROUP_EFFICIENCY_UPDATE);
                handlerContext.state().events().addEvent( ScheduleEvents::WELLGROUP_EFFICIENCY_UPDATE );
                handlerContext.state().groups.update(std::move(new_group));
            }
        }
    }
}

void handleGPMAINT(HandlerContext& handlerContext)
{
    for (const auto& record : handlerContext.keyword) {
        const std::string& groupNamePattern = record.getItem("GROUP").getTrimmedString(0);
        const auto group_names = handlerContext.groupNames(groupNamePattern);
        if (group_names.empty()) {
            handlerContext.invalidNamePattern(groupNamePattern);
        }

        const auto& target_string = record.getItem<ParserKeywords::GPMAINT::FLOW_TARGET>().get<std::string>(0);

        for (const auto& group_name : group_names) {
            auto new_group = handlerContext.state().groups.get(group_name);
            if (target_string == "NONE") {
                new_group.set_gpmaint();
            } else {
                GPMaint gpmaint(handlerContext.currentStep, record);
                new_group.set_gpmaint(std::move(gpmaint));
            }
            handlerContext.state().groups.update( std::move(new_group) );
        }
    }
}

void handleGRUPTREE(HandlerContext& handlerContext)
{
    for (const auto& record : handlerContext.keyword) {
        const std::string& childName = trim_wgname(handlerContext.keyword,
                                                   record.getItem("CHILD_GROUP").get<std::string>(0),
                                                   handlerContext.parseContext,
                                                   handlerContext.errors);
        const std::string& parentName = trim_wgname(handlerContext.keyword,
                                                    record.getItem("PARENT_GROUP").get<std::string>(0),
                                                    handlerContext.parseContext,
                                                    handlerContext.errors);

        if (!handlerContext.state().groups.has(childName)) {
            handlerContext.addGroup(childName);
        }

        if (!handlerContext.state().groups.has(parentName)) {
            handlerContext.addGroup(parentName);
        }

        handlerContext.addGroupToGroup(parentName, childName);
    }
}

}

std::vector<std::pair<std::string,KeywordHandlers::handler_function>>
getGroupHandlers()
{
    return {
        { "GCONINJE", &handleGCONINJE },
        { "GCONPROD", &handleGCONPROD },
        { "GCONSALE", &handleGCONSALE },
        { "GCONSUMP", &handleGCONSUMP },
        { "GSATPROD", &handleGSATPROD },
        { "GECON"   , &handleGECON    },
        { "GEFAC"   , &handleGEFAC    },
        { "GPMAINT" , &handleGPMAINT  },
        { "GRUPTREE", &handleGRUPTREE },
    };
}

}
