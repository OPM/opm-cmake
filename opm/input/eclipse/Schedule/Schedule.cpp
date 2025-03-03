/*
  Copyright 2013 Statoil ASA.

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

#include <opm/input/eclipse/Schedule/Schedule.hpp>

#include <opm/io/eclipse/rst/state.hpp>

#include <opm/common/OpmLog/LogUtil.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/utility/OpmInputError.hpp>
#include <opm/input/eclipse/Parser/InputErrorAction.hpp>
#include <opm/common/utility/String.hpp>
#include <opm/common/utility/numeric/cmp.hpp>
#include <opm/common/utility/shmatch.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/TracerConfig.hpp>

#include <opm/input/eclipse/Python/Python.hpp>

#include <opm/input/eclipse/Schedule/Action/ActionResult.hpp>
#include <opm/input/eclipse/Schedule/Action/Actions.hpp>
#include <opm/input/eclipse/Schedule/Action/ActionX.hpp>
#include <opm/input/eclipse/Schedule/Action/State.hpp>
#include <opm/input/eclipse/Schedule/Action/SimulatorUpdate.hpp>
#include <opm/input/eclipse/Schedule/Group/GConSale.hpp>
#include <opm/input/eclipse/Schedule/Group/GConSump.hpp>
#include <opm/input/eclipse/Schedule/Group/GSatProd.hpp>
#include <opm/input/eclipse/Schedule/Group/GroupEconProductionLimits.hpp>
#include <opm/input/eclipse/Schedule/Group/GTNode.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRateConfig.hpp>
#include <opm/input/eclipse/Schedule/GasLiftOpt.hpp>
#include <opm/input/eclipse/Schedule/MSW/SegmentMatcher.hpp>
#include <opm/input/eclipse/Schedule/MSW/SICD.hpp>
#include <opm/input/eclipse/Schedule/MSW/Valve.hpp>
#include <opm/input/eclipse/Schedule/MSW/WellSegments.hpp>
#include <opm/input/eclipse/Schedule/Network/Balance.hpp>
#include <opm/input/eclipse/Schedule/Network/ExtNetwork.hpp>
#include <opm/input/eclipse/Schedule/Network/Node.hpp>
#include <opm/input/eclipse/Schedule/OilVaporizationProperties.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/ReservoirCouplingInfo.hpp>
#include <opm/input/eclipse/Schedule/RFTConfig.hpp>
#include <opm/input/eclipse/Schedule/RPTConfig.hpp>
#include <opm/input/eclipse/Schedule/ScheduleGrid.hpp>
#include <opm/input/eclipse/Schedule/SummaryState.hpp>
#include <opm/input/eclipse/Schedule/Tuning.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQActive.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQConfig.hpp>
#include <opm/input/eclipse/Schedule/Well/WList.hpp>
#include <opm/input/eclipse/Schedule/Well/WListManager.hpp>
#include <opm/input/eclipse/Schedule/Well/WellBrineProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>
#include <opm/input/eclipse/Schedule/Well/WellEnums.hpp>
#include <opm/input/eclipse/Schedule/Well/WellFoamProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/WellMatcher.hpp>
#include <opm/input/eclipse/Schedule/Well/WellMICPProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/WellPolymerProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTestConfig.hpp>

#include <opm/input/eclipse/Units/Dimension.hpp>
#include <opm/input/eclipse/Units/UnitSystem.hpp>
#include <opm/input/eclipse/Units/Units.hpp>

#include <opm/input/eclipse/Parser/ErrorGuard.hpp>
#include <opm/input/eclipse/Parser/ParseContext.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/A.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/B.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/C.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/E.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/W.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Deck/DeckItem.hpp>
#include <opm/input/eclipse/Deck/DeckKeyword.hpp>
#include <opm/input/eclipse/Deck/DeckRecord.hpp>
#include <opm/input/eclipse/Deck/DeckSection.hpp>

#include "HandlerContext.hpp"
#include "KeywordHandlers.hpp"
#include "MSW/Compsegs.hpp"
#include "MSW/WelSegsSet.hpp"
#include "Well/injection.hpp"

#include <algorithm>
#include <ctime>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include <fmt/chrono.h>
#include <fmt/format.h>

namespace {

    bool name_match_any(const std::unordered_set<std::string>& patterns,
                        const std::string& name)
    {
        return std::any_of(patterns.begin(), patterns.end(),
                           [&name](const auto& pattern)
                           { return Opm::shmatch(pattern, name); });
    }
}

namespace Opm {

    Schedule::Schedule( const Deck& deck,
                        const EclipseGrid& ecl_grid,
                        const FieldPropsManager& fp,
                        const Runspec &runspec,
                        const ParseContext& parseContext,
                        ErrorGuard& errors,
                        std::shared_ptr<const Python> python,
                        const bool lowActionParsingStrictness,
                        const bool slave_mode,
                        bool keepKeywords,
                        const std::optional<int>& output_interval,
                        const RestartIO::RstState * rst,
                        const TracerConfig * tracer_config)
    try :
        m_static(python, ScheduleRestartInfo(rst, deck), deck, runspec,
                 output_interval, parseContext, errors, slave_mode)
        , m_sched_deck(TimeService::from_time_t(runspec.start_time()), deck, m_static.rst_info)
        , completed_cells(ecl_grid.getNX(), ecl_grid.getNY(), ecl_grid.getNZ())
        , m_lowActionParsingStrictness(lowActionParsingStrictness)
    {
        this->restart_output.resize(this->m_sched_deck.size());
        this->restart_output.clearRemainingEvents(0);
        this->simUpdateFromPython = std::make_shared<SimulatorUpdate>();

        //const ScheduleGridWrapper gridWrapper { grid } ;
        ScheduleGrid grid(ecl_grid, fp, this->completed_cells);

        if (!keepKeywords) {
            const auto& section = SCHEDULESection(deck);
            keepKeywords = section.has_keyword("ACTIONX") ||
                           section.has_keyword("PYACTION");
        }

        if (rst) {
            if (!tracer_config) {
                throw std::logic_error("Bug: when loading from restart a valid TracerConfig object must be supplied");
            }

            if (!keepKeywords) {
                keepKeywords = !rst->actions.empty();
            }

            auto restart_step = this->m_static.rst_info.report_step;
            this->iterateScheduleSection(0, restart_step, parseContext, errors,
                                         grid, nullptr, "", keepKeywords);
            this->load_rst(*rst, *tracer_config, grid, fp);
            if (! this->restart_output.writeRestartFile(restart_step))
                this->restart_output.addRestartOutput(restart_step);
            this->iterateScheduleSection(restart_step, this->m_sched_deck.size(),
                                         parseContext, errors, grid, nullptr, "", keepKeywords);
            // Events added during restart reading well be added to previous step, but need to be active at the
            // restart step to ensure well potentials and guide rates are available at the first step.
            const auto prev_step = std::max(static_cast<int>(restart_step-1), 0);
            this->snapshots[restart_step].wellgroup_events().merge(this->snapshots[prev_step].wellgroup_events());
            this->snapshots[restart_step].events().merge(this->snapshots[prev_step].events());
        } else {
            this->iterateScheduleSection(0, this->m_sched_deck.size(),
                                         parseContext, errors, grid, nullptr, "", keepKeywords);
        }
    }
    catch (const OpmInputError& opm_error) {
        OpmLog::error(opm_error.what());
        throw;
    }
    catch (const std::exception& std_error) {
        OpmLog::error(fmt::format("An error occurred while creating the reservoir schedule\n"
                                  "Internal error: {}", std_error.what()));
        throw;
    }

    template <typename T>
    Schedule::Schedule( const Deck& deck,
                        const EclipseGrid& grid,
                        const FieldPropsManager& fp,
                        const Runspec &runspec,
                        const ParseContext& parseContext,
                        T&& errors,
                        std::shared_ptr<const Python> python,
                        const bool lowActionParsingStrictness,
                        const bool slave_mode,
                        const bool keepKeywords,
                        const std::optional<int>& output_interval,
                        const RestartIO::RstState * rst,
                        const TracerConfig* tracer_config)
        : Schedule(deck,
                   grid,
                   fp,
                   runspec,
                   parseContext,
                   errors,
                   python,
                   lowActionParsingStrictness,
                   slave_mode,
                   keepKeywords,
                   output_interval,
                   rst,
                   tracer_config)
    {}

    Schedule::Schedule( const Deck& deck,
                        const EclipseGrid& grid,
                        const FieldPropsManager& fp,
                        const Runspec &runspec,
                        std::shared_ptr<const Python> python,
                        const bool lowActionParsingStrictness,
                        const bool slave_mode,
                        const bool keepKeywords,
                        const std::optional<int>& output_interval,
                        const RestartIO::RstState * rst,
                        const TracerConfig* tracer_config)
        : Schedule(deck,
                   grid,
                   fp,
                   runspec,
                   ParseContext(),
                   ErrorGuard(),
                   python,
                   lowActionParsingStrictness,
                   slave_mode,
                   keepKeywords,
                   output_interval,
                   rst,
                   tracer_config)
    {}

    Schedule::Schedule(const Deck& deck,
                       const EclipseState& es,
                       const ParseContext& parse_context,
                       ErrorGuard& errors,
                       std::shared_ptr<const Python> python,
                       bool lowActionParsingStrictness,
                       const bool slave_mode,
                       const bool keepKeywords,
                       const std::optional<int>& output_interval,
                       const RestartIO::RstState * rst)
        : Schedule(deck,
                   es.getInputGrid(),
                   es.fieldProps(),
                   es.runspec(),
                   parse_context,
                   errors,
                   python,
                   lowActionParsingStrictness,
                   slave_mode,
                   keepKeywords,
                   output_interval,
                   rst,
                   &es.tracer())
    {}

    template <typename T>
    Schedule::Schedule(const Deck& deck,
                       const EclipseState& es,
                       const ParseContext& parse_context,
                       T&& errors,
                       std::shared_ptr<const Python> python,
                       bool lowActionParsingStrictness,
                       const bool slave_mode,
                       const bool keepKeywords,
                       const std::optional<int>& output_interval,
                       const RestartIO::RstState * rst)
        : Schedule(deck,
                   es.getInputGrid(),
                   es.fieldProps(),
                   es.runspec(),
                   parse_context,
                   errors,
                   python,
                   lowActionParsingStrictness,
                   slave_mode,
                   keepKeywords,
                   output_interval,
                   rst,
                   &es.tracer())
    {}


    Schedule::Schedule(const Deck& deck,
                       const EclipseState& es,
                       std::shared_ptr<const Python> python,
                       bool lowActionParsingStrictness,
                       const bool slave_mode,
                       const bool keepKeywords,
                       const std::optional<int>& output_interval,
                       const RestartIO::RstState * rst)
        : Schedule(deck,
                   es,
                   ParseContext(),
                   ErrorGuard(),
                   python,
                   lowActionParsingStrictness,
                   slave_mode,
                   keepKeywords,
                   output_interval,
                   rst)
    {}

    Schedule::Schedule(const Deck& deck,
                       const EclipseState& es,
                       const std::optional<int>& output_interval,
                       const RestartIO::RstState * rst)
        : Schedule(deck,
                   es,
                   ParseContext(),
                   ErrorGuard(),
                   std::make_shared<const Python>(),
                   false,
                   /*slave_mode=*/false,
                   true,
                   output_interval,
                   rst)
        {}

    Schedule::Schedule(std::shared_ptr<const Python> python_handle) :
        m_static( python_handle )
    {
    }

    /*
      In general the serializationTestObject() instances are used as targets for
      deserialization, i.e. the serialized buffer is unpacked into this
      instance. However the Schedule object is a top level object, and the
      simulator will instantiate and manage a Schedule object to unpack into, so
      the instance created here is only for testing.
    */
    Schedule Schedule::serializationTestObject()
    {
        Schedule result;
        result.m_treat_critical_as_non_critical = false;
        result.m_static = ScheduleStatic::serializationTestObject();
        result.m_sched_deck = ScheduleDeck::serializationTestObject();
        result.action_wgnames = Action::WGNames::serializationTestObject();
        result.potential_wellopen_patterns = std::unordered_set<std::string> {"W1"};
        result.exit_status = EXIT_FAILURE;
        result.snapshots = { ScheduleState::serializationTestObject() };
        result.restart_output = WriteRestartFileEvents::serializationTestObject();
        result.completed_cells = CompletedCells::serializationTestObject();
        result.current_report_step = 0;
        result.m_lowActionParsingStrictness = false;
        result.simUpdateFromPython = std::make_shared<SimulatorUpdate>(SimulatorUpdate::serializationTestObject());

        return result;
    }

    std::time_t Schedule::getStartTime() const {
        return this->posixStartTime( );
    }

    std::time_t Schedule::posixStartTime() const {
        return std::chrono::system_clock::to_time_t(this->m_sched_deck[0].start_time());
    }

    std::time_t Schedule::posixEndTime() const {
        // This should indeed access the start_time() property of the last
        // snapshot.
        if (this->snapshots.size() > 0)
            return std::chrono::system_clock::to_time_t(this->snapshots.back().start_time());
        else
            return this->posixStartTime( );
    }


    void Schedule::handleKeyword(std::size_t currentStep,
                                 const ScheduleBlock& block,
                                 const DeckKeyword& keyword,
                                 const ParseContext& parseContext,
                                 ErrorGuard& errors,
                                 const ScheduleGrid& grid,
                                 const Action::Result::MatchingEntities& matches,
                                 const bool action_mode,
                                 SimulatorUpdate* sim_update,
                                 const std::unordered_map<std::string, double>* target_wellpi,
                                 std::unordered_map<std::string, double>& wpimult_global_factor,
                                 WelSegsSet* welsegs_wells,
                                 std::set<std::string>* compsegs_wells)
    {
        HandlerContext handlerContext { *this, block, keyword, grid, currentStep,
                                        matches, action_mode,
                                        parseContext, errors, sim_update, target_wellpi,
                                        wpimult_global_factor, welsegs_wells, compsegs_wells};

        if (!KeywordHandlers::getInstance().handleKeyword(handlerContext)) {
            OpmLog::warning(fmt::format("No handler registered for keyword {} "
                                        "in file {} line {}",
                            keyword.name(), keyword.location().filename,
                            keyword.location().lineno));
        }
    }

namespace {

class ScheduleLogger
{
public:
    enum class Stream
    {
        Info, Note, Debug,
    };

    ScheduleLogger(const Stream stream, const std::string& prefix, const KeywordLocation& location)
        : stream_      {stream}
        , prefix_      {prefix}
        , current_file_{location.filename}
    {
        switch (this->stream_) {
        case Stream::Info:
            this->log_function_ = &OpmLog::info;
            break;

        case Stream::Note:
            this->log_function_ = &OpmLog::note;
            break;

        case Stream::Debug:
            this->log_function_ = &OpmLog::debug;
            break;
        }
    }

    void operator()(const std::string& msg)
    {
        this->log_function_(this->format_message(msg));
    }

    void operator()(const std::vector<std::string>& msg_list)
    {
        std::for_each(msg_list.begin(), msg_list.end(),
                      [this](const std::string& record)
                      {
                          (*this)(record);
                      });
    }

    void info(const std::string& msg)
    {
        OpmLog::info(this->format_message(msg));
    }

    void complete_step(const std::string& msg)
    {
        (*this)(msg);

        this->step_count_ += 1;
        if (this->step_count_ == this->max_print_) {
            this->redirect(&OpmLog::note, {"Report limit reached, see PRT-file for remaining Schedule initialization.", ""});
        }
        else {
            // Blank line
            (*this)("");
        }
    };

    void restart()
    {
        this->step_count_ = 0;
        this->redirect(&OpmLog::info);
    }

    void location(const KeywordLocation& location)
    {
        if (this->current_file_ == location.filename) {
            return;
        }

        (*this)(fmt::format("Reading from: {} line {}", location.filename, location.lineno));
        this->current_file_ = location.filename;
    }

    static Stream select_stream(const bool log_to_debug,
                                const bool restart_skip);

private:
    using LogFunction = void (*) (const std::string&);

    Stream stream_{Stream::Info};
    std::size_t step_count_ = 0;
    std::size_t max_print_  = 5;
    std::string prefix_{};
    std::string current_file_{};

    LogFunction log_function_{};

    void redirect(const LogFunction               new_stream,
                  const std::vector<std::string>& messages = {})
    {
        if (this->stream_ == Stream::Debug) {
            // If we're writing to OpmLog::debug then continue to do so.
            return;
        }

        this->log_function_ = new_stream;
        (*this)(messages);
    }

    std::string format_message(const std::string& message) const
    {
        return fmt::format("{}{}", this->prefix_, message);
    }
};

    ScheduleLogger::Stream
    ScheduleLogger::select_stream(const bool log_to_debug,
                                  const bool restart_skip)
    {
        if (log_to_debug) {
            return ScheduleLogger::Stream::Debug;
        }

        if (restart_skip) {
            return ScheduleLogger::Stream::Note;
        }

        return ScheduleLogger::Stream::Info;
    }
}
} // end namespace Opm

namespace
{
/// \brief Check whether each MS well has COMPSEGS entry andissue error if not.
/// \param welsegs All wells with a WELSEGS entry together with the location.
/// \param compegs All wells with a COMPSEGS entry
void check_compsegs_consistency(const Opm::WelSegsSet& welsegs,
                                const std::set<std::string>& compsegs,
                                const std::vector<::Opm::Well>& wells)
{
    const auto difference = welsegs.difference(compsegs, wells);
    
    if (!difference.empty()) {
        std::string well_str = "well";
        if (difference.size() > 1) {
            well_str.append("s");
        }
        well_str.append(":");

        for(const auto& [name, location] : difference) {
            well_str.append(fmt::format("\n   {} in {} at line {}",
                                        name, location.filename, location.lineno));
        }
        auto msg = fmt::format("Missing COMPSEGS keyword for the following multisegment {}.", well_str);
        throw Opm::OpmInputError(msg, std::get<1>(difference[0]));
    }
}
}// end anonymous namespace

namespace Opm
{
void Schedule::iterateScheduleSection(std::size_t load_start, std::size_t load_end,
                                      const ParseContext& parseContext,
                                      ErrorGuard& errors,
                                      const ScheduleGrid& grid,
                                      const std::unordered_map<std::string, double> * target_wellpi,
                                      const std::string& prefix,
                                      const bool keepKeywords,
                                      const bool log_to_debug)
{
        std::vector<std::pair< const DeckKeyword* , std::size_t> > rftProperties;
        std::string time_unit = this->m_static.m_unit_system.name(UnitSystem::measure::time);
        auto deck_time = [this](double seconds) { return this->m_static.m_unit_system.from_si(UnitSystem::measure::time, seconds); };
        /*
          The keywords in the skiprest_whitelist set are loaded from the
          SCHEDULE section even though the SKIPREST keyword is in action. The
          full list includes some additional keywords which we do not support at
          all.
        */
        std::unordered_set<std::string> skiprest_whitelist = {"VFPPROD", "VFPINJ", "RPTSCHED", "RPTRST", "TUNING", "MESSAGES"};
        /*
          The behavior of variable restart_skip is more lenient than the
          SKIPREST keyword. If this is a restarted[1] run the loop iterating
          over keywords will skip the all keywords[2] until DATES keyword with
          the restart date is encountered - irrespective of whether the SKIPREST
          keyword is present in the deck or not.

          [1]: opm/flow can restart in a mode where all the keywords from the
               historical part of the Schedule section is internalized, and only
               the solution fields are read from the restart file. In this case
               we will have Schedule::restart_offset() == 0.

          [2]: With the exception of the keywords in the skiprest_whitelist;
               these keywords will be assigned to report step 0.
        */

        const auto restart_skip = load_start < this->m_static.rst_info.report_step;
        ScheduleLogger logger(ScheduleLogger::select_stream(log_to_debug, restart_skip),
                              prefix, this->m_sched_deck.location());
        {
            const auto& location = this->m_sched_deck.location();
            logger({"", "Processing dynamic information from", fmt::format("{} line {}", location.filename, location.lineno)});
            if (restart_skip && !log_to_debug) {
                logger.info(fmt::format("This is a restarted run - skipping "
                                        "until report step {} at {}",
                                        this->m_static.rst_info.report_step,
                                        Schedule::formatDate(this->m_static.rst_info.time)));
            }

            logger(fmt::format("Initializing report step {}/{} at {} {} {} line {}",
                               load_start,
                               this->m_sched_deck.size() - 1,
                               Schedule::formatDate(this->getStartTime()),
                               deck_time(this->m_sched_deck.seconds(load_start)),
                               time_unit,
                               location.lineno));
        }

        std::set<std::string> compsegs_wells;
        WelSegsSet welsegs_wells;

        const auto matches = Action::Result { false }.matches();

        for (auto report_step = load_start; report_step < load_end; report_step++) {
            std::size_t keyword_index = 0;
            auto& block = this->m_sched_deck[report_step];
            auto time_type = block.time_type();
            if (time_type == ScheduleTimeType::DATES || time_type == ScheduleTimeType::TSTEP) {
                const auto& start_date = Schedule::formatDate(std::chrono::system_clock::to_time_t(block.start_time()));
                const auto& days = deck_time(this->stepLength(report_step - 1));
                const auto& days_total = deck_time(this->seconds(report_step - 1));
                logger.complete_step(fmt::format("Complete report step {0} ({1} {2}) at {3} ({4} {2})",
                                                 report_step,
                                                 days,
                                                 time_unit,
                                                 start_date,
                                                 days_total));

                if (report_step < (load_end - 1)) {
                    logger.location(block.location());
                    logger(fmt::format("Initializing report step {}/{} at {} ({} {}) line {}",
                                       report_step + 1,
                                       this->m_sched_deck.size() - 1,
                                       start_date,
                                       days_total,
                                       time_unit,
                                       block.location().lineno));
                }
            }
            this->create_next(block);

            std::unordered_map<std::string, double> wpimult_global_factor;

            while (true) {
                if (keyword_index == block.size())
                    break;

                const auto& keyword = block[keyword_index];
                const auto& location = keyword.location();
                logger.location(location);

                if (keyword.is<ParserKeywords::ACTIONX>()) {
                    auto [action, condition_errors] =
                        Action::parseActionX(keyword,
                                              this->m_static.m_runspec.actdims(),
                                              std::chrono::system_clock::to_time_t(this->snapshots[report_step].start_time()));

                    for(const auto& [ marker, msg]: condition_errors) {
                        parseContext.handleError(marker, msg, keyword.location(), errors);
                    }

                    while (true) {
                        keyword_index++;
                        if (keyword_index == block.size())
                            throw OpmInputError("Missing keyword ENDACTIO", keyword.location());

                        const auto& action_keyword = block[keyword_index];
                        if (action_keyword.is<ParserKeywords::ENDACTIO>())
                            break;

                        bool valid = Action::ActionX::valid_keyword(action_keyword.name());
                        if (this->m_lowActionParsingStrictness or valid){
                            if (this->m_lowActionParsingStrictness and !valid) {
                                logger(fmt::format("The keyword {} is not supported in the ACTIONX block, but you have set --action-parsing-strictness = low, so flow will try to apply the keyword still.", action_keyword.name()));
                            }
                            action.addKeyword(action_keyword);
                            this->prefetchPossibleFutureConnections(grid, action_keyword, parseContext, errors);
                            this->store_wgnames(action_keyword);
                        }
                        else {
                            std::string msg_fmt = fmt::format("The keyword {} is not supported in the ACTIONX block", action_keyword.name());
                            parseContext.handleError( ParseContext::ACTIONX_ILLEGAL_KEYWORD, msg_fmt, action_keyword.location(), errors);
                        }
                    }
                    this->addACTIONX(action);
                    keyword_index++;
                    continue;
                }

                logger(fmt::format("Processing keyword {} at line {}", location.keyword, location.lineno));
                this->handleKeyword(report_step,
                                    block,
                                    keyword,
                                    parseContext,
                                    errors,
                                    grid,
                                    matches,
                                    /* actionx_mode = */ false,
                                    /* sim_update = */ nullptr,
                                    target_wellpi,
                                    wpimult_global_factor,
                                    &welsegs_wells,
                                    &compsegs_wells);
                keyword_index++;
            }

            check_compsegs_consistency(welsegs_wells, compsegs_wells, this->getWells(report_step));
            this->applyGlobalWPIMULT(wpimult_global_factor);
            this->end_report(report_step);

            if (this->must_write_rst_file(report_step)) {
                this->restart_output.addRestartOutput(report_step);
            }

            if (!keepKeywords) {
                this->m_sched_deck.clearKeywords(report_step);
            }
        } // for (auto report_step = load_start
    }

    void Schedule::applyGlobalWPIMULT( const std::unordered_map<std::string, double>& wpimult_global_factor) {
        for (const auto& [well_name, factor] : wpimult_global_factor) {
            auto well = this->snapshots.back().wells(well_name);
            if (well.applyGlobalWPIMULT(factor)) {
                this->snapshots.back().wells.update(std::move(well));
            }
        }
    }

    void Schedule::addACTIONX(const Action::ActionX& action) {
        auto new_actions = this->snapshots.back().actions.get();
        new_actions.add( action );
        this->snapshots.back().actions.update( std::move(new_actions) );
    }


    void Schedule::store_wgnames(const DeckKeyword& keyword) {
        if(keyword.is<ParserKeywords::WELSPECS>()) {
            for (const auto& record : keyword) {
                const auto& wname = record.getItem<ParserKeywords::WELSPECS::WELL>().get<std::string>(0);
                const auto& gname = record.getItem<ParserKeywords::WELSPECS::GROUP>().get<std::string>(0);
                this->action_wgnames.add_well(wname);
                this->action_wgnames.add_group(gname);
            }
        } else if (keyword.is<ParserKeywords::WELOPEN>()  ||
                   keyword.is<ParserKeywords::WCONHIST>() ||
                   keyword.is<ParserKeywords::WCONPROD>() ||
                   keyword.is<ParserKeywords::WCONINJH>() ||
                   keyword.is<ParserKeywords::WCONINJE>()) {      // Add any other keywords that can open a well..

            for (const auto& record : keyword) {
                const std::string& wname_pattern = record.getItem("WELL").getTrimmedString(0);
                this->potential_wellopen_patterns.insert(wname_pattern);
            }
        }
    }


    void Schedule::prefetchPossibleFutureConnections(const ScheduleGrid& grid,
                                                     const DeckKeyword&  keyword,
                                                     const ParseContext& parseContext,
                                                     ErrorGuard&         errors)
    {
        if (keyword.is<ParserKeywords::COMPDAT>()) {
            for (const auto& record : keyword) {
                const auto& itemI = record.getItem("I");
                const auto& itemJ = record.getItem("J");

                const auto defaulted_I = itemI.defaultApplied(0) || (itemI.get<int>(0) == 0);
                const auto defaulted_J = itemJ.defaultApplied(0) || (itemJ.get<int>(0) == 0);

                if (defaulted_I || defaulted_J) {
                    const auto msg_fmt = std::string { R"(Problem with COMPDAT in ACTIONX
In {{file}} line {{line}}
Defaulted grid coordinates is not allowed for COMPDAT as part of ACTIONX)"
                    };

                    parseContext.handleError(ParseContext::SCHEDULE_COMPDAT_INVALID,
                                             msg_fmt, keyword.location(), errors);
                }

                const auto I = itemI.get<int>(0) - 1;
                const auto J = itemJ.get<int>(0) - 1;

                const auto K1 = record.getItem("K1").get<int>(0) - 1;
                const auto K2 = record.getItem("K2").get<int>(0) - 1;

                const auto wellName = record.getItem("WELL").getTrimmedString(0);

                // Retrieve or create the set of future connections for the well
                auto& currentSet = this->possibleFutureConnections[wellName];
                for (int k = K1; k <= K2; k++) {
                    try {
                        // Adds this cell to the "active cells" of the
                        // schedule grid by calling grid.get_cell(I, J, k)
                        const auto& cell = grid.get_cell(I, J, k);

                        // Insert the global id of the cell into the
                        // possible future connections of the well
                        currentSet.insert(cell.global_index);
                    }
                    catch (const std::invalid_argument& e) {
                        const std::string msg_fmt =
                            fmt::format("Problem with COMPDAT in ACTIONX\n"
                                        "In {{file}} line {{line}}\n"
                                        "Cell ({}, {}, {}) of well {} is not part of the grid ({}).",
                                        I + 1, J + 1, k + 1, wellName, e.what());

                        parseContext.handleError(ParseContext::SCHEDULE_COMPDAT_INVALID,
                                                 msg_fmt, keyword.location(), errors);
                    }
                }
            }
        }

        if (keyword.is<ParserKeywords::COMPSEGS>()) {
            bool first_record = true;

            for (const auto& record : keyword) {
                if (first_record) {
                    first_record = false;
                    continue;
                }

                const int I = record.getItem("I").get<int>(0) - 1;
                const int J = record.getItem("J").get<int>(0) - 1;
                const int K = record.getItem("K").get<int>(0) - 1;

                try {
                    const auto& cell = grid.get_cell(I, J, K);
                    static_cast<void>(cell);
                }
                catch (const std::invalid_argument& e) {
                    const auto msg_fmt =
                        fmt::format("Problem with COMPSEGs in ACTIONX\n"
                                    "In {{file}} line {{line}}\n"
                                    "Cell ({}, {}, {}) is not part of the grid ({}).",
                                    I + 1, J + 1, K + 1, e.what());

                    parseContext.handleError(ParseContext::SCHEDULE_COMPSEGS_INVALID,
                                             msg_fmt, keyword.location(), errors);
                }
            }
        }
    }

    void Schedule::shut_well(const std::string& well_name, std::size_t report_step) {
        this->internalWELLSTATUSACTIONXFromPYACTION(well_name, report_step, "SHUT");
    }
    void Schedule::shut_well(const std::string& well_name) {
        this->shut_well(well_name, this->current_report_step);
    }

    void Schedule::open_well(const std::string& well_name, std::size_t report_step) {
        this->internalWELLSTATUSACTIONXFromPYACTION(well_name, report_step, "OPEN");
    }
    void Schedule::open_well(const std::string& well_name) {
        this->open_well(well_name, this->current_report_step);
    }

    void Schedule::stop_well(const std::string& well_name, std::size_t report_step) {
        this->internalWELLSTATUSACTIONXFromPYACTION(well_name, report_step, "STOP");
    }
    void Schedule::stop_well(const std::string& well_name) {
        this->stop_well(well_name, this->current_report_step);
    }

    void Schedule::internalWELLSTATUSACTIONXFromPYACTION(const std::string& well_name, std::size_t report_step, const std::string& wellStatus) {
        if (report_step < this->current_report_step) {
            throw std::invalid_argument(fmt::format("Well status change for past report step {} requested, current report step is {}.", report_step, this->current_report_step));
        } else if (report_step >= this->m_sched_deck.size()) {
            throw std::invalid_argument(fmt::format("Well status change for report step {} requested, this exceeds the total number of report steps, being {}.", report_step, this->m_sched_deck.size() - 1));
        }
        std::time_t start_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::from_time_t(0));
        Opm::Action::ActionX action("openwell", 1, 0.0, start_time);
        DeckItem wellItem("WELL", std::string());
        wellItem.push_back(well_name);
        DeckItem statusItem("STATUS", std::string());
        statusItem.push_back(wellStatus);
        DeckRecord deckRecord;
        deckRecord.addItem(std::move(wellItem));
        deckRecord.addItem(std::move(statusItem));
        ParserKeyword parserKeyword("WELOPEN", KeywordSize(SLASH_TERMINATED));

        DeckKeyword action_keyword(parserKeyword);
        action_keyword.addRecord(std::move(deckRecord));
        action.addKeyword(action_keyword);
        SimulatorUpdate delta = this->applyAction(report_step, action,
                                                  /* matches = */ Action::Result{false}.matches(),
                                                  std::unordered_map<std::string,double>{}/*target_wellpi*/);
        this->simUpdateFromPython->append(delta);
    }

    /*
      Function is quite dangerous - because if this is called while holding a
      Well pointer that will go stale and needs to be refreshed.
    */
    bool Schedule::updateWellStatus( const std::string& well_name, std::size_t reportStep , Well::Status status, std::optional<KeywordLocation> location) {
        if (status != Well::Status::SHUT) {
            this->potential_wellopen_patterns.insert(well_name);
        }
        auto well2 = this->snapshots[reportStep].wells.get(well_name);
        if (well2.getConnections().empty() && status == Well::Status::OPEN) {
            if (location) {
                auto msg = fmt::format("Problem with {}\n"
                                       "In {} line{}\n"
                                       "Well {} has no connections to grid and will remain SHUT", location->keyword, location->filename, location->lineno, well_name);
                OpmLog::warning(msg);
            } else
                OpmLog::warning(fmt::format("Well {} has no connections to grid and will remain SHUT", well_name));
            return false;
        }

        auto old_status = well2.getStatus();
        bool update = false;
        if (well2.updateStatus(status)) {
            if (status == Well::Status::OPEN) {
                auto new_rft = this->snapshots.back().rft_config().well_open(well_name);
                if (new_rft.has_value())
                    this->snapshots.back().rft_config.update( std::move(*new_rft) );
            }

            /*
              The Well::updateStatus() will always return true because a new
              WellStatus object should be created. But the new object might have
              the same value as the previous object; therefor we need to check
              for an actual status change before we emit a WELL_STATUS_CHANGE
              event.
            */
            if (old_status != status) {
                this->snapshots.back().events().addEvent( ScheduleEvents::WELL_STATUS_CHANGE);
                this->snapshots.back().wellgroup_events().addEvent( well2.name(), ScheduleEvents::WELL_STATUS_CHANGE);
            }
            this->snapshots[reportStep].wells.update( std::move(well2) );
            update = true;
        }
        return update;
    }

    void Schedule::clear_event(ScheduleEvents::Events event, std::size_t report_step) {
        auto events = this->snapshots[report_step].events();
        events.clearEvent(event);
        this->snapshots[report_step].update_events(events);
    }

    void Schedule::add_event(ScheduleEvents::Events event, std::size_t report_step)
    {
        auto events = this->snapshots[report_step].events();
        events.addEvent(event);
        this->snapshots[report_step].update_events(events);
    }


    bool Schedule::updateWPAVE(const std::string& wname, std::size_t report_step, const PAvg& pavg) {
        const auto& well = this->getWell(wname, report_step);
        if (well.pavg() != pavg) {
            auto new_well = this->snapshots[report_step].wells.get(wname);
            new_well.updateWPAVE( pavg );
            this->snapshots[report_step].wells.update( std::move(new_well) );
            return true;
        }
        return false;
    }




    std::optional<std::size_t> Schedule::first_RFT() const {
        for (std::size_t report_step = 0; report_step < this->snapshots.size(); report_step++) {
            if (this->snapshots[report_step].rft_config().active())
                return report_step;
        }
        return {};
    }

    GTNode Schedule::groupTree(const std::string& root_node, std::size_t report_step, std::size_t level, const std::optional<std::string>& parent_name) const {
        auto root_group = this->getGroup(root_node, report_step);
        GTNode tree(root_group, level, parent_name);

        for (const auto& wname : root_group.wells()) {
            const auto& well = this->getWell(wname, report_step);
            tree.add_well(well);
        }

        for (const auto& gname : root_group.groups()) {
            auto child_group = this->groupTree(gname, report_step, level + 1, root_node);
            tree.add_group(child_group);
        }

        return tree;
    }

    GTNode Schedule::groupTree(const std::string& root_node, std::size_t report_step) const {
        return this->groupTree(root_node, report_step, 0, {});
    }

    GTNode Schedule::groupTree(std::size_t report_step) const {
        return this->groupTree("FIELD", report_step);
    }

    void Schedule::addWell(const std::string& wellName,
                           const DeckRecord& record,
                           std::size_t timeStep,
                           Connection::Order wellConnectionOrder)
    {
        // We change from eclipse's 1 - n, to a 0 - n-1 solution
        int headI = record.getItem("HEAD_I").get< int >(0) - 1;
        int headJ = record.getItem("HEAD_J").get< int >(0) - 1;
        Phase preferredPhase;
        {
            const std::string phaseStr = record.getItem("PHASE").getTrimmedString(0);
            if (phaseStr == "LIQ") {
                // We need a workaround in case the preferred phase is "LIQ",
                // which is not proper phase and will cause the get_phase()
                // function to throw. In that case we choose to treat it as OIL.
                preferredPhase = Phase::OIL;
                OpmLog::warning("LIQ_PREFERRED_PHASE",
                                "LIQ preferred phase not supported for well " + wellName + ", using OIL instead");
            } else {
                preferredPhase = get_phase(phaseStr);
            }
        }
        const auto& refDepthItem = record.getItem("REF_DEPTH");
        std::optional<double> ref_depth;
        if (refDepthItem.hasValue( 0 ))
            ref_depth = refDepthItem.getSIDouble( 0 );

        double drainageRadius = record.getItem( "D_RADIUS" ).getSIDouble(0);

        bool allowCrossFlow = true;
        const std::string& allowCrossFlowStr = record.getItem<ParserKeywords::WELSPECS::CROSSFLOW>().getTrimmedString(0);
        if (allowCrossFlowStr == "NO")
            allowCrossFlow = false;

        bool automaticShutIn = true;
        const std::string& automaticShutInStr = record.getItem<ParserKeywords::WELSPECS::AUTO_SHUTIN>().getTrimmedString(0);
        if (automaticShutInStr == "STOP") {
            automaticShutIn = false;
        }

        const std::string& group = record.getItem<ParserKeywords::WELSPECS::GROUP>().getTrimmedString(0);
        auto pvt_table = record.getItem<ParserKeywords::WELSPECS::P_TABLE>().get<int>(0);
        auto gas_inflow = WellGasInflowEquationFromString(record.getItem<ParserKeywords::WELSPECS::INFLOW_EQ>().get<std::string>(0));

        this->addWell(wellName,
                      group,
                      headI,
                      headJ,
                      preferredPhase,
                      ref_depth,
                      drainageRadius,
                      allowCrossFlow,
                      automaticShutIn,
                      pvt_table,
                      gas_inflow,
                      timeStep,
                      wellConnectionOrder);
    }

    void Schedule::addWell(Well well) {
        const std::string wname = well.name();
        auto& sched_state = this->snapshots.back();

        sched_state.events().addEvent( ScheduleEvents::NEW_WELL );
        sched_state.wellgroup_events().addWell( wname );
        {
            auto wo = sched_state.well_order.get();
            wo.add( wname );
            sched_state.well_order.update( std::move(wo) );
        }
        well.setInsertIndex(sched_state.wells.size());
        sched_state.wells.update( std::move(well) );
    }

    void Schedule::addWell(const std::string& wellName,
                           const std::string& group,
                           int headI,
                           int headJ,
                           Phase preferredPhase,
                           const std::optional<double>& ref_depth,
                           double drainageRadius,
                           bool allowCrossFlow,
                           bool automaticShutIn,
                           int pvt_table,
                           Well::GasInflowEquation gas_inflow,
                           std::size_t timeStep,
                           Connection::Order wellConnectionOrder) {

        const auto& sched_state = this->operator[](timeStep);
        Well well(wellName,
                  group,
                  timeStep,
                  0,
                  headI, headJ,
                  ref_depth,
                  WellType(preferredPhase),
                  sched_state.whistctl(),
                  wellConnectionOrder,
                  this->m_static.m_unit_system,
                  this->getUDQConfig(timeStep).params().undefinedValue(),
                  drainageRadius,
                  allowCrossFlow,
                  automaticShutIn,
                  pvt_table,
                  gas_inflow,
                  this->m_static.m_runspec.temp());

        this->addWell( std::move(well) );

        const auto& ts = this->operator[](timeStep);
        this->updateWPAVE( wellName, timeStep, ts.pavg.get() );
    }


    std::size_t Schedule::numWells() const {
        return this->snapshots.back().wells.size();
    }

    std::size_t Schedule::numWells(std::size_t timestep) const {
        auto well_names = this->wellNames(timestep);
        return well_names.size();
    }

    bool Schedule::hasWell(const std::string& wellName) const {
        return this->snapshots.back().wells.has(wellName);
    }

    bool Schedule::hasWell(const std::string& wellName, std::size_t timeStep) const {
        return this->snapshots[timeStep].wells.has(wellName);
    }

    bool Schedule::hasGroup(const std::string& groupName, std::size_t timeStep) const {
        return this->snapshots[timeStep].groups.has(groupName);
    }

    /*
      This function will return a list of wells which have changed
      *structurally* in the last report_step; wells where only production
      settings have changed will not be included.
    */
    std::vector<std::string> Schedule::changed_wells(std::size_t report_step) const
    {
        std::vector<std::string> wells;
        const auto& state = this->snapshots[report_step];
        const auto& all_wells = state.wells();

        if (report_step == 0)
            std::transform( all_wells.begin(), all_wells.end(),
                           std::back_inserter(wells),
                           [](const auto& well_ref) { return well_ref.get().name(); });
        else {
            const auto& prev_state = this->snapshots[report_step - 1];
            for (const auto& well_ref : all_wells) {
                const auto& wname = well_ref.get().name();
                if (prev_state.wells.has(wname)) {
                    const auto& prev_well = prev_state.wells.get( wname );
                    if (!prev_well.cmp_structure(well_ref.get()))
                        wells.push_back( wname );
                } else
                    wells.push_back( wname );
            }
        }

        return wells;
    }


    std::vector<Well> Schedule::getWells(std::size_t timeStep) const
    {
        auto wells = std::vector<Well>{};

        if (timeStep >= this->snapshots.size()) {
            throw std::invalid_argument {
                fmt::format("timeStep {} exceeds simulation run's "
                            "number of report steps ({})",
                            timeStep, this->snapshots.size())
            };
        }

        const auto& well_order = this->snapshots[timeStep].well_order();
        std::transform(well_order.begin(), well_order.end(),
                       std::back_inserter(wells),
                       [&wells = this->snapshots[timeStep].wells]
                       (const auto& wname) -> decltype(auto)
                       { return wells.get(wname); });

        return wells;
    }

    std::vector<Well> Schedule::getWellsatEnd() const {
        return this->getWells(this->snapshots.size() - 1);
    }

    std::vector<Well> Schedule::getActiveWellsAtEnd() const {
        std::vector<Well> wells;
        const auto lastStep = this->snapshots.size() - 1;
        const auto& well_order = this->snapshots[lastStep].well_order();

        for (const auto& wname : well_order) {
            const auto& well = this->snapshots[lastStep].wells.get(wname);
            if (well.hasProduced() || well.hasInjected() || name_match_any(this->potential_wellopen_patterns, wname))
                wells.push_back(well);
        }

        return wells;
    }

    std::vector<std::string> Schedule::getInactiveWellNamesAtEnd() const {
        std::vector<std::string> well_names;
        const auto lastStep = this->snapshots.size() - 1;
        const auto& well_order = this->snapshots[lastStep].well_order();

        for (const auto& wname : well_order) {
            const auto& well = this->snapshots[lastStep].wells.get(wname);
            if (well.hasProduced() || well.hasInjected() || name_match_any(this->potential_wellopen_patterns, wname))
                continue;
            well_names.push_back(wname);
        }

        return well_names;
    }


    const Well& Schedule::getWellatEnd(const std::string& well_name) const {
        return this->getWell(well_name, this->snapshots.size() - 1);
    }

    const std::unordered_map<std::string, std::set<int>>&
    Schedule::getPossibleFutureConnections() const
    {
        return this->possibleFutureConnections;
    }

    std::unordered_set<int> Schedule::getAquiferFluxSchedule() const {
        std::unordered_set<int> ids;
        for (const auto& snapshot : this->snapshots) {
            const auto& aquflux = snapshot.aqufluxs;
            for ([[maybe_unused]] const auto& [id, aqu]  : aquflux) {
                ids.insert(id);
            }
        }
        return ids;
    }

    const Well& Schedule::getWell(const std::string& wellName, std::size_t timeStep) const {
        return this->snapshots[timeStep].wells.get(wellName);
    }

    const Well& Schedule::getWell(std::size_t well_index, std::size_t timeStep) const {
        const auto find_pred = [well_index] (const auto& well_pair) -> bool
        {
            return well_pair.second->seqIndex() == well_index;
        };

        auto well_ptr = this->snapshots[timeStep].wells.find( find_pred );
        if (well_ptr == nullptr)
            throw std::invalid_argument(fmt::format("There is no well with well_index:{} at report_step:{}", well_index, timeStep));

        return *well_ptr;
    }

    const Group& Schedule::getGroup(const std::string& groupName, std::size_t timeStep) const {
        return this->snapshots[timeStep].groups.get(groupName);
    }

    void Schedule::updateGuideRateModel(const GuideRateModel& new_model, std::size_t report_step) {
        auto new_config = this->snapshots[report_step].guide_rate();
        if (new_config.update_model(new_model))
            this->snapshots[report_step].guide_rate.update( std::move(new_config) );
    }

    // There are many SCHEDULE keywords which operate on well names.  In
    // addition to fully qualified names like 'W1', there are shell-style
    // wildcard patterns like 'W*' or 'PROD?'.  Similarly, you can request
    // all wells in a well list '*WL'[1] and the well name '?', when used in
    // an ACTIONX keyword block, matches all wells which trigger the
    // condition in the same ACTIONX keyword[2].  This function is intended
    // to be the final arbiter for well names matching these kinds of
    // patterns.  The time step argument filters out wells which do not
    // exist at that time level (i.e., zero-based report step index).
    //
    // [1]: The leading '*' in a WLIST name should not be treated as a
    //      pattern matching wildcard.  On the other hand, the pattern
    //      '\*WL' matches all wells whose names end in 'WL'.  In this case,
    //      the leading backslash "escapes" the initial asterisk, thus
    //      disambiguating it as a normal wildcard.
    //
    // [2]: A leading '?' character can be escaped as '\?' in order to not
    //      be misconstrued as the '?' pattern.  Thus, the pattern '\?????'
    //      matches all wells whose names consist of exactly five
    //      characters.
    std::vector<std::string>
    Schedule::wellNames(const std::string&              pattern,
                        const std::size_t               timeStep,
                        const std::vector<std::string>& matching_wells) const
    {
        const auto wm = this->wellMatcher(timeStep);

        return (pattern == "?")
            ? wm.sort(matching_wells) // ACTIONX handler
            : wm.wells(pattern);      // Normal well name pattern matching
    }



    std::vector<std::string>
    Schedule::wellNames(const std::string&    pattern,
                        const HandlerContext& context,
                        const bool            allowEmpty)
    {
        auto names = this->wellNames(pattern, context.currentStep,
                                     context.matches.wells().asVector());

        if (names.empty() && !allowEmpty) {
            if (this->action_wgnames.has_well(pattern)) {
                const auto& location = context.keyword.location();

                const auto msg = fmt::format(R"(Well: {} not yet defined for keyword {}.
Expecting well to be defined with WELSPECS in ACTIONX before actual use.
File {} line {}.)", pattern, location.keyword, location.filename, location.lineno);

                OpmLog::warning(msg);
            }
            else {
                context.invalidNamePattern(pattern);
            }
        }

        return names;
    }

    WellMatcher Schedule::wellMatcher(const std::size_t report_step) const
    {
        const auto& schedState = (report_step < this->snapshots.size())
            ? this->snapshots[report_step]
            : this->snapshots.back();

        return { &schedState.well_order(), schedState.wlist_manager() };
    }

    std::function<std::unique_ptr<SegmentMatcher>()>
    Schedule::segmentMatcherFactory(const std::size_t report_step) const
    {
        return {
            [report_step, this]() -> std::unique_ptr<SegmentMatcher>
            {
                return std::make_unique<SegmentMatcher>((*this)[report_step]);
            }
        };
    }

    std::vector<std::string> Schedule::wellNames(const std::string& pattern) const
    {
        return this->wellNames(pattern, this->size() - 1);
    }

    std::vector<std::string> Schedule::wellNames(std::size_t timeStep) const
    {
        return this->snapshots[timeStep].well_order().names();
    }

    std::vector<std::string> Schedule::wellNames() const
    {
        return this->snapshots.back().well_order().names();
    }

    std::vector<std::string> Schedule::groupNames(const std::string& pattern,
                                                  const std::size_t timeStep) const
    {
        return this->snapshots[timeStep].group_order().names(pattern);
    }

    const std::vector<std::string>& Schedule::groupNames(std::size_t timeStep) const
    {
        return this->snapshots[timeStep].group_order().names();
    }

    std::vector<std::string> Schedule::groupNames(const std::string& pattern) const
    {
        return this->groupNames(pattern, this->snapshots.size() - 1);
    }

    const std::vector<std::string>& Schedule::groupNames() const
    {
        return this->snapshots.back().group_order().names();
    }

    std::vector<const Group*> Schedule::restart_groups(std::size_t timeStep) const
    {
        const auto restart_groups = this->snapshots[timeStep].group_order().restart_groups();

        std::vector<const Group*> rst_groups(restart_groups.size(), nullptr);
        for (std::size_t restart_index = 0;
             restart_index < restart_groups.size(); ++restart_index)
        {
            const auto& group_name = restart_groups[restart_index];

            if (group_name.has_value()) {
                rst_groups[restart_index] = &this->getGroup(group_name.value(), timeStep);
            }
        }

        return rst_groups;
    }


    void Schedule::addGroup(Group group) {
        std::string group_name = group.name();
        auto& sched_state = this->snapshots.back();
        sched_state.groups.update(std::move(group) );
        sched_state.events().addEvent( ScheduleEvents::NEW_GROUP );
        sched_state.wellgroup_events().addGroup(group_name);
        {
            auto go = sched_state.group_order.get();
            go.add( group_name );
            sched_state.group_order.update( std::move(go) );
        }

        // All newly created groups are attached to the field group,
        // can then be relocated with the GRUPTREE keyword.
        if (group_name != "FIELD")
            this->addGroupToGroup("FIELD", group_name);
    }


    void Schedule::addGroup(const std::string& groupName, std::size_t timeStep) {
        auto udq_undefined = this->getUDQConfig(timeStep).params().undefinedValue();
        const auto& sched_state = this->snapshots.back();
        auto insert_index = sched_state.groups.size();
        this->addGroup( Group(groupName, insert_index, udq_undefined, this->m_static.m_unit_system) );
    }


    void Schedule::addGroup(const RestartIO::RstGroup& rst_group, std::size_t timeStep) {
        auto udq_undefined = this->getUDQConfig(timeStep).params().undefinedValue();

        const auto insert_index = this->snapshots.back().groups.size();
        auto new_group = Group(rst_group, insert_index, udq_undefined, this->m_static.m_unit_system);
        if (rst_group.name != "FIELD") {
            // Common case.  Add new group.
            this->addGroup( std::move(new_group) );
            return;
        }

        // If we get here we're updating the FIELD group to incorporate any
        // applicable field-wide GCONPROD and/or GCONINJE settings stored in
        // the restart file.  Happens at most once per run.

        auto& field = this->snapshots.back().groups.get("FIELD");
        if (new_group.isProductionGroup())
            // Initialise field-wide GCONPROD settings from restart.
            field.updateProduction(new_group.productionProperties());
        for (const auto phase : { Phase::GAS, Phase::WATER })
            if (new_group.hasInjectionControl(phase))
                // Initialise field-wide GCONINJE settings (phase) from restart.
                field.updateInjection(new_group.injectionProperties(phase));
    }


    void Schedule::addGroupToGroup( const std::string& parent_name, const std::string& child_name) {
        auto parent_group = this->snapshots.back().groups.get(parent_name);
        if (parent_group.addGroup(child_name))
            this->snapshots.back().groups.update( std::move(parent_group) );

        // Check and update backreference in child
        const auto& child_group = this->snapshots.back().groups.get(child_name);
        if (child_group.parent() != parent_name) {
            auto old_parent = this->snapshots.back().groups.get(child_group.parent());
            old_parent.delGroup(child_group.name());
            this->snapshots.back().groups.update( std::move(old_parent) );

            auto new_child_group = Group{ child_group };
            new_child_group.updateParent(parent_name);
            this->snapshots.back().groups.update( std::move(new_child_group) );
        }

        // Update standard network if required
        auto network = this->snapshots.back().network.get();
        if (!network.is_standard_network())
            return;
        if (network.has_node(child_name)) {
            auto old_branch = network.uptree_branch(child_name);
            if (old_branch.has_value()) {
                auto new_branch = old_branch.value();
                new_branch.set_uptree_node(parent_name);
                network.add_or_replace_branch(new_branch);
                this->snapshots.back().network.update( std::move(network));
            }
            // If no previous uptree branch the child is a fixed-pressure node, so no need to update network
        }
    }

    void Schedule::addWellToGroup( const std::string& group_name, const std::string& well_name , std::size_t timeStep) {
        auto well = this->getWell(well_name, timeStep);
        const auto old_gname = well.groupName();
        if (old_gname != group_name) {
            well.updateGroup(group_name);
            this->snapshots.back().wells.update( std::move(well) );
            this->snapshots.back().wellgroup_events().addEvent( well_name, ScheduleEvents::WELL_WELSPECS_UPDATE );

            // Remove well child reference from previous group
            auto group = this->snapshots.back().groups.get( old_gname );
            group.delWell(well_name);
            this->snapshots.back().groups.update( std::move(group) );
        }

        // Add well child reference to new group
        auto group = this->snapshots.back().groups.get( group_name );
        group.addWell(well_name);
        this->snapshots.back().groups.update( std::move(group) );
        this->snapshots.back().events().addEvent( ScheduleEvents::GROUP_CHANGE );
    }



    Well::ProducerCMode Schedule::getGlobalWhistctlMmode(std::size_t timestep) const {
        return this->operator[](timestep).whistctl();
    }



    void Schedule::checkIfAllConnectionsIsShut(std::size_t timeStep) {
        const auto& well_names = this->wellNames(timeStep);
        for (const auto& wname : well_names) {
            const auto& well = this->getWell(wname, timeStep);
            const auto& connections = well.getConnections();
            if (connections.allConnectionsShut() && well.getStatus() != Well::Status::SHUT) {
                auto days = unit::convert::to(seconds(timeStep), unit::day);
                auto msg = fmt::format("All completions in well {} is shut at {} days\n"
                                       "The well is therefore also shut", well.name(), days);
                OpmLog::note(msg);
                this->updateWellStatus( well.name(), timeStep, Well::Status::SHUT);
            }
        }
    }

    void Schedule::end_report(std::size_t report_step) {
        this->checkIfAllConnectionsIsShut(report_step);
    }


    void Schedule::filterConnections(const ActiveGridCells& grid) {
        for (auto& sched_state : this->snapshots) {
            for (auto& well : sched_state.wells()) {
                well.get().filterConnections(grid);
            }
        }
    }


    const UDQConfig& Schedule::getUDQConfig(std::size_t timeStep) const {
        return this->snapshots[timeStep].udq.get();
    }

    std::optional<int> Schedule::exitStatus() const {
        return this->exit_status;
    }

    std::size_t Schedule::size() const {
        return this->snapshots.size();
    }


    double Schedule::seconds(std::size_t timeStep) const {
        if (this->snapshots.empty())
            return 0;

        if (timeStep >= this->snapshots.size())
            throw std::logic_error(fmt::format("seconds({}) - invalid timeStep. Valid range [0,{}>", timeStep, this->snapshots.size()));

        auto elapsed = this->snapshots[timeStep].start_time() - this->snapshots[0].start_time();
        using DurationInSeconds = std::chrono::duration<double>; // Tick is 1 second, stored in double.
        return DurationInSeconds(elapsed).count();
    }

    std::time_t Schedule::simTime(std::size_t timeStep) const {
        return std::chrono::system_clock::to_time_t( this->snapshots[timeStep].start_time() );
    }

    double Schedule::stepLength(std::size_t timeStep) const {
        const auto start_time = this->snapshots[timeStep].start_time();
        const auto end_time = this->snapshots[timeStep].end_time();
        if (start_time > end_time) {
            throw std::invalid_argument {
                    fmt::format(" Report step {} has start time after end time,\n"
                                "   * Start time = {:%d-%b-%Y %H:%M:%S}\n"
                                "   * End time   = {:%d-%b-%Y %H:%M:%S}.\n"
                                " Possibly due to inconsistent RESTART/SKIPREST settings.",
                                timeStep + 1,
                                fmt::gmtime(TimeService::to_time_t(start_time)),
                                fmt::gmtime(TimeService::to_time_t(end_time))) };
        }
        using DurationInSeconds = std::chrono::duration<double>; // Tick is 1 second, stored in double.
        return DurationInSeconds(end_time - start_time).count();
    }

    void Schedule::applyKeywords(std::vector<std::unique_ptr<DeckKeyword>>& keywords, std::unordered_map<std::string, double>& target_wellpi, bool action_mode)
    {
        Schedule::applyKeywords(keywords, target_wellpi, action_mode, this->current_report_step);
    }

    void Schedule::applyKeywords(std::vector<std::unique_ptr<DeckKeyword>>& keywords, std::unordered_map<std::string, double>& target_wellpi,
                                 bool action_mode, const std::size_t reportStep)
    {
        if (reportStep < this->current_report_step) {
            throw std::invalid_argument {
                fmt::format("Insert keyword for past report step {} "
                            "requested, current report step is {}.",
                            reportStep, this->current_report_step)
            };
        }
        else if (reportStep >= this->m_sched_deck.size()) {
            throw std::invalid_argument {
                fmt::format("Insert keyword for report step {} requested "
                            "which exceeds the total number of report steps {}.",
                            reportStep, this->m_sched_deck.size() - 1)
            };
        }

        ParseContext parseContext;
        ErrorGuard errors;
        ScheduleGrid grid(this->completed_cells);
        SimulatorUpdate sim_update;
        std::unordered_map<std::string, double> wpimult_global_factor;
        const auto matches = Action::Result{false}.matches();
        const std::string prefix = "| "; // logger prefix string

        this->snapshots.resize(reportStep + 1);

        auto& input_block = this->m_sched_deck[reportStep];
        ScheduleLogger logger(ScheduleLogger::select_stream(true, false), // will log to OpmLog::debug
                              prefix, this->m_sched_deck.location());

        for (const auto& keyword : keywords) {
            const auto valid = Action::PyAction::valid_keyword(keyword->name());

            if (valid || this->m_lowActionParsingStrictness) {
                if (!valid) {
                    logger(fmt::format("The keyword {} is not supported for insertion "
                                       "from Python into a simulation, but you have set "
                                       "--action-parsing-strictness = low, so flow will "
                                       "try to apply the keyword still.", keyword->name()));
                }

                input_block.push_back(*keyword);

                this->handleKeyword(reportStep,
                                    input_block,
                                    *keyword,
                                    parseContext,
                                    errors,
                                    grid,
                                    matches,
                                    action_mode,
                                    &sim_update,
                                    &target_wellpi,
                                    wpimult_global_factor);    
            }
            else {
                const std::string msg_fmt =
                    fmt::format("The keyword {} is not supported for insertion "
                                "from Python into a simulation", keyword->name());

                parseContext.handleError(ParseContext::PYACTION_ILLEGAL_KEYWORD,
                                         msg_fmt, keyword->location(), errors);
            }
        }

        this->applyGlobalWPIMULT(wpimult_global_factor);
        this->end_report(reportStep);

        if (reportStep < this->m_sched_deck.size() - 1) {
            this->iterateScheduleSection(reportStep + 1,
                                         this->m_sched_deck.size(),
                                         parseContext,
                                         errors,
                                         grid,
                                         &target_wellpi,
                                         prefix,
                                         /* keepKeywords = */ true,
                                         /* log_to_debug */ true);
        }

        this->simUpdateFromPython->append(sim_update);
    }

    namespace {

        std::unordered_map<std::string, double>
        convertToDoubleMap(const std::unordered_map<std::string, float>& target_wellpi)
        {
            return { target_wellpi.begin(), target_wellpi.end() };
        }

    } // Anonymous namespace

    SimulatorUpdate
    Schedule::applyAction(const std::size_t reportStep,
                          const Action::ActionX& action,
                          const Action::Result::MatchingEntities& matches,
                          const std::unordered_map<std::string, float>& target_wellpi)
    {
        return this->applyAction(reportStep, action, matches, convertToDoubleMap(target_wellpi));
    }


    SimulatorUpdate
    Schedule::applyAction(std::size_t reportStep,
                          const Action::ActionX& action,
                          const Action::Result::MatchingEntities& matches,
                          const std::unordered_map<std::string, double>& target_wellpi)
    {
        const std::string prefix = "| ";
        ParseContext parseContext;
        if (this->m_treat_critical_as_non_critical) { // Continue with invalid names if parsing strictness is set to low
            parseContext.update(ParseContext::SCHEDULE_INVALID_NAME, InputErrorAction::WARN);
        }

        ErrorGuard errors;
        SimulatorUpdate sim_update;
        ScheduleGrid grid(this->completed_cells);

        OpmLog::debug("/----------------------------------------------------------------------");
        OpmLog::debug(fmt::format("{0}Action {1} triggered. Will add action "
                                  "keywords and\n{0}rerun Schedule section.\n{0}",
                                  prefix, action.name()));

        this->snapshots.resize(reportStep + 1);
        auto& input_block = this->m_sched_deck[reportStep];

        std::unordered_map<std::string, double> wpimult_global_factor;
        for (const auto& keyword : action) {
            input_block.push_back(keyword);

            const auto& location = keyword.location();
            OpmLog::debug(fmt::format("{}Processing keyword {} from {} line {}", prefix,
                                      location.keyword, location.filename, location.lineno));

            this->handleKeyword(reportStep,
                                input_block,
                                keyword,
                                parseContext,
                                errors,
                                grid,
                                matches,
                                /* actionx_mode = */ true,
                                &sim_update,
                                &target_wellpi,
                                wpimult_global_factor);
        }

        this->applyGlobalWPIMULT(wpimult_global_factor);

        this->end_report(reportStep);

        if (! sim_update.affected_wells.empty()) {
            this->snapshots.back().events()
                .addEvent(ScheduleEvents::ACTIONX_WELL_EVENT);

            auto& wgEvents = this->snapshots.back().wellgroup_events();

            for (const auto& well: sim_update.affected_wells) {
                wgEvents.addEvent(well, ScheduleEvents::ACTIONX_WELL_EVENT);
            }
        }

        if (reportStep < this->m_sched_deck.size() - 1) {
            const auto keepKeywords = true;
            const auto log_to_debug = true;
            this->iterateScheduleSection(reportStep + 1, this->m_sched_deck.size(),
                                         parseContext, errors, grid, &target_wellpi,
                                         prefix, keepKeywords, log_to_debug);
        }

        OpmLog::debug("\\----------------------------------------------------------------------");

        return sim_update;
    }

    SimulatorUpdate
    Schedule::modifyCompletions(const std::size_t reportStep,
                                const std::map<std::string, std::vector<Connection>>& extraConns)
    {
        SimulatorUpdate sim_update{};

        this->snapshots.resize(reportStep + 1);
        for (const auto& [well, newConns] : extraConns) {
            if (newConns.empty()) { continue; }

            // Note: We go through map_member::get() here rather than
            // map_member::operator() because the latter returns a reference
            // to const whereas we need a reference to mutable.
            auto& conns = this->snapshots[reportStep]
                .wells.get(well).getConnections();

            auto allConnsExist = true;
            for (const auto& newConn : newConns) {
                auto* existingConn = conns
                    .maybeGetFromGlobalIndex(newConn.global_index());

                if (existingConn != nullptr) {
                    // Connection 'newConn' already exists in 'conns'.
                    // Change existing CTF if needed.
                    if (newConn.CF() > existingConn->CF()) {
                        existingConn->setCF(newConn.CF());
                    }
                }
                else {
                    // 'newConn' does not already exist in 'conns'.  Add to
                    // collection.
                    allConnsExist = false;

                    const auto seqIndex = conns.size();
                    conns.addConnection(newConn.getI(),
                                        newConn.getJ(),
                                        newConn.getK(),
                                        newConn.global_index(),
                                        newConn.state(),
                                        newConn.depth(),
                                        newConn.ctfProperties(),
                                        /* satTableID = */ 1,
                                        newConn.dir(),
                                        newConn.kind(),
                                        seqIndex,
                                        /* defaultSatTableID = */ false);
                }
            }

            if (allConnsExist) {
                sim_update.welpi_wells.insert(well);
            }
            else {
                sim_update.well_structure_changed = true;
            }
        }

        if (reportStep < this->m_sched_deck.size() - 1) {
            ParseContext parseContext{};
            if (this->m_treat_critical_as_non_critical) {
                // Continue with invalid names if parsing strictness is set
                // to low.
                parseContext.update(ParseContext::SCHEDULE_INVALID_NAME,
                                    InputErrorAction::WARN);
            }

            ErrorGuard errors{};
            ScheduleGrid grid(this->completed_cells);

            const std::string prefix = "| "; /* logger prefix string */

            const auto keepKeywords = true;
            const auto log_to_debug = true;
            this->iterateScheduleSection(reportStep + 1, this->m_sched_deck.size(),
                                         parseContext, errors, grid,
                                         /* target_wellpi = */ nullptr,
                                         prefix, keepKeywords, log_to_debug);
        }

        return sim_update;
    }


    // This function will typically be called from the apply_action_callback()
    // which is invoked in a PYACTION plugin, i.e. the arguments here are
    // supplied by the user in a script - can very well be wrong.
    SimulatorUpdate Schedule::applyAction(const std::size_t reportStep,
                                          const std::string& action_name,
                                          const std::vector<std::string>& matching_wells)
    {
        const auto& actions = this->snapshots[reportStep].actions();
        if (actions.has(action_name)) {
            std::vector<std::string> well_names;
            for (const auto& wname : matching_wells) {
                if (this->hasWell(wname, reportStep)) {
                    well_names.push_back(wname);
                }
                else {
                    OpmLog::error(fmt::format("Tried to apply action {} on non-"
                                              "existing well '{}'", action_name, wname));
                }
            }

            return this->applyAction(reportStep, actions[action_name],
                                     Action::Result{true}.wells(well_names).matches(),
                                     std::unordered_map<std::string,double>{});
        }
        else {
            OpmLog::error(fmt::format("Tried to apply unknown action: '{}'", action_name));
            return {};
        }
    }


    /*
        The runPyAction() method is a utility to run PYACTION keywords. The
        PYACTION keywords contain a link to a file with Python code that will be executed, as
        documented in https://opm.github.io/opm-python-documentation/master/index.html.

        For backwards compatibility, we have kept the *old* way of using the PyAction keyword,
        where the Python code needs to contain a run function with the signature
        def run(ecl_state, schedule, report_step, summary_state, actionx_callback).

        The ecl_state, schedule, report_step and summary_state objects can be accessed as
        documented in https://opm.github.io/opm-python-documentation/master/index.html.

        The lambda 'apply_action_callback' has been kept for backwards compatibility and
        can be used from Python to call back to the C++ code in order to run ACTIONX keywords
        defined in the .DATA file. The sequence of calls is:

        1. The simulator calls the method Schedule::runPyAction()

        2. The Schedule::runPyAction() method creates a SimulatorUpdate
            instance and captures a reference to that in the lambda
            apply_action_callback. When calling the pyaction.run() method the
            apply_action_callback is passed as a callable all the way down to
            the python run() method.

        3. In python the apply_action_callback comes in as the parameter
            'actionx_callback' in the run() function. If the python code decides
            to invoke the keywords from an actionx it will be like:

            def run(ecl_state, schedule, report_step, summary_state, actionx_callback):
                ...
                ...
                wells = ["W1", "W2"]
                actionx_callback("ACTION_NAME", wells)

            Observe that the wells argument must be a Python lvalue!

        4. The callable will go back into C++ and eventually reach the
            Schedule::applyAction() which will invoke the
            Schedule::iterateScheduleSection() and return an update
            SimulatorUpdate which will be assigned to the instance in the
            Schedule::runPyAction().

        5. When the pyaction.run() method returns the Schedule structure and
            the sim_update variable have been correctly updated, and the
            sim_update is returned to the simulator.

        For the apply_action_callback to work, three different systems must be aligned:

        1. The Python code executes normally, and will possibly decide to
            apply an ACTIONX keyword.

        2. When an AXTIONX keyword is applied the Schedule implementation will
            need to add the new keywords to the correct ScheduleBlock and
            reiterate the Schedule section.

        3. As part of the Schedule iteration we record which changes which must
            be taken into account in the simulator afterwards. These changes are
            recorded in a Action::SimulatorUpdate instance.
    */

    SimulatorUpdate Schedule::runPyAction(std::size_t reportStep,
                                          const Action::PyAction& pyaction,
                                          Action::State& action_state,
                                          EclipseState& ecl_state,
                                          SummaryState& summary_state) {
        return this->runPyAction(reportStep, pyaction, action_state, ecl_state, summary_state, std::unordered_map<std::string, double>{});
    }

    SimulatorUpdate Schedule::runPyAction(std::size_t reportStep,
                                          const Action::PyAction& pyaction,
                                          Action::State& action_state,
                                          EclipseState& ecl_state,
                                          SummaryState& summary_state,
                                          const std::unordered_map<std::string, float>& target_wellpi) {
        return this->runPyAction(reportStep, pyaction, action_state, ecl_state, summary_state, convertToDoubleMap(target_wellpi));
    }

    SimulatorUpdate Schedule::runPyAction(std::size_t reportStep,
                                          const Action::PyAction& pyaction,
                                          Action::State& action_state,
                                          EclipseState& ecl_state,
                                          SummaryState& summary_state,
                                          const std::unordered_map<std::string, double>& target_wellpi)
    {
        // Reset simUpdateFromPython, pyaction.run(...) will run through the PyAction script, the calls that trigger a simulator update will append this to simUpdateFromPython.
        this->simUpdateFromPython->reset();
        // Set the current_report_step to the report step in which this PyAction was triggered.
        this->current_report_step = reportStep;

        // Set up the actionx_callback - simulator updates from this also get appended to simUpdateFromPython.
        auto apply_action_callback = [&reportStep, this](const std::string& action_name, const std::vector<std::string>& matching_wells) {
            SimulatorUpdate simUpdateFromActionXCallback = this->applyAction(reportStep, action_name, matching_wells);
            this->simUpdateFromPython->append(simUpdateFromActionXCallback);
        };

        auto result = pyaction.run(ecl_state, *this, reportStep, summary_state, apply_action_callback, target_wellpi);
        action_state.add_run(pyaction, result);

        // The whole pyaction script was executed, now the simUpdateFromPython is returned.
        return *(this->simUpdateFromPython);
    }

    void Schedule::applyWellProdIndexScaling(const std::string& well_name, const std::size_t reportStep, const double newWellPI) {
        if (reportStep >= this->snapshots.size())
            return;

        if (!this->snapshots[reportStep].wells.has(well_name))
            return;

        std::vector<Well *> unique_wells;
        for (std::size_t step = reportStep; step < this->snapshots.size(); step++) {
            auto& well = this->snapshots[step].wells.get(well_name);
            if (unique_wells.empty() || (!(*unique_wells.back() == well)))
                unique_wells.push_back( &well );
        }

        std::vector<bool> scalingApplicable;
        const auto targetPI = this->snapshots[reportStep].target_wellpi.at(well_name);
        auto prev_well = unique_wells[0];
        auto scalingFactor = prev_well->convertDeckPI(targetPI) / newWellPI;
        prev_well->applyWellProdIndexScaling(scalingFactor, scalingApplicable);

        for (std::size_t well_index = 1; well_index < unique_wells.size(); well_index++) {
            auto wellPtr = unique_wells[well_index];
            if (! wellPtr->hasSameConnectionsPointers(*prev_well)) {
                wellPtr->applyWellProdIndexScaling(scalingFactor, scalingApplicable);
                prev_well = wellPtr;
            }
        }
    }

    bool Schedule::write_rst_file(const std::size_t report_step) const
    {
        return this->restart_output.writeRestartFile(report_step) || this->operator[](report_step).save();
    }

    bool Schedule::must_write_rst_file(const std::size_t report_step) const
    {
        if (this->m_static.output_interval.has_value() &&
            (*this->m_static.output_interval > 0) &&
            (report_step > 0))
        {
            return (report_step % *this->m_static.output_interval) == 0;
        }

        if (report_step == 0) {
            return this->m_static.rst_config.write_rst_file.has_value()
                && *this->m_static.rst_config.write_rst_file;
        }

        const auto previous_restart_output_step =
            this->restart_output.lastRestartEventBefore(report_step);

        // Previous output event time or start of simulation if no previous
        // event recorded
        const auto previous_output = previous_restart_output_step.has_value()
            ? this->snapshots[previous_restart_output_step.value()].start_time()
            : this->snapshots[0].start_time();

        const auto& rst_config = this->snapshots[report_step - 1].rst_config();
        return this->snapshots[report_step].rst_file(rst_config, previous_output);
    }

    bool Schedule::isWList(std::size_t report_step, const std::string& pattern) const
    {
        const ScheduleState * sched_state;

        if (report_step < this->snapshots.size())
            sched_state = &this->snapshots[report_step];
        else
            sched_state = &this->snapshots.back();

        return sched_state->wlist_manager.get().hasList(pattern);
    }

    const std::map< std::string, int >& Schedule::rst_keywords( size_t report_step ) const {
        if (report_step == 0)
            return this->m_static.rst_config.keywords;

        const auto& keywords = this->snapshots[report_step - 1].rst_config().keywords;
        return keywords;
    }

    bool Schedule::operator==(const Schedule& data) const {
        // If this has a simUpdateFromPython pointer and data does not
        // (or the other way round), then they are *not* equal.
        if ((this->simUpdateFromPython && !data.simUpdateFromPython) ||
            (!this->simUpdateFromPython && data.simUpdateFromPython))
            return false;

        bool simUpdateFromPythonIsEqual = !this->simUpdateFromPython ||
             (*(this->simUpdateFromPython) == *(data.simUpdateFromPython));

        return this->m_static == data.m_static
            && this->m_treat_critical_as_non_critical == data.m_treat_critical_as_non_critical
            && this->m_sched_deck == data.m_sched_deck
            && this->action_wgnames == data.action_wgnames
            && this->potential_wellopen_patterns == data.potential_wellopen_patterns
            && this->exit_status == data.exit_status
            && this->snapshots == data.snapshots
            && this->restart_output == data.restart_output
            && this->completed_cells == data.completed_cells
            && this->current_report_step == data.current_report_step
            && this->m_lowActionParsingStrictness == data.m_lowActionParsingStrictness
            && simUpdateFromPythonIsEqual
            ;
     }


    std::string Schedule::formatDate(std::time_t t) {
        const auto ts { TimeStampUTC(t) } ;
        return fmt::format("{:04d}-{:02d}-{:02d}" , ts.year(), ts.month(), ts.day());
    }

    std::string Schedule::simulationDays(std::size_t currentStep) const {
        const double sim_time { this->m_static.m_unit_system.from_si(UnitSystem::measure::time, simTime(currentStep)) } ;
        return fmt::format("{} {}", sim_time, this->m_static.m_unit_system.name(UnitSystem::measure::time));
    }

namespace {

    // Duplicated from Well.cpp
    Connection::Order order_from_int(int int_value) {
        switch(int_value) {
        case 0:
            return Connection::Order::TRACK;
        case 1:
            return Connection::Order::DEPTH;
        case 2:
            return Connection::Order::INPUT;
        default:
            throw std::invalid_argument("Invalid integer value: " + std::to_string(int_value) + " encountered when determining connection ordering");
        }
    }
}

    void Schedule::load_rst(const RestartIO::RstState& rst_state,
                            const TracerConfig&        tracer_config,
                            const ScheduleGrid&        grid,
                            const FieldPropsManager&   fp)
    {
        const auto report_step = rst_state.header.report_step - 1;

        std::map<int, std::string> rst_group_names;
        for (const auto& rst_group : rst_state.groups) {
            this->addGroup(rst_group, report_step);
            const auto& group = this->snapshots.back().groups.get( rst_group.name );
            rst_group_names[group.insert_index()] = rst_group.name;
            if (group.isProductionGroup()) {
                // Was originally at report_step + 1
                this->snapshots.back().events().addEvent(ScheduleEvents::GROUP_PRODUCTION_UPDATE );
                this->snapshots.back().wellgroup_events().addEvent(rst_group.name, ScheduleEvents::GROUP_PRODUCTION_UPDATE);
            }

            if (group.isInjectionGroup()) {
                // Was originally at report_step + 1
                this->snapshots.back().events().addEvent(ScheduleEvents::GROUP_INJECTION_UPDATE );
                this->snapshots.back().wellgroup_events().addEvent(rst_group.name, ScheduleEvents::GROUP_INJECTION_UPDATE);
            }

            OpmLog::info(fmt::format("Adding group {} from restart file", rst_group.name));
        }

        //! \todo{ Restart GCONSUMP when consumption/import is defined via UDQs. }
        //! \todo{ Restart GCONSUMP with network node name defined. }
        auto udq_undefined = this->getUDQConfig(report_step).params().undefinedValue();
        auto new_gconsump = this->snapshots.back().gconsump.get();
        for (const auto& rst_group : rst_state.groups) {
            const auto crate = rst_group.gas_consumption_rate;
            const auto irate = rst_group.gas_import_rate;
            if (crate != 0 || irate != 0) {
                // UDAs stored in output unit by convention
                const auto dim = this->m_static.m_unit_system.getDimension(UnitSystem::measure::gas_surface_rate);
                new_gconsump.add(rst_group.name, UDAValue(crate, dim), UDAValue(irate, dim), "", udq_undefined, this->m_static.m_unit_system);
            }
        }
        this->snapshots.back().gconsump.update( std::move(new_gconsump) );

        auto Glo = this->snapshots.back().glo();
        Glo.all_newton(rst_state.header.glift_all_nupcol);
        Glo.min_wait(rst_state.header.glift_min_wait);
        Glo.min_eco_gradient(rst_state.header.glift_min_eco_grad);
        Glo.gaslift_increment(rst_state.header.glift_rate_delta);

        for (std::size_t group_index = 0; group_index < rst_state.groups.size(); group_index++) {
            const auto& rst_group = rst_state.groups[group_index];

            if (rst_group.parent_group == 0)
                continue;

            if (rst_group.parent_group == rst_state.header.max_groups_in_field)
                continue;

            const auto& parent_group = rst_state.groups[rst_group.parent_group - 1];
            this->addGroupToGroup(parent_group.name, rst_group.name);

            if (GasLiftGroup::active(rst_group))
                Glo.add_group(GasLiftGroup(rst_group));
        }

        const auto& vfpprod = this->snapshots.back().vfpprod;
        for (const auto& rst_well : rst_state.wells) {
            std::optional<VFPProdTable::ALQ_TYPE> alq_type = vfpprod.has(rst_well.vfp_table)
                                                             ? std::optional(vfpprod.get(rst_well.vfp_table).getALQType())
                                                             : std::nullopt;
            auto well = Well {
                rst_well,
                report_step,
                rst_state.header.histctl_override,
                tracer_config,
                this->m_static.m_unit_system,
                rst_state.header.udq_undefined,
                alq_type
            };

            auto rst_connections = std::vector<Connection> {};
            std::transform(rst_well.connections.begin(), rst_well.connections.end(),
                           std::back_inserter(rst_connections),
                           [&grid, &fp](const auto& rst_conn)
                           {
                               return Connection{rst_conn, grid, fp};
                           });

            if (rst_well.segments.empty()) {
                auto connections = std::make_shared<WellConnections>
                    (order_from_int(rst_well.completion_ordering),
                     rst_well.ij[0], rst_well.ij[1], rst_connections);

                well.updateConnections(std::move(connections), grid);
            }
            else {
                auto rst_segments = std::unordered_map<int, Segment>{};
                for (const auto& rst_segment : rst_well.segments) {
                    rst_segments.try_emplace(rst_segment.segment, rst_segment, rst_well.name);
                }

                const auto& [connections, segments] =
                    Compsegs::rstUpdate(rst_well, rst_connections, rst_segments);

                well.updateConnections(std::make_shared<WellConnections>(connections), grid);
                well.updateSegments(std::make_shared<WellSegments>(segments));
            }

            this->addWell(well);
            this->addWellToGroup(well.groupName(), well.name(), report_step);

            OpmLog::info(fmt::format("Adding well {} from restart file", rst_well.name));

            if (GasLiftWell::active(rst_well)) {
                Glo.add_well(GasLiftWell(rst_well));
            }
        }

        this->snapshots.back().glo.update( std::move(Glo) );
        this->snapshots.back().update_tuning(rst_state.tuning);
        this->snapshots.back().events().addEvent( ScheduleEvents::TUNING_CHANGE );

        this->snapshots.back().update_oilvap(rst_state.oilvap);

        {
            const auto& header = rst_state.header;
            // A NONE target written to .UNRST may indicate no GUIDERAT (i.e., during history)
            const auto target = GuideRateModel::TargetFromRestart(header.guide_rate_nominated_phase);
            if ((target != GuideRateModel::Target::NONE) &&
                GuideRateModel::rst_valid(header.guide_rate_delay,
                                          header.guide_rate_a,
                                          header.guide_rate_b,
                                          header.guide_rate_c,
                                          header.guide_rate_d,
                                          header.guide_rate_e,
                                          header.guide_rate_f,
                                          header.guide_rate_damping))
            {
                const bool allow_increase = true;
                const bool use_free_gas = false;

                const auto guide_rate_model = GuideRateModel {
                    header.guide_rate_delay,
                    target,
                    header.guide_rate_a,
                    header.guide_rate_b,
                    header.guide_rate_c,
                    header.guide_rate_d,
                    header.guide_rate_e,
                    header.guide_rate_f,
                    allow_increase,
                    header.guide_rate_damping,
                    use_free_gas
                };

                this->updateGuideRateModel(guide_rate_model, report_step);
            }
        }

        if (rst_state.header.histctl_override > 0)
            this->snapshots.back().update_whistctl(WellProducerCModeFromInt(rst_state.header.histctl_override));

        for (const auto& rst_group : rst_state.groups) {
            auto& group = this->snapshots.back().groups.get( rst_group.name );
            if (group.isProductionGroup()) {
                auto new_config = this->snapshots.back().guide_rate();
                new_config.update_production_group(group);
                this->snapshots.back().guide_rate.update(std::move(new_config));
            }
            if (group.isInjectionGroup()) {
                // Set name of VREP group if different than default
                if (static_cast<int>(group.insert_index()) != rst_group.voidage_group_index) {
                    for (const auto& [phase, orig_inj_prop] : group.injectionProperties()) {
                        Group::GroupInjectionProperties inj_prop(orig_inj_prop);
                        inj_prop.voidage_group = rst_group_names[rst_group.voidage_group_index];
                        group.updateInjection(inj_prop);
                    }
                }
             }
        }

        this->snapshots.back().udq.update( UDQConfig(this->m_static.m_runspec.udqParams(), rst_state) );
        const auto& uda_records = UDQActive::load_rst( this->m_static.m_unit_system, this->snapshots.back().udq(), rst_state, this->wellNames(report_step), this->groupNames(report_step));
        if (!uda_records.empty()) {
            const auto& udq_config = this->snapshots.back().udq();
            auto udq_active = this->snapshots.back().udq_active();

            for (const auto& [control, value, wgname, ig_phase] : uda_records) {
                if (UDQ::well_control(control)) {
                    auto& well = this->snapshots.back().wells.get(wgname);

                    if (UDQ::is_well_injection_control(control, well.isInjector())) {
                        auto injection_properties = std::make_shared<Well::WellInjectionProperties>(well.getInjectionProperties());
                        injection_properties->update_uda(udq_config, udq_active, control, value);
                        well.updateInjection(std::move(injection_properties));
                    }

                    if (UDQ::is_well_production_control(control, well.isProducer())) {
                        auto production_properties = std::make_shared<Well::WellProductionProperties>(well.getProductionProperties());
                        production_properties->update_uda(udq_config, udq_active, control, value);
                        well.updateProduction(std::move(production_properties));
                    }
                } else {
                    auto& group = this->snapshots.back().groups.get(wgname);
                    if (UDQ::is_group_injection_control(control)) {
                        auto injection_properties = group.injectionProperties(ig_phase.value());
                        injection_properties.update_uda(udq_config, udq_active, control, value);
                        group.updateInjection(injection_properties);
                    }

                    if (UDQ::is_group_production_control(control)) {
                        auto production_properties = group.productionProperties();
                        production_properties.update_uda(udq_config, udq_active, control, value);
                        group.updateProduction(production_properties);
                    }
                }
            }
            this->snapshots.back().udq_active.update( std::move(udq_active) );
        }

        if (!rst_state.actions.empty()) {
            auto actions = this->snapshots.back().actions();
            for (const auto& rst_action : rst_state.actions)
                actions.add( Action::ActionX(rst_action) );
            this->snapshots.back().actions.update( std::move(actions) );
        }

        this->snapshots.back().wtest_config.update( WellTestConfig{rst_state, report_step});

        this->snapshots.back().network_balance.update(Network::Balance { rst_state.netbalan });

        for (const auto& aquflux : rst_state.aquifers.constantFlux()) {
            auto aqPos =  this->snapshots.back().aqufluxs
                .insert_or_assign(aquflux.aquiferID,
                                  SingleAquiferFlux { aquflux.aquiferID });

            aqPos.first->second.flux = aquflux.flow_rate;
            aqPos.first->second.active = true;
        }

        if (!rst_state.wlists.empty()) {
            this->snapshots.back().wlist_manager.update(WListManager(rst_state));
        }

        if (rst_state.network.isActive()) {
            auto network = this->snapshots.back().network();

            // Note: We presently support only the default value of BRANPROP(4).
            const auto alq_value =
                ParserKeywords::BRANPROP::ALQ::defaultValue;

            const auto& rst_nodes = rst_state.network.nodes();
            for (const auto& rst_branch : rst_state.network.branches()) {
                if ((rst_branch.down < 0) || (rst_branch.up < 0)) {
                    // Prune branches to non-existent nodes.
                    continue;
                }

                const auto& downtree_node = rst_nodes[rst_branch.down].name;
                const auto& uptree_node = rst_nodes[rst_branch.up].name;

                network.add_branch({ downtree_node, uptree_node, rst_branch.vfp, alq_value });
            }

            for (const auto& rst_node : rst_nodes) {
                auto node = Network::Node { rst_node.name };

                if (rst_node.terminal_pressure.has_value()) {
                    node.terminal_pressure(rst_node.terminal_pressure.value());
                }

                if (rst_node.as_choke.has_value()) {
                    node.as_choke(rst_node.as_choke.value());
                }

                node.add_gas_lift_gas(rst_node.add_lift_gas);

                network.update_node(std::move(node));
            }

            this->snapshots.back().network.update(std::move(network));
        }
    }

    std::shared_ptr<const Python> Schedule::python() const
    {
        return this->m_static.m_python_handle;
    }


    const GasLiftOpt& Schedule::glo(std::size_t report_step) const {
        return this->snapshots[report_step].glo();
    }

namespace {
/*
  The insane trickery here (thank you Stackoverflow!) is to be able to provide a
  simple templated comparison function

     template <typename T>
     int not_equal(const T& arg1, const T& arg2, const std::string& msg);

  which will print arg1 and arg2 on stderr *if* T supports operator<<, otherwise
  it will just print the typename of T.
*/


template<typename T, typename = int>
struct cmpx
{
    int neq(const T& arg1, const T& arg2, const std::string& msg) {
        if (arg1 == arg2)
            return 0;

        std::cerr << "Error when comparing <" << typeid(arg1).name() << ">: " << msg << std::endl;
        return 1;
    }
};

template <typename T>
struct cmpx<T, decltype(std::cout << T(), 0)>
{
    int neq(const T& arg1, const T& arg2, const std::string& msg) {
        if (arg1 == arg2)
            return 0;

        std::cerr << "Error when comparing: " << msg << " " << arg1 << " != " << arg2 << std::endl;
        return 1;
    }
};


template <typename T>
int not_equal(const T& arg1, const T& arg2, const std::string& msg) {
    return cmpx<T>().neq(arg1, arg2, msg);
}


template <>
int not_equal(const double& arg1, const double& arg2, const std::string& msg) {
    if (Opm::cmp::scalar_equal(arg1, arg2))
        return 0;

    std::cerr << "Error when comparing: " << msg << " " << arg1 << " != " << arg2 << std::endl;
    return 1;
}

template <>
int not_equal(const UDAValue& arg1, const UDAValue& arg2, const std::string& msg) {
    if (arg1.is<double>())
        return not_equal( arg1.get<double>(), arg2.get<double>(), msg);
    else
        return not_equal( arg1.get<std::string>(), arg2.get<std::string>(), msg);
}


std::string well_msg(const std::string& well, const std::string& msg) {
    return "Well: " + well + " " + msg;
}

std::string well_segment_msg(const std::string& well, int segment_number, const std::string& msg) {
    return "Well: " + well + " Segment: " + std::to_string(segment_number) + " " + msg;
}

std::string well_connection_msg(const std::string& well, const Connection& conn, const std::string& msg) {
    return "Well: " + well + " Connection: " + std::to_string(conn.getI()) + ", " + std::to_string(conn.getJ()) + ", " + std::to_string(conn.getK()) + "  " + msg;
}

}

bool Schedule::cmp(const Schedule& sched1, const Schedule& sched2, std::size_t report_step) {
    int count = not_equal(sched1.wellNames(report_step), sched2.wellNames(report_step), "Wellnames");
    if (count != 0)
        return false;
    {
        //if (sched1.size() != sched2.size())
        //    return false;

        //for (std::size_t step=0; step < sched1.size(); step++) {
        //    auto start1 = sched1[step].start_time();
        //    auto start2 = sched2[step].start_time();
        //    if (start1 != start2)
        //        return false;

        //    if (step < sched1.size() - 1) {
        //        auto end1 = sched1[step].end_time();
        //        auto end2 = sched2[step].end_time();
        //        if (end1 != end2)
        //            return false;
        //    }
        //}
    }

    for (const auto& wname : sched1.wellNames(report_step)) {
        const auto& well1 = sched1.getWell(wname, report_step);
        const auto& well2 = sched2.getWell(wname, report_step);
        int well_count = 0;
        {
            const auto& connections2 = well2.getConnections();
            const auto& connections1 = well1.getConnections();

            well_count += not_equal( connections1.ordering(), connections2.ordering(), well_msg(well1.name(), "Connection: ordering"));
            for (std::size_t icon = 0; icon < connections1.size(); icon++) {
                const auto& conn1 = connections1[icon];
                const auto& conn2 = connections2[icon];
                well_count += not_equal( conn1.getI(), conn2.getI(), well_connection_msg(well1.name(), conn1, "I"));
                well_count += not_equal( conn1.getJ() , conn2.getJ() , well_connection_msg(well1.name(), conn1, "J"));
                well_count += not_equal( conn1.getK() , conn2.getK() , well_connection_msg(well1.name(), conn1, "K"));
                well_count += not_equal( conn1.state() , conn2.state(), well_connection_msg(well1.name(), conn1, "State"));
                well_count += not_equal( conn1.dir() , conn2.dir(), well_connection_msg(well1.name(), conn1, "dir"));
                well_count += not_equal( conn1.complnum() , conn2.complnum(), well_connection_msg(well1.name(), conn1, "complnum"));
                well_count += not_equal( conn1.segment() , conn2.segment(), well_connection_msg(well1.name(), conn1, "segment"));
                well_count += not_equal( conn1.kind() , conn2.kind(), well_connection_msg(well1.name(), conn1, "CFKind"));
                well_count += not_equal( conn1.sort_value(), conn2.sort_value(), well_connection_msg(well1.name(), conn1, "sort_value"));


                well_count += not_equal( conn1.CF(), conn2.CF(), well_connection_msg(well1.name(), conn1, "CF"));
                well_count += not_equal( conn1.Kh(), conn2.Kh(), well_connection_msg(well1.name(), conn1, "Kh"));
                well_count += not_equal( conn1.rw(), conn2.rw(), well_connection_msg(well1.name(), conn1, "rw"));
                well_count += not_equal( conn1.depth(), conn2.depth(), well_connection_msg(well1.name(), conn1, "depth"));

                //well_count += not_equal( conn1.r0(), conn2.r0(), well_connection_msg(well1.name(), conn1, "r0"));
                well_count += not_equal( conn1.skinFactor(), conn2.skinFactor(), well_connection_msg(well1.name(), conn1, "skinFactor"));

            }
        }

        if (not_equal(well1.isMultiSegment(), well2.isMultiSegment(), well_msg(well1.name(), "Is MSW")))
            return false;

        if (well1.isMultiSegment()) {
            const auto& segments1 = well1.getSegments();
            const auto& segments2 = well2.getSegments();
            if (not_equal(segments1.size(), segments2.size(), "Segments: size"))
                return false;

            for (std::size_t iseg=0; iseg < segments1.size(); iseg++) {
                const auto& segment1 = segments1[iseg];
                const auto& segment2 = segments2[iseg];
                //const auto& segment2 = segments2.getFromSegmentNumber(segment1.segmentNumber());
                well_count += not_equal(segment1.segmentNumber(), segment2.segmentNumber(), well_segment_msg(well1.name(), segment1.segmentNumber(), "segmentNumber"));
                well_count += not_equal(segment1.branchNumber(), segment2.branchNumber(), well_segment_msg(well1.name(), segment1.segmentNumber(), "branchNumber"));
                well_count += not_equal(segment1.outletSegment(), segment2.outletSegment(), well_segment_msg(well1.name(), segment1.segmentNumber(), "outletSegment"));
                well_count += not_equal(segment1.totalLength(), segment2.totalLength(), well_segment_msg(well1.name(), segment1.segmentNumber(), "totalLength"));
                well_count += not_equal(segment1.depth(), segment2.depth(), well_segment_msg(well1.name(), segment1.segmentNumber(), "depth"));
                well_count += not_equal(segment1.internalDiameter(), segment2.internalDiameter(), well_segment_msg(well1.name(), segment1.segmentNumber(), "internalDiameter"));
                well_count += not_equal(segment1.roughness(), segment2.roughness(), well_segment_msg(well1.name(), segment1.segmentNumber(), "roughness"));
                well_count += not_equal(segment1.crossArea(), segment2.crossArea(), well_segment_msg(well1.name(), segment1.segmentNumber(), "crossArea"));
                well_count += not_equal(segment1.volume(), segment2.volume(), well_segment_msg(well1.name(), segment1.segmentNumber(), "volume"));
            }
        }

        well_count += not_equal(well1.getStatus(), well2.getStatus(), well_msg(well1.name(), "status"));
        {
            const auto& prod1 = well1.getProductionProperties();
            const auto& prod2 = well2.getProductionProperties();
            well_count += not_equal(prod1.name, prod2.name , well_msg(well1.name(), "Prod: name"));
            well_count += not_equal(prod1.OilRate, prod2.OilRate, well_msg(well1.name(), "Prod: OilRate"));
            well_count += not_equal(prod1.GasRate, prod2.GasRate, well_msg(well1.name(), "Prod: GasRate"));
            well_count += not_equal(prod1.WaterRate, prod2.WaterRate, well_msg(well1.name(), "Prod: WaterRate"));
            well_count += not_equal(prod1.LiquidRate, prod2.LiquidRate, well_msg(well1.name(), "Prod: LiquidRate"));
            well_count += not_equal(prod1.ResVRate, prod2.ResVRate, well_msg(well1.name(), "Prod: ResVRate"));
            well_count += not_equal(prod1.BHPTarget, prod2.BHPTarget, well_msg(well1.name(), "Prod: BHPTarget"));
            well_count += not_equal(prod1.THPTarget, prod2.THPTarget, well_msg(well1.name(), "Prod: THPTarget"));
            well_count += not_equal(prod1.VFPTableNumber, prod2.VFPTableNumber, well_msg(well1.name(), "Prod: VFPTableNumber"));
            well_count += not_equal(prod1.ALQValue, prod2.ALQValue, well_msg(well1.name(), "Prod: ALQValue"));
            well_count += not_equal(prod1.predictionMode, prod2.predictionMode, well_msg(well1.name(), "Prod: predictionMode"));
            if (!prod1.predictionMode) {
                well_count += not_equal(prod1.bhp_hist_limit, prod2.bhp_hist_limit, well_msg(well1.name(), "Prod: bhp_hist_limit"));
                well_count += not_equal(prod1.thp_hist_limit, prod2.thp_hist_limit, well_msg(well1.name(), "Prod: thp_hist_limit"));
                well_count += not_equal(prod1.BHPH, prod2.BHPH, well_msg(well1.name(), "Prod: BHPH"));
                well_count += not_equal(prod1.THPH, prod2.THPH, well_msg(well1.name(), "Prod: THPH"));
            }
            well_count += not_equal(prod1.productionControls(), prod2.productionControls(), well_msg(well1.name(), "Prod: productionControls"));
            if (well1.getStatus() == Well::Status::OPEN)
                well_count += not_equal(prod1.controlMode, prod2.controlMode, well_msg(well1.name(), "Prod: controlMode"));
            well_count += not_equal(prod1.whistctl_cmode, prod2.whistctl_cmode, well_msg(well1.name(), "Prod: whistctl_cmode"));
        }
        {
            const auto& inj1 = well1.getInjectionProperties();
            const auto& inj2 = well2.getInjectionProperties();

            well_count += not_equal(inj1.name, inj2.name, well_msg(well1.name(), "Well::Inj: name"));
            well_count += not_equal(inj1.surfaceInjectionRate, inj2.surfaceInjectionRate, well_msg(well1.name(), "Well::Inj: surfaceInjectionRate"));
            well_count += not_equal(inj1.reservoirInjectionRate, inj2.reservoirInjectionRate, well_msg(well1.name(), "Well::Inj: reservoirInjectionRate"));
            well_count += not_equal(inj1.BHPTarget, inj2.BHPTarget, well_msg(well1.name(), "Well::Inj: BHPTarget"));
            well_count += not_equal(inj1.THPTarget, inj2.THPTarget, well_msg(well1.name(), "Well::Inj: THPTarget"));
            well_count += not_equal(inj1.bhp_hist_limit, inj2.bhp_hist_limit, well_msg(well1.name(), "Well::Inj: bhp_hist_limit"));
            well_count += not_equal(inj1.thp_hist_limit, inj2.thp_hist_limit, well_msg(well1.name(), "Well::Inj: thp_hist_limit"));
            well_count += not_equal(inj1.BHPH, inj2.BHPH, well_msg(well1.name(), "Well::Inj: BHPH"));
            well_count += not_equal(inj1.THPH, inj2.THPH, well_msg(well1.name(), "Well::Inj: THPH"));
            well_count += not_equal(inj1.VFPTableNumber, inj2.VFPTableNumber, well_msg(well1.name(), "Well::Inj: VFPTableNumber"));
            well_count += not_equal(inj1.predictionMode, inj2.predictionMode, well_msg(well1.name(), "Well::Inj: predictionMode"));
            well_count += not_equal(inj1.injectionControls, inj2.injectionControls, well_msg(well1.name(), "Well::Inj: injectionControls"));
            well_count += not_equal(inj1.injectorType, inj2.injectorType, well_msg(well1.name(), "Well::Inj: injectorType"));
            well_count += not_equal(inj1.controlMode, inj2.controlMode, well_msg(well1.name(), "Well::Inj: controlMode"));
        }

        {
            well_count += well2.firstTimeStep() > report_step;
            well_count += not_equal( well1.groupName(), well2.groupName(), well_msg(well1.name(), "Well: groupName"));
            well_count += not_equal( well1.getHeadI(), well2.getHeadI(), well_msg(well1.name(), "Well: getHeadI"));
            well_count += not_equal( well1.getHeadJ(), well2.getHeadJ(), well_msg(well1.name(), "Well: getHeadJ"));
            well_count += not_equal( well1.getRefDepth(), well2.getRefDepth(), well_msg(well1.name(), "Well: getRefDepth"));
            well_count += not_equal( well1.isMultiSegment(), well2.isMultiSegment() , well_msg(well1.name(), "Well: isMultiSegment"));
            well_count += not_equal( well1.isAvailableForGroupControl(), well2.isAvailableForGroupControl() , well_msg(well1.name(), "Well: isAvailableForGroupControl"));
            well_count += not_equal( well1.getGuideRate(), well2.getGuideRate(), well_msg(well1.name(), "Well: getGuideRate"));
            well_count += not_equal( well1.getGuideRatePhase(), well2.getGuideRatePhase(), well_msg(well1.name(), "Well: getGuideRatePhase"));
            well_count += not_equal( well1.getGuideRateScalingFactor(), well2.getGuideRateScalingFactor(), well_msg(well1.name(), "Well: getGuideRateScalingFactor"));
            well_count += not_equal( well1.predictionMode(), well2.predictionMode(), well_msg(well1.name(), "Well: predictionMode"));
            well_count += not_equal( well1.isProducer(), well2.isProducer(), well_msg(well1.name(), "Well: isProducer"));
            well_count += not_equal( well1.isInjector(), well2.isInjector(), well_msg(well1.name(), "Well: isInjector"));
            if (well1.isInjector())
                well_count += not_equal( well1.injectorType(), well2.injectorType(), well_msg(well1.name(), "Well1: injectorType"));
            well_count += not_equal( well1.seqIndex(), well2.seqIndex(), well_msg(well1.name(), "Well: seqIndex"));
            well_count += not_equal( well1.getAutomaticShutIn(), well2.getAutomaticShutIn(), well_msg(well1.name(), "Well: getAutomaticShutIn"));
            well_count += not_equal( well1.getAllowCrossFlow(), well2.getAllowCrossFlow(), well_msg(well1.name(), "Well: getAllowCrossFlow"));
            well_count += not_equal( well1.getSolventFraction(), well2.getSolventFraction(), well_msg(well1.name(), "Well: getSolventFraction"));
            well_count += not_equal( well1.getStatus(), well2.getStatus(), well_msg(well1.name(), "Well: getStatus"));
            //well_count += not_equal( well1.getInjectionProperties(), well2.getInjectionProperties(), "Well: getInjectionProperties");


            if (well1.isProducer())
                well_count += not_equal( well1.getPreferredPhase(), well2.getPreferredPhase(), well_msg(well1.name(), "Well: getPreferredPhase"));
            well_count += not_equal( well1.getDrainageRadius(), well2.getDrainageRadius(), well_msg(well1.name(), "Well: getDrainageRadius"));
            well_count += not_equal( well1.getEfficiencyFactor(), well2.getEfficiencyFactor(), well_msg(well1.name(), "Well: getEfficiencyFactor"));
        }
        count += well_count;
        if (well_count > 0)
            std::cerr << std::endl;
    }
    return (count == 0);
}

const ScheduleState& Schedule::back() const {
    return this->snapshots.back();
}

const ScheduleState& Schedule::operator[](std::size_t index) const {
    return this->snapshots.at(index);
}

std::vector<ScheduleState>::const_iterator Schedule::begin() const {
    return this->snapshots.begin();
}

std::vector<ScheduleState>::const_iterator Schedule::end() const {
    return this->snapshots.end();
}

void Schedule::create_first(const time_point& start_time, const std::optional<time_point>& end_time) {
    if (end_time.has_value())
        this->snapshots.emplace_back( start_time, end_time.value() );
    else
        this->snapshots.emplace_back(start_time);

    const auto& run_spec = this->m_static.m_runspec;
    auto& sched_state = snapshots.back();
    sched_state.init_nupcol(run_spec.nupcol());
    if (this->m_static.oilVap.has_value()) {
        sched_state.update_oilvap(*this->m_static.oilVap);
    } else {
        sched_state.update_oilvap( OilVaporizationProperties(run_spec.tabdims().getNumPVTTables()));
    }
    sched_state.update_message_limits( this->m_static.m_deck_message_limits );
    sched_state.pavg.update( PAvg() );
    sched_state.wtest_config.update( WellTestConfig() );
    sched_state.gconsale.update( GConSale() );
    sched_state.gconsump.update( GConSump() );
    sched_state.gsatprod.update( GSatProd() );
    sched_state.gecon.update( GroupEconProductionLimits() );
    sched_state.wlist_manager.update( WListManager() );
    sched_state.network.update( Network::ExtNetwork() );
    sched_state.rescoup.update( ReservoirCoupling::CouplingInfo() );
    sched_state.rpt_config.update( RPTConfig() );
    sched_state.actions.update( Action::Actions() );
    sched_state.udq_active.update( UDQActive() );
    sched_state.well_order.update( NameOrder() );
    sched_state.group_order.update( GroupOrder(run_spec.wellDimensions().maxGroupsInField()));
    sched_state.udq.update( UDQConfig(run_spec.udqParams()));
    sched_state.glo.update( GasLiftOpt() );
    sched_state.guide_rate.update( GuideRateConfig() );
    sched_state.rft_config.update( RFTConfig() );
    sched_state.rst_config.update( RSTConfig::first( this->m_static.rst_config ) );
    sched_state.network_balance.update(Network::Balance{run_spec.networkDimensions().active()});
    sched_state.update_sumthin(this->m_static.sumthin);
    sched_state.rptonly(this->m_static.rptonly);
    sched_state.bhp_defaults.update( ScheduleState::BHPDefaults() );
    sched_state.source.update( Source() );
    sched_state.wcycle.update( WCYCLE() );
    //sched_state.update_date( start_time );
    this->addGroup("FIELD", 0);
}

void Schedule::create_next(const time_point& start_time, const std::optional<time_point>& end_time) {
    if (this->snapshots.empty())
        this->create_first(start_time, end_time);
    else {
        const auto& last = this->snapshots.back();
        if (end_time.has_value())
            this->snapshots.emplace_back( last, start_time, end_time.value() );
        else
            this->snapshots.emplace_back( last, start_time );
    }
}


void Schedule::create_next(const ScheduleBlock& block) {
    const auto& start_time = block.start_time();
    const auto& end_time = block.end_time();
    this->create_next(start_time, end_time);
}

void Schedule::dump_deck(std::ostream& os) const
{
    this->m_sched_deck.dump_deck(os, this->getUnits());
}

std::ostream& operator<<(std::ostream& os, const Schedule& sched)
{
    sched.dump_deck(os);
    return os;
}

}
