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
#ifndef SCHEDULE_HPP
#define SCHEDULE_HPP

#include <cstddef>
#include <ctime>
#include <functional>
#include <iosfwd>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <opm/input/eclipse/Schedule/Action/ActionResult.hpp>
#include <opm/input/eclipse/Schedule/Action/SimulatorUpdate.hpp>
#include <opm/input/eclipse/Schedule/Action/WGNames.hpp>
#include <opm/input/eclipse/Schedule/CompletedCells.hpp>
#include <opm/input/eclipse/Schedule/Group/Group.hpp>
#include <opm/input/eclipse/Schedule/ScheduleDeck.hpp>
#include <opm/input/eclipse/Schedule/ScheduleState.hpp>
#include <opm/input/eclipse/Schedule/ScheduleStatic.hpp>
#include <opm/input/eclipse/Schedule/Well/PAvg.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/Connection.hpp>
#include <opm/input/eclipse/Schedule/WriteRestartFileEvents.hpp>
#include <opm/input/eclipse/Units/UnitSystem.hpp>

namespace Opm {
    class ActiveGridCells;
    class Deck;
    class DeckKeyword;
    class DeckRecord;
    enum class ConnectionOrder;
    class EclipseGrid;
    class EclipseState;
    class ErrorGuard;
    class FieldPropsManager;
    class GasLiftOpt;
    class GTNode;
    class GuideRateConfig;
    class GuideRateModel;
    class HandlerContext;
    enum class InputErrorAction;
    class ParseContext;
    class Python;
    namespace ReservoirCoupling {
        class CouplingInfo;
    }
    class Runspec;
    class RPTConfig;
    class ScheduleGrid;
    class SCHEDULESection;
    class SegmentMatcher;
    class SummaryState;
    class TracerConfig;
    class UDQConfig;
    class Well;
    enum class WellGasInflowEquation;
    class WellMatcher;
    enum class WellProducerCMode;
    enum class WellStatus;
    class WelSegsSet;
    class WellTestConfig;
} // namespace Opm

namespace Opm::Action {
    class ActionX;
    class PyAction;
    class State;
} // namespace Opm::Action

namespace Opm::RestartIO {
    struct RstState;
} // namespace Opm::RestartIO

namespace Opm {
    class Schedule
    {
    public:
        Schedule() = default;

        explicit Schedule(std::shared_ptr<const Python> python_handle);

        /*! \brief Construct a Schedule object from a deck.
         *  \param deck Deck to construct Schedule from
         *  \param fp Field property manager
         *  \param runspec Run specification parameters to use
         *  \param parseContext Parsing context
         *  \param errors Error configuration
         *  \param python Python interpreter to use
         *  \param lowActionParsingStrictness Reduce parsing strictness for actions
         *  \param slave_mode Slave mode flag
         *  \param keepKeywords Keep the schdule keywords even if there are no actions
         *  \param output_interval Output interval to use
         *  \param rst Restart state to use
         *  \param tracer_config Tracer configuration to use
         */
        Schedule(const Deck& deck,
                 const EclipseGrid& grid,
                 const FieldPropsManager& fp,
                 const Runspec &runspec,
                 const ParseContext& parseContext,
                 ErrorGuard& errors,
                 std::shared_ptr<const Python> python,
                 const bool lowActionParsingStrictness = false,
                 const bool slave_mode = false,
                 const bool keepKeywords = true,
                 const std::optional<int>& output_interval = {},
                 const RestartIO::RstState* rst = nullptr,
                 const TracerConfig* tracer_config = nullptr);

        template<typename T>
        Schedule(const Deck& deck,
                 const EclipseGrid& grid,
                 const FieldPropsManager& fp,
                 const Runspec &runspec,
                 const ParseContext& parseContext,
                 T&& errors,
                 std::shared_ptr<const Python> python,
                 const bool lowActionParsingStrictness = false,
                 const bool slave_mode = false,
                 const bool keepKeywords = true,
                 const std::optional<int>& output_interval = {},
                 const RestartIO::RstState* rst = nullptr,
                 const TracerConfig* tracer_config = nullptr);

        Schedule(const Deck& deck,
                 const EclipseGrid& grid,
                 const FieldPropsManager& fp,
                 const Runspec &runspec,
                 std::shared_ptr<const Python> python,
                 const bool lowActionParsingStrictness = false,
                 const bool slave_mode = false,
                 const bool keepKeywords = true,
                 const std::optional<int>& output_interval = {},
                 const RestartIO::RstState* rst = nullptr,
                 const TracerConfig* tracer_config = nullptr);

        Schedule(const Deck& deck,
                 const EclipseState& es,
                 const ParseContext& parseContext,
                 ErrorGuard& errors,
                 std::shared_ptr<const Python> python,
                 const bool lowActionParsingStrictness = false,
                 const bool slave_mode = false,
                 const bool keepKeywords = true,
                 const std::optional<int>& output_interval = {},
                 const RestartIO::RstState* rst = nullptr);

        template <typename T>
        Schedule(const Deck& deck,
                 const EclipseState& es,
                 const ParseContext& parseContext,
                 T&& errors,
                 std::shared_ptr<const Python> python,
                 const bool lowActionParsingStrictness = false,
                 const bool slave_mode = false,
                 const bool keepKeywords = true,
                 const std::optional<int>& output_interval = {},
                 const RestartIO::RstState* rst = nullptr);

        Schedule(const Deck& deck,
                 const EclipseState& es,
                 std::shared_ptr<const Python> python,
                 const bool lowActionParsingStrictness = false,
                 const bool slave_mode = false,
                 const bool keepKeywords = true,
                 const std::optional<int>& output_interval = {},
                 const RestartIO::RstState* rst = nullptr);

        // The constructor *without* the Python arg should really only be used from Python itself
        Schedule(const Deck& deck,
                 const EclipseState& es,
                 const std::optional<int>& output_interval = {},
                 const RestartIO::RstState* rst = nullptr);

        ~Schedule() = default;

        static Schedule serializationTestObject();

        /*
         * If the input deck does not specify a start time, Eclipse's 1. Jan
         * 1983 is defaulted
         */
        std::time_t getStartTime() const;
        std::time_t posixStartTime() const;
        std::time_t posixEndTime() const;
        std::time_t simTime(std::size_t timeStep) const;
        double seconds(std::size_t timeStep) const;
        double stepLength(std::size_t timeStep) const;
        std::optional<int> exitStatus() const;
        const UnitSystem& getUnits() const { return this->m_static.m_unit_system; }
        const Runspec& runspec() const { return this->m_static.m_runspec; }

        std::size_t numWells() const;
        std::size_t numWells(std::size_t timestep) const;
        bool hasWell(const std::string& wellName) const;
        bool hasWell(const std::string& wellName, std::size_t timeStep) const;

        WellMatcher wellMatcher(std::size_t report_step) const;
        std::function<std::unique_ptr<SegmentMatcher>()> segmentMatcherFactory(std::size_t report_step) const;
        std::vector<std::string> wellNames(const std::string& pattern, std::size_t timeStep, const std::vector<std::string>& matching_wells = {}) const;
        std::vector<std::string> wellNames(const std::string& pattern) const;
        std::vector<std::string> wellNames(std::size_t timeStep) const;
        std::vector<std::string> wellNames() const;
        //const ScheduleGrid& get_schedulegrid() const{ return this->grid; };
        /// Query for group existence at particular time
        ///
        /// \param[in] groupName Fully specified group name.
        ///
        /// \param[in] timeStep Zero-based report step index.
        ///
        /// \return Whether or not group \p groupName exists at time \p
        /// timeStep.
        bool hasGroup(const std::string& groupName, std::size_t timeStep) const;

        /// Retrieve names of all groups at particular time whose names
        /// match a pattern.
        ///
        /// \param[in] pattern Group name pattern.  Expected to be a group
        /// name root such as 'PR*' or a fully specified group name such as
        /// 'PROD'.
        ///
        /// \param[in] timeStep Zero-based report step index.
        ///
        /// \return Names of all groups at \p timeStep whose names match the
        /// \p pattern.
        std::vector<std::string>
        groupNames(const std::string& pattern, std::size_t timeStep) const;

        /// Retrieve names of all groups at a particular time.
        ///
        /// \param[in] timeStep Zero-based report step index.
        ///
        /// \return Names of all groups in simulation which are active at
        /// time \p timeStep.
        const std::vector<std::string>& groupNames(std::size_t timeStep) const;

        /// Retrieve names of all groups matching a pattern
        ///
        /// \param[in] pattern Group name pattern.  Expected to be a group
        /// name root such as 'PR*' or a fully specified group name such as
        /// 'PROD'.
        ///
        /// \return Names of all groups in simulation run matching the \p
        /// pattern.
        std::vector<std::string> groupNames(const std::string& pattern) const;

        /// Retrieve names of all groups in model
        ///
        /// Includes FIELD group.
        const std::vector<std::string>& groupNames() const;

        /// Retrieve collection of group objects suiteable for restart file output
        ///
        /// The collection is sized and ordered as follows
        ///
        ///  -# The number of elements is WELLDIMS::MAXGROUPS + 1
        ///  -# The elements are sorted according to group.insert_index().
        ///  -# Nullptr represents an unused slot.
        ///  -# The very last element corresponds to the FIELD group.
        ///
        /// \param[in] timeStep Zero-based report step index.
        std::vector<const Group*> restart_groups(std::size_t timeStep) const;

        std::vector<std::string> changed_wells(std::size_t reportStep) const;
        const Well& getWell(std::size_t well_index, std::size_t timeStep) const;
        const Well& getWell(const std::string& wellName, std::size_t timeStep) const;
        const Well& getWellatEnd(const std::string& well_name) const;
        // get the list of the constant flux aquifer specified in the whole schedule
        std::unordered_set<int> getAquiferFluxSchedule() const;
        std::vector<Well> getWells(std::size_t timeStep) const;
        std::vector<Well> getWellsatEnd() const;

        // Get wells that have been active any time during simulation
        std::vector<Well> getActiveWellsAtEnd() const;

        // Get well names of wells that have never been active
        std::vector<std::string> getInactiveWellNamesAtEnd() const;

        const std::unordered_map<std::string, std::set<int>>&
        getPossibleFutureConnections() const;

        void shut_well(const std::string& well_name, std::size_t report_step);
        void shut_well(const std::string& well_name);
        void stop_well(const std::string& well_name, std::size_t report_step);
        void stop_well(const std::string& well_name);
        void open_well(const std::string& well_name, std::size_t report_step);
        void open_well(const std::string& well_name);
        void clear_event(ScheduleEvents::Events, std::size_t report_step);
        void add_event(ScheduleEvents::Events, std::size_t report_step);
        void applyWellProdIndexScaling(const std::string& well_name, const std::size_t reportStep, const double scalingFactor);

        WellProducerCMode getGlobalWhistctlMmode(std::size_t timestep) const;

        const UDQConfig& getUDQConfig(std::size_t timeStep) const;
        void evalAction(const SummaryState& summary_state, std::size_t timeStep);

        GTNode groupTree(std::size_t report_step) const;
        GTNode groupTree(const std::string& root_node, std::size_t report_step) const;
        const Group& getGroup(const std::string& groupName, std::size_t timeStep) const;

        std::optional<std::size_t> first_RFT() const;
        /*
          Will remove all completions which are connected to cell which is not
          active. Will scan through all wells and all timesteps.
        */
        void filterConnections(const ActiveGridCells& grid);
        std::size_t size() const;

        bool write_rst_file(std::size_t report_step) const;
        const std::map< std::string, int >& rst_keywords( size_t timestep ) const;

        // The applyAction() member function is invoked from the simulator
        // *after* an ACTIONX has triggered.  Its return value is a small
        // structure with 'information' which the simulator should take into
        // account when updating internal datastructures after the ACTIONX
        // keywords have been applied.
        SimulatorUpdate applyAction(std::size_t reportStep,
                                    const Action::ActionX& action,
                                    const Action::Result::MatchingEntities& matches,
                                    const std::unordered_map<std::string, double>& wellpi);

        SimulatorUpdate applyAction(std::size_t reportStep,
                                    const Action::ActionX& action,
                                    const Action::Result::MatchingEntities& matches,
                                    const std::unordered_map<std::string, float>& wellpi);
        /*
          The runPyAction() will run the Python script in a PYACTION keyword. In
          the case of Schedule updates the recommended way of doing that from
          PYACTION is to invoke a "normal" ACTIONX keyword internally from the
          Python code. he return value from runPyAction() comes from such a
          internal ACTIONX.
        */
        SimulatorUpdate runPyAction(std::size_t reportStep,
                                    const Action::PyAction& pyaction,
                                    Action::State& action_state,
                                    EclipseState& ecl_state,
                                    SummaryState& summary_state);

        SimulatorUpdate runPyAction(std::size_t reportStep,
                                    const Action::PyAction& pyaction,
                                    Action::State& action_state,
                                    EclipseState& ecl_state,
                                    SummaryState& summary_state,
                                    const std::unordered_map<std::string, double>& target_wellpi);

        SimulatorUpdate runPyAction(std::size_t reportStep,
                                    const Action::PyAction& pyaction,
                                    Action::State& action_state,
                                    EclipseState& ecl_state,
                                    SummaryState& summary_state,
                                    const std::unordered_map<std::string, float>& target_wellpi);

        SimulatorUpdate modifyCompletions(const std::size_t reportStep,
                                          const std::map<std::string, std::vector<Connection>>& extraConns);

        const GasLiftOpt& glo(std::size_t report_step) const;

        bool operator==(const Schedule& data) const;
        std::shared_ptr<const Python> python() const;


        const ScheduleState& back() const;
        const ScheduleState& operator[](std::size_t index) const;
        std::vector<ScheduleState>::const_iterator begin() const;
        std::vector<ScheduleState>::const_iterator end() const;
        void create_next(const time_point& start_time, const std::optional<time_point>& end_time);
        void create_next(const ScheduleBlock& block);
        void create_first(const time_point& start_time, const std::optional<time_point>& end_time);

        void treat_critical_as_non_critical(bool value) { this->m_treat_critical_as_non_critical = value; }

        /*
          The cmp() function compares two schedule instances in a context aware
          manner. Floating point numbers are compared with a tolerance. The
          purpose of this comparison function is to implement regression tests
          for the schedule instances created by loading a restart file.
        */
        static bool cmp(const Schedule& sched1, const Schedule& sched2, std::size_t report_step);
        void applyKeywords(std::vector<std::unique_ptr<DeckKeyword>>& keywords, std::unordered_map<std::string, double>& target_wellpi, bool action_mode, std::size_t report_step);
        void applyKeywords(std::vector<std::unique_ptr<DeckKeyword>>& keywords, std::unordered_map<std::string, double>& target_wellpi, bool action_mode);

        template<class Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(this->m_static);
            serializer(this->m_sched_deck);
            serializer(this->action_wgnames);
            serializer(this->potential_wellopen_patterns);
            serializer(this->exit_status);
            serializer(this->snapshots);
            serializer(this->restart_output);
            serializer(this->completed_cells);
            serializer(this->completed_cells_lgr);
            serializer(this->m_treat_critical_as_non_critical);
            serializer(this->current_report_step);
            serializer(this->m_lowActionParsingStrictness);
            serializer(this->simUpdateFromPython);

            // If we are deserializing we need to setup the pointer to the
            // unit system since this is process specific. This is safe
            // because we set the same value in all well instances.
            // We do some redundant assignments as these are shared_ptr's
            // with multiple pointers to any given instance, but it is not
            // significant so let's keep it simple.
            if (!serializer.isSerializing()) {
                for (auto& snapshot : snapshots) {
                    for (auto& well : snapshot.wells) {
                        well.second->updateUnitSystem(&m_static.m_unit_system);
                    }
                }
            }
        }

        template <typename T>
        std::vector<std::pair<std::size_t,  T>> unique() const
        {
            std::vector<std::pair<std::size_t, T>> values;
            for (std::size_t index = 0; index < this->snapshots.size(); index++) {
                const auto& member = this->snapshots[index].get<T>();
                const auto& value = member.get();
                if (values.empty() || !(value == values.back().second))
                    values.push_back( std::make_pair(index, value));
            }
            return values;
        }

        friend std::ostream& operator<<(std::ostream& os, const Schedule& sched);
        void dump_deck(std::ostream& os) const;

    private:
        friend class HandlerContext;

        // Please update the member functions
        //   - operator==(const Schedule&) const
        //   - serializationTestObject()
        //   - serializeOp(Serializer&)
        // when you update/change this list of data members.
        bool m_treat_critical_as_non_critical = false;
        ScheduleStatic m_static{};
        ScheduleDeck m_sched_deck{};
        Action::WGNames action_wgnames{};
        std::unordered_set<std::string> potential_wellopen_patterns{}; // Set of well name patterns that potentially can open
        std::optional<int> exit_status{};
        std::vector<ScheduleState> snapshots{};
        WriteRestartFileEvents restart_output{};
        CompletedCells completed_cells{};
        std::vector<CompletedCells> completed_cells_lgr{};

        // Boolean indicating the strictness of parsing process for ActionX and PyAction.
        // If lowActionParsingStrictness is true, the simulator tries to apply unsupported
        // keywords, if lowActionParsingStrictness is false, the simulator only applies
        // supported keywords.
        bool m_lowActionParsingStrictness = false;

        // This unordered_map contains possible future connections of wells that might get added through an ACTIONX.
        // For parallel runs, this unordered_map is retrieved by the grid partitioner to ensure these connections
        // end up on the same partition.
        std::unordered_map<std::string, std::set<int>> possibleFutureConnections;

        // The current_report_step is set to the current report step when a PYACTION call is executed.
        // This is needed since the Schedule object does not know the current report step of the simulator and
        // we only allow PYACTIONS for the current and future report steps. 
        std::size_t current_report_step = 0;
        // The simUpdateFromPython points to a SimulatorUpdate collecting all updates from one PYACTION call.
        // The SimulatorUpdate is reset before a new PYACTION call is executed.
        // It is a shared_ptr, so a Schedule can be constructed using the copy constructor sharing the simUpdateFromPython.
        // The copy constructor is needed for creating a mocked simulator (msim).
        std::shared_ptr<SimulatorUpdate> simUpdateFromPython{};

        void init_completed_cells_lgr(const EclipseGrid& ecl_grid);
        void load_rst(const RestartIO::RstState& rst,
                      const TracerConfig& tracer_config,
                      const ScheduleGrid& grid,
                      const FieldPropsManager& fp);
        void addWell(Well well);
        void addWell(const std::string& wellName,
                     const std::string& group,
                     int headI,
                     int headJ,
                     Phase preferredPhase,
                     const std::optional<double>& refDepth,
                     double drainageRadius,
                     bool allowCrossFlow,
                     bool automaticShutIn,
                     int pvt_table,
                     WellGasInflowEquation gas_inflow,
                     std::size_t timeStep,
                     ConnectionOrder wellConnectionOrder);
        bool updateWPAVE(const std::string& wname, std::size_t report_step, const PAvg& pavg);

        void updateGuideRateModel(const GuideRateModel& new_model, std::size_t report_step);
        GTNode groupTree(const std::string& root_node, std::size_t report_step, std::size_t level, const std::optional<std::string>& parent_name) const;
        bool checkGroups(const ParseContext& parseContext, ErrorGuard& errors);
        bool updateWellStatus( const std::string& well, std::size_t reportStep, WellStatus status, std::optional<KeywordLocation> = {});
        void addWellToGroup( const std::string& group_name, const std::string& well_name , std::size_t timeStep);
        void iterateScheduleSection(std::size_t load_start,
                                    std::size_t load_end,
                                    const ParseContext& parseContext,
                                    ErrorGuard& errors,
                                    const ScheduleGrid& grid,
                                    const std::unordered_map<std::string, double> * target_wellpi,
                                    const std::string& prefix,
                                    const bool keepKeywords,
                                    const bool log_to_debug = false);
        void addACTIONX(const Action::ActionX& action);
        void addGroupToGroup( const std::string& parent_group, const std::string& child_group);
        void addGroup(const std::string& groupName , std::size_t timeStep);
        void addGroup(Group group);
        void addGroup(const RestartIO::RstGroup& rst_group, std::size_t timeStep);
        void addWell(const std::string& wellName, const DeckRecord& record,
                    std::size_t timeStep, ConnectionOrder connection_order);
        void checkIfAllConnectionsIsShut(std::size_t currentStep);
        void end_report(std::size_t report_step);
        /// \param welsegs_wells All wells with a WELSEGS entry for checks.
        /// \param compegs_wells All wells with a COMPSEGS entry for checks.
        void handleKeyword(std::size_t currentStep,
                           const ScheduleBlock& block,
                           const DeckKeyword& keyword,
                           const ParseContext& parseContext,
                           ErrorGuard& errors,
                           const ScheduleGrid& grid,
                           const Action::Result::MatchingEntities& matches,
                           bool action_mode,
                           SimulatorUpdate* sim_update,
                           const std::unordered_map<std::string, double>* target_wellpi,
                           std::unordered_map<std::string, double>& wpimult_global_factor,
                           WelSegsSet* welsegs_wells = nullptr,
                           std::set<std::string>* compsegs_wells = nullptr);

        void internalWELLSTATUSACTIONXFromPYACTION(const std::string& well_name, std::size_t report_step, const std::string& wellStatus);
        void prefetchPossibleFutureConnections(const ScheduleGrid& grid, const DeckKeyword& keyword,
                                               const ParseContext& parseContext, ErrorGuard& errors);
        void store_wgnames(const DeckKeyword& keyword);
        std::vector<std::string> wellNames(const std::string& pattern,
                                           const HandlerContext& context,
                                           bool allowEmpty = false);
        std::vector<std::string> wellNames(const std::string& pattern, std::size_t timeStep, const std::vector<std::string>& matching_wells, InputErrorAction error_action, ErrorGuard& errors, const KeywordLocation& location) const;
        static std::string formatDate(std::time_t t);
        std::string simulationDays(std::size_t currentStep) const;
        void applyGlobalWPIMULT( const std::unordered_map<std::string, double>& wpimult_global_factor);

        bool must_write_rst_file(std::size_t report_step) const;

        bool isWList(std::size_t report_step, const std::string& pattern) const;

        SimulatorUpdate applyAction(std::size_t reportStep, const std::string& action_name, const std::vector<std::string>& matching_wells);
    };
}

#endif
