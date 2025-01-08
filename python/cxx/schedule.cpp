#include <ctime>
#include <chrono>
#include <map>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/input/eclipse/Parser/InputErrorAction.hpp>
#include <opm/input/eclipse/Parser/ParseContext.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Deck/DeckKeyword.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>

#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>

#include <fmt/format.h>
#include <pybind11/stl.h>
#include <pybind11/chrono.h>
#include "export.hpp"

#include <python/cxx/OpmCommonPythonDoc.hpp>

namespace {

    using system_clock = std::chrono::system_clock;


    /*
      timezones - the stuff that make you wonder why didn't do social science in
      university. The situation here is as follows:

      1. In the C++ code Eclipse style string literals like "20. NOV 2017" are
         converted to time_t values using the utc based function timegm() which
         does not take timezones into account.

      2. Here we use the function gmtime( ) to convert back from a time_t value
         to a broken down struct tm representation.

      3. The broken down representation is then converted to a time_t value
         using the timezone aware function mktime().

      4. The time_t value is converted to a std::chrono::system_clock value.

      Finally std::chrono::system_clock value is automatically converted to a
      python datetime object as part of the pybind11 process. This latter
      conversion *is* timezone aware, that is the reason we must go through
      these hoops.
    */
    system_clock::time_point datetime( std::time_t utc_time) {
        struct tm utc_tm;
        time_t local_time;

        gmtime_r(&utc_time, &utc_tm);
        local_time = mktime(&utc_tm);

        return system_clock::from_time_t(local_time);
    }

    const Well& get_well( const Schedule& sch, const std::string& name, const size_t& timestep ) try {
        return sch.getWell( name, timestep );
    } catch( const std::invalid_argument& e ) {
        throw py::key_error( name );
    }

    double zero_if_undefined(const UDAValue& val) {
        return val.is_numeric() ? val.get<double>() : 0.0;
    };

    std::map<std::string, double> get_production_properties(
        const Schedule& sch, const std::string& well_name, const size_t& report_step)
    {
        const Well* well = nullptr;
        try{
            well = &(sch.getWell( well_name, report_step ));
        } catch( const std::out_of_range& e ) {
            throw py::index_error( fmt::format("well {} is not defined", well_name ));
        }
        if (well->isProducer()) {
            auto& prod_prop = well->getProductionProperties();
            return {
                { "oil_rate", zero_if_undefined(prod_prop.OilRate) },
                { "gas_rate", zero_if_undefined(prod_prop.GasRate) },
                { "water_rate", zero_if_undefined(prod_prop.WaterRate) },
                { "liquid_rate", zero_if_undefined(prod_prop.LiquidRate) },
                { "resv_rate", zero_if_undefined(prod_prop.ResVRate) },
                { "bhp_target", zero_if_undefined(prod_prop.BHPTarget) },
                { "thp_target", zero_if_undefined(prod_prop.THPTarget) },
                { "alq_value", zero_if_undefined(prod_prop.ALQValue) },
            };
        }
        else {
            throw py::key_error( fmt::format("well {} is not a producer", well_name) );
        }
    }

    std::map<std::string, double> get_injection_properties(
        const Schedule& sch, const std::string& well_name, const size_t& report_step)
    {
        const Well* well = nullptr;
        try{
            well = &(sch.getWell( well_name, report_step ));
        } catch( const std::out_of_range& e ) {
            throw py::index_error( fmt::format("well {}: invalid well name", well_name ));
        }

        if (well->isInjector()) {
            auto& inj_prop = well->getInjectionProperties();
            return {
                { "surf_inj_rate", zero_if_undefined(inj_prop.surfaceInjectionRate) },
                { "resv_inj_rate", zero_if_undefined(inj_prop.reservoirInjectionRate) },
                { "bhp_target", zero_if_undefined(inj_prop.BHPTarget) },
                { "thp_target", zero_if_undefined(inj_prop.THPTarget) },
            };
        }
        else {
            throw py::key_error( fmt::format("well {} is not an injector", well_name) );
        }
    }

    system_clock::time_point get_start_time( const Schedule& s ) {
        return datetime(s.posixStartTime());
    }

    system_clock::time_point get_end_time( const Schedule& s ) {
        return datetime(s.posixEndTime());
    }

    std::vector<system_clock::time_point> get_timesteps( const Schedule& s ) {
        std::vector< system_clock::time_point > v;

        for( size_t i = 0; i < s.size(); ++i )
            v.push_back( datetime( std::chrono::system_clock::to_time_t(s[i].start_time() )));

        return v;
    }

    std::vector<Group> get_groups( const Schedule& sch, size_t timestep ) {
        std::vector< Group > groups;
        for( const auto& group_name : sch.groupNames())
            groups.push_back( sch.getGroup(group_name, timestep) );

        return groups;
    }

    bool has_well( const Schedule& sch, const std::string& wellName) {
        return sch.hasWell( wellName );
    }

    const ScheduleState& getitem(const Schedule& sch, std::size_t report_step) {
        return sch[report_step];
    }

    std::vector<std::unique_ptr<DeckKeyword>> parseKeywords(const std::string& deck_string, const UnitSystem& unit_system)
    {
        Parser parser;
        ParseContext parseContext;
        // The next line will suppress PARSE_INVALID_KEYWORD_COMBINATION errors. This is needed when a keyword requires other keywords, since the
        // required keywords might be in the original .DATA file, to which we do not have access here.
        // When inserting keywords that can cause a PARSE_INVALID_KEYWORD_COMBINATION error, i.e. keywords with required or prohibited keywords,
        // a warning will be added to the OpmLog.
        parseContext.update(ParseContext::PARSE_INVALID_KEYWORD_COMBINATION , InputErrorAction::IGNORE );
        std::string str {unit_system.deck_name() + "\n\n" + deck_string};
        auto deck = parser.parseString(str, parseContext);
        std::vector<std::unique_ptr<DeckKeyword>> keywords;
        for (auto &keyword : deck) {
            auto parserKeyword = parser.getKeyword(keyword.name());
            std::for_each(parserKeyword.requiredKeywords().begin(), parserKeyword.requiredKeywords().end(), [&keyword](const auto& requiredKeyword) {
                Opm::OpmLog::warning("Attention, the keyword " + keyword.name() + " needs the keywords " + requiredKeyword + " before.");
            });
            std::for_each(parserKeyword.prohibitedKeywords().begin(), parserKeyword.prohibitedKeywords().end(), [&keyword](const auto& prohibitedKeyword) {
                Opm::OpmLog::warning("Attention, the keyword " + keyword.name() + " is incompatible with the keyword " + prohibitedKeyword + ".");
            });
            keywords.push_back(std::make_unique<DeckKeyword>(keyword));
        }
        return keywords;
    }
    void insert_keywords(Schedule& sch, const std::string& deck_string, std::size_t report_step, const UnitSystem& unit_system)
    {
        auto kws = parseKeywords(deck_string, unit_system);
        sch.applyKeywords(kws, report_step);
    }

    void insert_keywords(Schedule& sch, const std::string& deck_string, std::size_t report_step)
    {
        auto kws = parseKeywords(deck_string,sch.getUnits());
        sch.applyKeywords(kws, report_step);
    }

    void insert_keywords(Schedule& sch, const std::string& deck_string)
    {
        auto kws = parseKeywords(deck_string,sch.getUnits());
        sch.applyKeywords(kws);
    }

    // NOTE: this overload does currently not work, see PR #2833. The plan
    //  is to fix this in a later commit. For now, the overload insert_keywords()
    //  above taking a deck_string (std::string) instead of a list of DeckKeywords
    //  has to be used instead.
    void insert_keywords(
        Schedule& sch, py::list& deck_keywords, std::size_t report_step)
    {
        std::vector<std::unique_ptr<DeckKeyword>> keywords;
        for (py::handle item : deck_keywords) {
            DeckKeyword &keyword = item.cast<DeckKeyword&>();
            keywords.push_back(std::make_unique<DeckKeyword>(keyword));
        }
        sch.applyKeywords(keywords, report_step);
    }
}



void python::common::export_Schedule(py::module& module) {

    using namespace Opm::Common::DocStrings;

    // Note: In the below class we use std::shared_ptr as the holder type, see:
    //
    //  https://pybind11.readthedocs.io/en/stable/advanced/smart_ptrs.html
    //
    // this makes it possible to share the returned object with e.g. and
    //   opm.simulators.BlackOilSimulator Python object
    //
    py::class_< Schedule, std::shared_ptr<Schedule> >( module, "Schedule", ScheduleClass_docstring)
    .def(py::init<const Deck&, const EclipseState& >(), py::arg("deck"), py::arg("eclipse_state"))
    .def("_groups", &get_groups, py::arg("report_step"), Schedule_groups_docstring)
    .def_property_readonly( "start",  &get_start_time )
    .def_property_readonly( "end",    &get_end_time )
    .def_property_readonly( "timesteps", &get_timesteps )
    .def("__len__", &Schedule::size)
    .def("__getitem__", &getitem, py::arg("report_step"), Schedule_getitem_docstring)
    .def("shut_well", py::overload_cast<const std::string&, std::size_t>(&Schedule::shut_well), py::arg("well_name"), py::arg("step"), Schedule_shut_well_well_name_step_docstring)
    .def("shut_well", py::overload_cast<const std::string&>(&Schedule::shut_well), py::arg("well_name"), Schedule_shut_well_well_name_docstring)
    .def("open_well", py::overload_cast<const std::string&, std::size_t>(&Schedule::open_well), py::arg("well_name"), py::arg("step"), Schedule_open_well_well_name_step_docstring)
    .def("open_well", py::overload_cast<const std::string&>(&Schedule::open_well), py::arg("well_name"), Schedule_open_well_well_name_docstring)
    .def("stop_well", py::overload_cast<const std::string&, std::size_t>(&Schedule::stop_well), py::arg("well_name"), py::arg("step"), Schedule_stop_well_well_name_step_docstring)
    .def("stop_well", py::overload_cast<const std::string&>(&Schedule::stop_well), py::arg("well_name"), Schedule_stop_well_well_name_docstring)
    .def( "get_wells", &Schedule::getWells, py::arg("well_name_pattern"), Schedule_get_wells_docstring)
    .def( "get_injection_properties", &get_injection_properties, py::arg("well_name"), py::arg("report_step"), Schedule_get_injection_properties_docstring)
    .def( "get_production_properties", &get_production_properties, py::arg("well_name"), py::arg("report_step"), Schedule_get_production_properties_docstring)
    .def("well_names", py::overload_cast<const std::string&>(&Schedule::wellNames, py::const_), py::arg("well_name_pattern"))
    .def( "get_well", &get_well, py::arg("well_name"), py::arg("report_step"), Schedule_get_well_docstring)
    .def( "insert_keywords", py::overload_cast<Schedule&, py::list&, std::size_t>(&insert_keywords), py::arg("keywords"), py::arg("step"))
    .def( "insert_keywords", py::overload_cast<Schedule&, const std::string&, std::size_t, const UnitSystem&>(&insert_keywords), py::arg("data"), py::arg("step"), py::arg("unit_system"))
    .def( "insert_keywords", py::overload_cast<Schedule&, const std::string&, std::size_t>(&insert_keywords), py::arg("data"), py::arg("step"))
    .def( "insert_keywords", py::overload_cast<Schedule&, const std::string&>(&insert_keywords),py::arg("data"))
    .def("__contains__", &has_well, py::arg("well_name"))
    ;
}
