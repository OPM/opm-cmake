/*
  Copyright 2019  Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify it under the terms
  of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

  OPM is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
  A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along with
  OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <opm/input/eclipse/EclipseState/Grid/FieldProps.hpp>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/LogUtil.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/common/utility/OpmInputError.hpp>

#include <opm/input/eclipse/EclipseState/Aquifer/NumericalAquifer/NumericalAquifers.hpp>
#include <opm/input/eclipse/EclipseState/Grid/Box.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FieldPropsManager.hpp> // Layering violation.  Needed for apply_tran() function.
#include <opm/input/eclipse/EclipseState/Grid/Keywords.hpp>
#include <opm/input/eclipse/EclipseState/Grid/SatfuncPropertyInitializers.hpp>
#include <opm/input/eclipse/EclipseState/Runspec.hpp>
#include <opm/input/eclipse/EclipseState/Tables/RtempvdTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/input/eclipse/EclipseState/Util/OrderedMap.hpp>

#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>

#include <opm/input/eclipse/Parser/ParserKeywords/A.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/B.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/C.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/E.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/M.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/O.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/P.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/T.hpp>

#include "Operate.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <functional>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <fmt/format.h>

namespace {
    Opm::Box makeGlobalGridBox(const Opm::EclipseGrid* gridPtr,
                               const std::vector<int>* actnum = nullptr,
                               const std::unordered_map<int, int>* index = nullptr)
    {
        return Opm::Box {
            *gridPtr,
            [gridPtr, actnum](const std::size_t global_index)
            {
                if (!actnum || actnum->empty()) {
                    return gridPtr->cellActive(global_index);
                }
                if (global_index >= actnum->size()) {
                    return false;
                } else {
                    return (*actnum)[global_index] > 0;
                }
            },
            [gridPtr, index](const std::size_t global_index) -> std::size_t
            {
                if (!index || index->empty()) {
                    return gridPtr->activeIndex(global_index);
                }

                const auto it = index->find(global_index);
                assert(it != index->end());
                return it->second;
            }
        };
    }

    bool is_capillary_pressure(const std::string& keyword)
    {
        return (keyword == "PCW")  || (keyword == "PCG")
            || (keyword == "IPCG") || (keyword == "IPCW");
    }
} // Anonymous namespace

namespace Opm {

namespace Fieldprops { namespace keywords {

namespace {
    std::string get_keyword_from_alias(const std::string& name)
    {
        auto kwPos = ALIAS::aliased_keywords.find(name);

        return (kwPos == ALIAS::aliased_keywords.end())
            ? name : kwPos->second;
    }
} // Anonymous namespace

static const std::set<std::string> oper_keywords = {"ADD", "EQUALS", "MAXVALUE", "MINVALUE", "MULTIPLY"};
static const std::set<std::string> region_oper_keywords = {"MULTIREG", "ADDREG", "EQUALREG", "OPERATER"};
static const std::set<std::string> box_keywords = {"BOX", "ENDBOX"};

bool is_oper_keyword(const std::string& name)
{
    return (oper_keywords.find(name) != oper_keywords.end()
            || region_oper_keywords.find(name) != region_oper_keywords.end());
}

template <>
keyword_info<double>
global_kw_info(const std::string& name, const bool allow_unsupported)
{
    if (auto kwPos = GRID::double_keywords.find(name);
        kwPos != GRID::double_keywords.end())
    {
        return kwPos->second;
    }

    if (auto kwPos = EDIT::double_keywords.find(name);
        kwPos != EDIT::double_keywords.end())
    {
        return kwPos->second;
    }

    if (auto kwPos = PROPS::double_keywords.find(name);
        kwPos != PROPS::double_keywords.end())
    {
        return kwPos->second;
    }

    if (PROPS::satfunc.count(name)) {
        return keyword_info<double>{};
    }

    if (auto kwPos = SOLUTION::double_keywords.find(name);
        kwPos != SOLUTION::double_keywords.end())
    {
        return kwPos->second;
    }
    
    if (auto kwPos = SOLUTION::composition_keywords.find(name);
        kwPos != SOLUTION::composition_keywords.end())
    {
        return kwPos->second;
    }

    if (auto kwPos = SCHEDULE::double_keywords.find(name);
        kwPos != SCHEDULE::double_keywords.end())
    {
        return kwPos->second;
    }

    if (allow_unsupported) {
        return keyword_info<double>{};
    }

    throw std::out_of_range {
        fmt::format("INFO: '{}' is not a double precision property.", name)
    };
}

template <>
keyword_info<int>
global_kw_info(const std::string& name, bool)
{
    if (auto kwPos = GRID::int_keywords.find(name);
        kwPos != GRID::int_keywords.end())
    {
        return kwPos->second;
    }

    if (auto kwPos = EDIT::int_keywords.find(name);
        kwPos != EDIT::int_keywords.end())
    {
        return kwPos->second;
    }

    if (auto kwPos = PROPS::int_keywords.find(name);
        kwPos != PROPS::int_keywords.end())
    {
        return kwPos->second;
    }

    if (auto kwPos = REGIONS::int_keywords.find(name);
        kwPos != REGIONS::int_keywords.end())
    {
        return kwPos->second;
    }

    if (auto kwPos = SCHEDULE::int_keywords.find(name);
        kwPos != SCHEDULE::int_keywords.end())
    {
        return kwPos->second;
    }

    if (isFipxxx(name)) {
        return keyword_info<int>{}.init(1);
    }

    throw std::out_of_range {
        fmt::format("INFO: '{}' is not an integer property", name)
    };
}

}} // namespace Fieldprops::keywords

namespace {

// The EQUALREG, MULTREG, COPYREG, ... keywords are used to manipulate
// vectors based on region values; for instance the statement
//
//   EQUALREG
//      PORO  0.25  3    /   -- Region array not specified
//      PERMX 100   3  F /
//   /
//
// will set the PORO field to 0.25 for all cells in region 3 and the PERMX
// value to 100 mD for the same cells.  The fourth optional argument to the
// EQUALREG keyword is used to indicate which REGION array should be used
// for the selection.
//
// If the REGION array is not indicated (as in the PORO case) above, the
// default region to use in the xxxREG keywords depends on the GRIDOPTS
// keyword:
//
//   1. If GRIDOPTS is present, and the NRMULT item is greater than zero,
//      the xxxREG keywords will default to use the MULTNUM region.
//
//   2. If the GRIDOPTS keyword is not present - or the NRMULT item equals
//      zero, the xxxREG keywords will default to use the FLUXNUM keyword.
//
// This quite weird behaviour comes from reading the GRIDOPTS and MULTNUM
// documentation, and practical experience with ECLIPSE simulations.
// Unfortunately the documentation of the xxxREG keywords does not confirm
// this.
std::string default_region_keyword(const Deck& deck)
{
    if (deck.hasKeyword("GRIDOPTS")) {
        const auto& gridOpts = deck["GRIDOPTS"].back();
        const auto& record = gridOpts.getRecord(0);
        const auto& nrmult_item = record.getItem("NRMULT");

        if (nrmult_item.get<int>(0) > 0) {
            return "MULTNUM"; // GRIDOPTS and positive NRMULT
        }
    }

    return "FLUXNUM";
}

void reject_undefined_operation(const KeywordLocation& loc,
                                const std::size_t      numUnInit,
                                const std::size_t      numElements,
                                std::string_view       operation,
                                std::string_view       arrayName)
{
    const auto* plural = (numElements > 1) ? "s" : "";

    throw OpmInputError {
        fmt::format
        (R"({0} operation on array {1} would
generate an undefined result in {2} out of {3} block{4}.)",
         operation, arrayName, numUnInit, numElements, plural),
        loc
    };
}

template <typename T>
void verify_deck_data(const Fieldprops::keywords::keyword_info<T>& kw_info,
                      const DeckKeyword&                           keyword,
                      const std::vector<T>&                        deck_data,
                      const Box&                                   box)
{
    // there can be multiple values for each grid cell
    if (box.size() * kw_info.num_value != deck_data.size()) {
        const auto& location = keyword.location();
        std::string msg = "Fundamental error with keyword: " + keyword.name() +
            " at: " + location.filename + ", line: " + std::to_string(location.lineno) +
            " got " + std::to_string(deck_data.size()) + " elements - expected : " + std::to_string(box.size() * kw_info.num_value);
        throw std::invalid_argument(msg);
    }
}

void log_empty_region(const DeckKeyword& keyword,
                      const std::string& region_name,
                      const int          region_id,
                      const std::string& array_name)
{
    const auto message =
        fmt::format(R"(Region {} of {} has no active cells when processing operation {} on array {}.
Please check whether this is on purpose or if you did not properly define this region set.)",
                    region_id, region_name, keyword.name(), array_name);

    OpmLog::warning(Log::fileMessage(keyword.location(), message));
}

template <typename T>
void assign_deck(const Fieldprops::keywords::keyword_info<T>& kw_info,
                 const DeckKeyword& keyword,
                 Fieldprops::FieldData<T>& field_data,
                 const std::vector<T>& deck_data,
                 const std::vector<value::status>& deck_status,
                 const Box& box)
{
    verify_deck_data(kw_info, keyword, deck_data, box);

    for (const auto& cell_index : box.index_list()) {
        auto active_index = cell_index.active_index;
        auto data_index = cell_index.data_index;
        for (size_t i = 0; i < kw_info.num_value; ++i) {
            auto deck_data_index = i* box.size() + data_index;
            if (value::has_value(deck_status[deck_data_index])) {
                auto data_active_index = i * box.size() + active_index;
                if (deck_status[deck_data_index] == value::status::deck_value ||
                    field_data.value_status[data_active_index] == value::status::uninitialized) {
                    field_data.data[data_active_index] = deck_data[deck_data_index];
                    field_data.value_status[data_active_index] = deck_status[deck_data_index];
                }
            }
        }
    }

    if (kw_info.global) {
        auto& global_data = field_data.global_data.value();
        auto& global_status = field_data.global_value_status.value();
        const auto& index_list = box.global_index_list();

        for (const auto& cell : index_list) {
            if ((deck_status[cell.data_index] == value::status::deck_value) ||
                (global_status[cell.global_index] == value::status::uninitialized))
            {
                global_data[cell.global_index] = deck_data[cell.data_index];
                global_status[cell.global_index] = deck_status[cell.data_index];
            }
        }
    }
}

template <typename T>
void multiply_deck(const Fieldprops::keywords::keyword_info<T>& kw_info,
                   const DeckKeyword& keyword,
                   Fieldprops::FieldData<T>& field_data,
                   const std::vector<T>& deck_data,
                   const std::vector<value::status>& deck_status,
                   const Box& box)
{
    verify_deck_data(kw_info, keyword, deck_data, box);
    for (const auto& cell_index : box.index_list()) {
        auto active_index = cell_index.active_index;
        auto data_index = cell_index.data_index;

        if (value::has_value(deck_status[data_index]) &&
            value::has_value(field_data.value_status[active_index]))
        {
            field_data.data[active_index] *= deck_data[data_index];
            field_data.value_status[active_index] = deck_status[data_index];
        }
    }

    if (kw_info.global) {
        auto& global_data = field_data.global_data.value();
        auto& global_status = field_data.global_value_status.value();
        const auto& index_list = box.global_index_list();

        for (const auto& cell : index_list) {
            if ((deck_status[cell.data_index] == value::status::deck_value) ||
                (global_status[cell.global_index] == value::status::uninitialized))
            {
                global_data[cell.global_index] *= deck_data[cell.data_index];
                global_status[cell.global_index] = deck_status[cell.data_index];
            }
        }
    }
}

template <typename T>
void assign_scalar(std::vector<T>&                     data,
                   std::vector<value::status>&         value_status,
                   const T                             value,
                   const std::vector<Box::cell_index>& index_list)
{
    for (const auto& cell_index : index_list) {
        data[cell_index.active_index] = value;
        value_status[cell_index.active_index] = value::status::deck_value;
    }
}

template <typename T>
void multiply_scalar(const KeywordLocation&              loc,
                     std::string_view                    arrayName,
                     std::vector<T>&                     data,
                     std::vector<value::status>&         value_status,
                     const T                             value,
                     const std::vector<Box::cell_index>& index_list)
{
    auto unInit = 0;

    for (const auto& cell_index : index_list) {
        const auto ix = cell_index.active_index;

        if (value::has_value(value_status[ix])) {
            data[ix] *= value;
        }
        else {
            ++unInit;
        }
    }

    if (unInit > 0) {
        reject_undefined_operation(loc, unInit,
                                   index_list.size(),
                                   "Multiplication", arrayName);
    }
}

template <typename T>
void add_scalar(const KeywordLocation&              loc,
                std::string_view                    arrayName,
                std::vector<T>&                     data,
                std::vector<value::status>&         value_status,
                const T                             value,
                const std::vector<Box::cell_index>& index_list)
{
    auto unInit = 0;

    for (const auto& cell_index : index_list) {
        const auto ix = cell_index.active_index;

        if (value::has_value(value_status[ix])) {
            data[ix] += value;
        }
        else {
            ++unInit;
        }
    }

    if (unInit > 0) {
        reject_undefined_operation(loc, unInit,
                                   index_list.size(),
                                   "Addition", arrayName);
    }
}

template <typename T>
void min_value(const KeywordLocation&              loc,
               std::string_view                    arrayName,
               std::vector<T>&                     data,
               std::vector<value::status>&         value_status,
               const T                             value,
               const std::vector<Box::cell_index>& index_list)
{
    auto unInit = 0;

    for (const auto& cell_index : index_list) {
        const auto ix = cell_index.active_index;

        if (value::has_value(value_status[ix])) {
            data[ix] = std::max(data[ix], value);
        }
        else {
            ++unInit;
        }
    }

    if (unInit > 0) {
        reject_undefined_operation(loc, unInit,
                                   index_list.size(),
                                   "Minimum threshold",
                                   arrayName);
    }
}

template <typename T>
void max_value(const KeywordLocation&              loc,
               std::string_view                    arrayName,
               std::vector<T>&                     data,
               std::vector<value::status>&         value_status,
               const T                             value,
               const std::vector<Box::cell_index>& index_list)
{
    auto unInit = 0;

    for (const auto& cell_index : index_list) {
        const auto ix = cell_index.active_index;

        if (value::has_value(value_status[ix])) {
            data[ix] = std::min(data[ix], value);
        }
        else {
            ++unInit;
        }
    }

    if (unInit > 0) {
        reject_undefined_operation(loc, unInit,
                                   index_list.size(),
                                   "Maximum threshold",
                                   arrayName);
    }
}


template<typename T>
void update_global_from_local(Fieldprops::FieldData<T>& data,
                              const std::vector<Box::cell_index>& index_list)
{
    if(data.global_data)
    {
        auto& to = *data.global_data;
        auto to_st = *data.global_value_status;
        const auto& from = data.data;
        const auto& from_st = data.value_status;

        for (const auto& cell_index : index_list) {
            to[cell_index.global_index] = from[cell_index.active_index];
            to_st[cell_index.global_index] = from_st[cell_index.active_index];
        }
    }
}

template <typename T>
void apply(const Fieldprops::ScalarOperation   op,
           const KeywordLocation&              loc,
           std::string_view                    arrayName,
           std::vector<T>&                     data,
           std::vector<value::status>&         value_status,
           const T                             scalar_value,
           const std::vector<Box::cell_index>& index_list)
{
    switch (op) {
    case Fieldprops::ScalarOperation::EQUAL:
        assign_scalar(data, value_status, scalar_value, index_list);
        return;

    case Fieldprops::ScalarOperation::MUL:
        multiply_scalar(loc, arrayName, data, value_status, scalar_value, index_list);
        return;

    case Fieldprops::ScalarOperation::ADD:
        add_scalar(loc, arrayName, data, value_status, scalar_value, index_list);
        return;

    case Fieldprops::ScalarOperation::MIN:
        min_value(loc, arrayName, data, value_status, scalar_value, index_list);
        return;

    case Fieldprops::ScalarOperation::MAX:
        max_value(loc, arrayName, data, value_status, scalar_value, index_list);
        return;
    }

    throw std::invalid_argument {
        fmt::format("'{}' is not a known operation.", static_cast<int>(op))
    };
}

std::string make_region_name(const std::string& deck_value)
{
    if (deck_value == "O") { return "OPERNUM"; }
    if (deck_value == "F") { return "FLUXNUM"; }
    if (deck_value == "M") { return "MULTNUM"; }

    throw std::invalid_argument {
        fmt::format("Input string '{}' is not a valid "
                    "region set name. Expected 'O'/'F'/'M'", deck_value)
    };
}

Fieldprops::ScalarOperation fromString(const std::string& keyword)
{
    if ((keyword == ParserKeywords::ADD::keywordName) ||
        (keyword == ParserKeywords::ADDREG::keywordName))
    {
        return Fieldprops::ScalarOperation::ADD;
    }

    if ((keyword == ParserKeywords::EQUALS::keywordName) ||
        (keyword == ParserKeywords::EQUALREG::keywordName))
    {
        return Fieldprops::ScalarOperation::EQUAL;
    }

    if ((keyword == ParserKeywords::MULTIPLY::keywordName) ||
        (keyword == ParserKeywords::MULTIREG::keywordName))
    {
        return Fieldprops::ScalarOperation::MUL;
    }

    if (keyword == ParserKeywords::MINVALUE::keywordName) {
        return Fieldprops::ScalarOperation::MIN;
    }

    if (keyword == ParserKeywords::MAXVALUE::keywordName) {
        return Fieldprops::ScalarOperation::MAX;
    }

    throw std::invalid_argument {
        fmt::format("Keyword operation ({}) not recognized", keyword)
    };
}

void handle_box_keyword(const DeckKeyword& deckKeyword, Box& box)
{
    if (deckKeyword.name() == ParserKeywords::BOX::keywordName) {
        const auto& record = deckKeyword.getRecord(0);
        box.update(record);
    }
    else {
        box.reset();
    }
}

std::vector<double> extract_cell_volume(const EclipseGrid& grid)
{
    return grid.activeVolume();
}

std::vector<double> extract_cell_depth(const EclipseGrid& grid)
{
    std::vector<double> cell_depth(grid.getNumActive());

    for (std::size_t active_index = 0; active_index < grid.getNumActive(); ++active_index) {
        cell_depth[active_index] = grid.getCellDepth(grid.getGlobalIndex(active_index));
    }

    return cell_depth;
}

// The rst_compare_data function compares the main std::map<std::string,
// std::vector<T>> data containers. If one of the containers contains a keyword
// *which is fully defaulted* and the other container does not contain said
// keyword - the containers are considered to be equal.
template <typename T>
bool rst_compare_data(const std::unordered_map<std::string, Fieldprops::FieldData<T>>& data1,
                      const std::unordered_map<std::string, Fieldprops::FieldData<T>>& data2)
{
    std::unordered_set<std::string> keys;
    for (const auto& [key, _] : data1) {
        (void)_;
        keys.insert(key);
    }

    for (const auto& [key, _] : data2) {
        (void)_;
        keys.insert(key);
    }

    for (const auto& key : keys) {
        const auto& d1 = data1.find(key);
        const auto& d2 = data2.find(key);

        if (d1 == data1.end()) {
            if (!d2->second.valid_default())
                return false;
            continue;
        }

        if (d2 == data2.end()) {
            if (!d1->second.valid_default())
                return false;
            continue;
        }

        if (!(d1->second == d2->second))
            return false;
    }

    return true;
}

} // Anonymous namespace

bool FieldProps::operator==(const FieldProps& other) const
{
    return (this->unit_system == other.unit_system)
        && (this->nx == other.nx)
        && (this->ny == other.ny)
        && (this->nz == other.nz)
        && (this->m_phases == other.m_phases)
        && (this->m_satfuncctrl == other.m_satfuncctrl)
        && (this->m_actnum == other.m_actnum)
        && (this->cell_volume == other.cell_volume)
        && (this->cell_depth == other.cell_depth)
        && (this->m_default_region == other.m_default_region)
        && (this->m_rtep == other.m_rtep)
        && (this->tables == other.tables)
        && (this->multregp == other.multregp)
        && (this->int_data == other.int_data)
        && (this->double_data == other.double_data)
        && (this->fipreg_shortname_translation == other.fipreg_shortname_translation)
        && (this->tran == other.tran)
        ;
}

bool FieldProps::rst_cmp(const FieldProps& full_arg, const FieldProps& rst_arg)
{
    if (!rst_compare_data(full_arg.double_data, rst_arg.double_data)) {
        return false;
    }

    if (!rst_compare_data(full_arg.int_data, rst_arg.int_data)) {
        return false;
    }

    if (!UnitSystem::rst_cmp(full_arg.unit_system, rst_arg.unit_system)) {
        return false;
    }

    return (full_arg.nx == rst_arg.nx)
        && (full_arg.ny == rst_arg.ny)
        && (full_arg.nz == rst_arg.nz)
        && (full_arg.m_phases == rst_arg.m_phases)
        && (full_arg.m_satfuncctrl == rst_arg.m_satfuncctrl)
        && (full_arg.m_actnum == rst_arg.m_actnum)
        && (full_arg.cell_volume == rst_arg.cell_volume)
        && (full_arg.cell_depth == rst_arg.cell_depth)
        && (full_arg.m_default_region == rst_arg.m_default_region)
        && (full_arg.m_rtep == rst_arg.m_rtep)
        && (full_arg.tables == rst_arg.tables)
        && (full_arg.multregp == rst_arg.multregp)
        && (full_arg.tran == rst_arg.tran)
        ;
}

// init_get methods have to be specialized before their instantiation in the
// constructor below. Otherwise we get a compilation error.
template <>
Fieldprops::FieldData<double>&
FieldProps::init_get(const std::string& keyword_name,
                     const Fieldprops::keywords::keyword_info<double>& kw_info,
                     const bool multiplier_in_edit)
{
    if (multiplier_in_edit && !kw_info.scalar_init.has_value()) {
        OPM_THROW(std::logic_error, "Keyword " +  keyword_name +
                  " is a multiplier and should have a default initial value.");
    }

    const auto keyword = Fieldprops::keywords::get_keyword_from_alias(keyword_name);
    const auto mult_keyword = std::string(multiplier_in_edit ? getMultiplierPrefix() : "") + keyword;

    auto iter = this->double_data.find(mult_keyword);
    if (iter != this->double_data.end()) {
        return iter->second;
    }
    else if (multiplier_in_edit) {
        assert(keyword != ParserKeywords::PORV::keywordName);
        assert(keyword != ParserKeywords::TEMPI::keywordName);
        assert(Fieldprops::keywords::PROPS::satfunc.count(keyword) != 1);
        assert(!is_capillary_pressure(keyword));

        this->multiplier_kw_infos_.insert_or_assign(mult_keyword, kw_info);
    }

    auto elmDescr = this->double_data
        .try_emplace(mult_keyword, kw_info, this->active_size,
                     kw_info.global ? this->global_size : std::size_t{0});

    auto& propData = elmDescr.first->second;

    if (keyword == ParserKeywords::PORV::keywordName) {
        this->init_porv(propData);
    }

    if (keyword == ParserKeywords::TEMPI::keywordName) {
        this->init_tempi(propData);
    }

    if ((Fieldprops::keywords::PROPS::satfunc.count(keyword) == 1) ||
        is_capillary_pressure(keyword))
    {
        this->init_satfunc(keyword, propData);
    }

    return propData;
}

template <>
Fieldprops::FieldData<double>&
FieldProps::init_get(const std::string& keyword, const bool allow_unsupported)
{
    return this->init_get(keyword, Fieldprops::keywords::global_kw_info<double>(keyword, allow_unsupported));
}

template <>
Fieldprops::FieldData<int>&
FieldProps::init_get(const std::string&                             keyword,
                     const Fieldprops::keywords::keyword_info<int>& kw_info,
                     const bool)
{
    auto iter = this->int_data.find(keyword);
    if (iter != this->int_data.end()) {
        return iter->second;
    }

    return this->int_data
        .try_emplace(keyword, kw_info, this->active_size,
                     kw_info.global ? this->global_size : 0).first->second;
}

template <>
Fieldprops::FieldData<int>&
FieldProps::init_get(const std::string& keyword, bool)
{
    if (Fieldprops::keywords::isFipxxx(keyword)) {
        auto kw_info = Fieldprops::keywords::keyword_info<int>{};
        kw_info.init(1);

        return this->init_get(this->canonical_fipreg_name(keyword), kw_info);
    }

    return this->init_get(keyword, Fieldprops::keywords::global_kw_info<int>(keyword));
}


FieldProps::FieldProps(const Deck& deck,
                       const Phases& phases,
                       EclipseGrid& grid,
                       const TableManager& tables_arg,
                       const std::size_t ncomps)
    : active_size(grid.getNumActive())
    , global_size(grid.getCartesianSize())
    , unit_system(deck.getActiveUnitSystem())
    , nx(grid.getNX())
    , ny(grid.getNY())
    , nz(grid.getNZ())
    , m_phases(phases)
    , m_satfuncctrl(deck)
    , m_actnum(grid.getACTNUM())
    , cell_volume(extract_cell_volume(grid))
    , cell_depth(extract_cell_depth(grid))
    , m_default_region(default_region_keyword(deck))
    , grid_ptr(&grid)
    , tables(tables_arg)
{
    this->tran.emplace("TRANX", "TRANX");
    this->tran.emplace("TRANY", "TRANY");
    this->tran.emplace("TRANZ", "TRANZ");

    if (deck.hasKeyword<ParserKeywords::MULTREGP>()) {
        this->processMULTREGP(deck);
    }

    if (DeckSection::hasGRID(deck)) {
        if (grid.getMinpvMode() == MinpvMode::EclSTD) {
            // Intial values were set via MINPV/MINPORV in the grid
            // Set them here, too.
            auto kw_info = Fieldprops::keywords::global_kw_info<double>("MINPVV");
            auto& data = this->init_get<double>("MINPVV", kw_info);
            data.default_assign_global(grid.getMinpvVector());
            data.update_local_from_global([&grid](std::size_t i){ return grid.getGlobalIndex(i);});
        }
        this->scanGRIDSection(GRIDSection(deck));
    }

    if (DeckSection::hasEDIT(deck)) {
        this->scanEDITSection(EDITSection(deck));
    }

    grid.resetACTNUM(this->actnum());
    this->reset_actnum(grid.getACTNUM());

    if (DeckSection::hasREGIONS(deck)) {
        this->scanREGIONSSection(REGIONSSection(deck));
    }

    // Update PVTNUM/SATNUM for numerical aquifer cells
    {
        const auto& aqcell_tabnums = grid.getAquiferCellTabnums();

        const bool has_pvtnum = this->int_data.count("PVTNUM") != 0;
        const bool has_satnum = this->int_data.count("SATNUM") != 0;

        std::vector<int>* pvtnum = has_pvtnum ? &(this->int_data["PVTNUM"].data) : nullptr;
        std::vector<int>* satnum = has_satnum ? &(this->int_data["SATNUM"].data) : nullptr;
        for (const auto& [globCell, regionID] : aqcell_tabnums) {
            const auto aix = grid.activeIndex(globCell);
            if (has_pvtnum) { (*pvtnum)[aix] = std::max(regionID[0], (*pvtnum)[aix]); }
            if (has_satnum) { (*satnum)[aix] = std::max(regionID[1], (*satnum)[aix]); }
        }
    }

    if (DeckSection::hasPROPS(deck)) {
        this->scanPROPSSection(PROPSSection(deck));
    }

    if (DeckSection::hasSOLUTION(deck)) {
        this->scanSOLUTIONSection(SOLUTIONSection(deck), ncomps);
    }
}


// Special constructor ONLY used to get the correct ACTNUM.
// The grid argument should have all active cells.
FieldProps::FieldProps(const Deck& deck, const EclipseGrid& grid)
    : active_size(grid.getNumActive())
    , global_size(grid.getCartesianSize())
    , unit_system(deck.getActiveUnitSystem())
    , nx(grid.getNX())
    , ny(grid.getNY())
    , nz(grid.getNZ())
    , m_phases()
    , m_satfuncctrl(deck)
    , m_actnum(global_size, 1)  // NB! activates all at start!
    , cell_volume()             // NB! empty for this purpose.
    , cell_depth()              // NB! empty for this purpose.
    , m_default_region(default_region_keyword(deck))
    , grid_ptr(&grid)
    , tables()                  // NB! empty for this purpose.
{
    if (this->active_size != this->global_size) {
        throw std::logic_error {
            "Programmer error: FieldProps special case processing for ACTNUM "
            "called with grid object that already had deactivated cells."
        };
    }

    if (DeckSection::hasGRID(deck)) {
        this->scanGRIDSectionOnlyACTNUM(GRIDSection(deck));
    }
}

void FieldProps::deleteMINPVV()
{
    double_data.erase("MINPVV");
}

void FieldProps::reset_actnum(const std::vector<int>& new_actnum)
{
    if (this->global_size != new_actnum.size()) {
        throw std::logic_error {
            "reset_actnum() must be called with "
            "the same number of global cells"
        };
    }

    if (new_actnum == this->m_actnum) {
        return;
    }

    std::vector<bool> active_map(this->active_size, true);
    std::size_t active_index = 0;
    std::size_t new_active_size = 0;
    for (std::size_t g = 0; g < this->m_actnum.size(); ++g) {
        if (this->m_actnum[g] != 0) {
            if (new_actnum[g] == 0) {
                active_map[active_index] = false;
            }
            else {
                new_active_size += 1;
            }

            active_index += 1;
        }
        else if (new_actnum[g] != 0) {
            throw std::logic_error {
                "It is not possible to activate cells"
            };
        }
    }

    for (auto& data : this->double_data) {
        data.second.compress(active_map);
    }

    for (auto& data : this->int_data) {
        data.second.compress(active_map);
    }

    Fieldprops::compress(this->cell_volume, active_map);
    Fieldprops::compress(this->cell_depth, active_map);

    this->m_actnum = std::move(new_actnum);
    this->active_size = new_active_size;
}

void FieldProps::prune_global_for_schedule_run()
{
    for (auto& data : this->double_data) {
        if(data.second.kw_info.local_in_schedule) {
            data.second.global_data.reset();
            data.second.global_value_status.reset();
        }
    }

    for (auto&  data : this->int_data) {
        if(data.second.kw_info.local_in_schedule) {
            data.second.global_data.reset();
            data.second.global_value_status.reset();
        }
    }
}
void FieldProps::distribute_toplayer(Fieldprops::FieldData<double>& field_data,
                                     const std::vector<double>& deck_data,
                                     const Box& box)
{
    const std::size_t layer_size = this->nx * this->ny;
    Fieldprops::FieldData<double> toplayer(field_data.kw_info, layer_size, 0);
    for (const auto& cell_index : box.index_list()) {
        if (cell_index.global_index < layer_size) {
            toplayer.data[cell_index.global_index] = deck_data[cell_index.data_index];
            toplayer.value_status[cell_index.global_index] = value::status::deck_value;
        }
    }

    std::size_t active_index = 0;
    for (std::size_t k = 0; k < this->nz; k++) {
        for (std::size_t j = 0; j < this->ny; j++) {
            for (std::size_t i = 0; i < this->nx; i++) {
                std::size_t g = i + j*this->nx + k*this->nx*this->ny;
                if (this->m_actnum[g]) {
                    if (field_data.value_status[active_index] == value::status::uninitialized) {
                        std::size_t layer_index = i + j*this->nx;
                        if (toplayer.value_status[layer_index] == value::status::deck_value) {
                            field_data.data[active_index] = toplayer.data[layer_index];
                            field_data.value_status[active_index] = value::status::valid_default;
                        }
                    }
                    active_index += 1;
                }
            }
        }
    }
}

template <>
bool FieldProps::supported<double>(const std::string& keyword)
{
    if (Fieldprops::keywords::GRID::double_keywords.count(keyword) != 0) {
        return true;
    }

    if (Fieldprops::keywords::EDIT::double_keywords.count(keyword) != 0) {
        return true;
    }

    if (Fieldprops::keywords::PROPS::double_keywords.count(keyword) != 0) {
        return true;
    }

    if (Fieldprops::keywords::PROPS::satfunc.count(keyword) != 0) {
        return true;
    }

    if (Fieldprops::keywords::SOLUTION::double_keywords.count(keyword) != 0 ||
        Fieldprops::keywords::SOLUTION::composition_keywords.count(keyword) != 0) {
        return true;
    }

    return false;
}

template <>
bool FieldProps::supported<int>(const std::string& keyword)
{
    if (Fieldprops::keywords::REGIONS::int_keywords.count(keyword) != 0) {
        return true;
    }

    if (Fieldprops::keywords::GRID::int_keywords.count(keyword) != 0) {
        return true;
    }

    if (Fieldprops::keywords::SCHEDULE::int_keywords.count(keyword) != 0) {
        return true;
    }

    return Fieldprops::keywords::isFipxxx(keyword);
}

std::pair<std::vector<Box::cell_index>,bool>
FieldProps::region_index(const std::string& region_name, const int region_value)
{
    std::vector<Box::cell_index> index_list;
    bool all_active = true;
    const auto& region = this->init_get<int>(region_name);
    if (!region.valid()) {
        throw std::invalid_argument("Trying to work with invalid region: " + region_name);
    }

    std::size_t active_index = 0;
    for (std::size_t g = 0; g < this->m_actnum.size(); ++g) {
        if (this->m_actnum[g] != 0) {
            if (region.data[active_index] == region_value) {
                index_list.emplace_back(g, active_index, g);
            }

            active_index += 1;
        }else{
            all_active=false;
        }
    }

    return {index_list, all_active};
}

std::string FieldProps::region_name(const DeckItem& region_item) const
{
    return region_item.defaultApplied(0)
        ? this->m_default_region
        : make_region_name(region_item.get<std::string>(0));
}

template <>
bool FieldProps::has<double>(const std::string& keyword_name) const
{
    const auto keyword = Fieldprops::keywords::get_keyword_from_alias(keyword_name);

    return this->double_data.find(keyword) != this->double_data.end();
}

template <>
bool FieldProps::has<int>(const std::string& keyword) const
{
    const auto& kw = Fieldprops::keywords::isFipxxx(keyword)
        ? this->canonical_fipreg_name(keyword)
        : keyword;

    return this->int_data.find(kw) != this->int_data.end();
}

void FieldProps::apply_multipliers()
{
    // We need to manually search for PORV in the map here instead of using
    // the get method. The latter one will compute PORV from the cell volume
    // using MULTPV, NTG, and PORO. Our intend here in the EDIT section is
    // is to multiply exsisting MULTPV with new new ones and consistently
    // change PORV as well. If PORV has been created before this is as easy
    // as multiplying it and MULTPV with the additional MULTPV. In the other
    // case we do not create PORV at all but just change MULTPV as PORV will
    // be correctly computed from it.
    // Hence we need to know whether PORV is already there
    const bool hasPorvBefore =
        (this->double_data.find(ParserKeywords::PORV::keywordName) ==
         this->double_data.end());
    static const auto prefix = getMultiplierPrefix();

    for(const auto& [mult_keyword, kw_info]: multiplier_kw_infos_)
    {
        const std::string keyword = mult_keyword.substr(prefix.size());
        auto mult_iter = this->double_data.find(mult_keyword);
        assert(mult_iter != this->double_data.end());
        auto iter = this->double_data.find(keyword);
        if (iter == this->double_data.end()) {
            iter = this->double_data
                .emplace(std::piecewise_construct,
                         std::forward_as_tuple(keyword),
                         std::forward_as_tuple(kw_info, this->active_size, kw_info.global ? this->global_size : 0))
                .first;
        }

        std::transform(iter->second.data.begin(), iter->second.data.end(),
                       mult_iter->second.data.begin(), iter->second.data.begin(),
                       std::multiplies<>());

        // If data is global, then we also need to set the global_data. I think they should be the same at this stage, though!
        if (kw_info.global)
        {
            assert(mult_iter->second.global_data.has_value());
            assert(iter->second.global_data.has_value());
            std::transform(iter->second.global_data->begin(),
                           iter->second.global_data->end(),
                           mult_iter->second.global_data->begin(),
                           iter->second.global_data->begin(),
                           std::multiplies<>());
        }
        // If this is MULTPV we also need to apply the additional multiplier to PORV if that was initialized already.
        // Note that the check for PORV is essential as otherwise the value constructed durig init_get will already apply
        // current MULTPV.
        if (keyword == ParserKeywords::MULTPV::keywordName && !hasPorvBefore) {
            auto& porv = this->init_get<double>(ParserKeywords::PORV::keywordName);
            auto& porv_data = porv.data;
            std::transform(porv_data.begin(), porv_data.end(), mult_iter->second.data.begin(), porv_data.begin(), std::multiplies<>());
        }
        this->double_data.erase(mult_iter);
    }
    multiplier_kw_infos_.clear();
}

// The ACTNUM and PORV keywords are special cased with quite extensive
// postprocessing, and should be accessed through the special ::porv() and
// ::actnum() methods instead of the general ::get<T>( ) method. These two
// keywords are also hidden from the keys<T>() vectors.
//
// If there are TRAN? fields in the container, these will be transferred
// even if they are not fully defined. These TRAN? fields will ultimately be
// combined with the TRAN? values calculated from the simulator, and it does
// therefore not make sense to require fully defined fields.

template <>
std::vector<std::string> FieldProps::keys<double>() const
{
    std::vector<std::string> klist;

    for (const auto& [key, field] : this->double_data) {
        if (key.rfind("TRAN", 0) == 0) {
            klist.push_back(key);
            continue;
        }

        if (field.valid() && (key != "PORV")) {
            klist.push_back(key);
        }
    }

    return klist;
}

template <>
std::vector<std::string> FieldProps::keys<int>() const
{
    std::vector<std::string> klist;

    for (const auto& [key, field] : this->int_data) {
        if (field.valid() && (key != "ACTNUM")) {
            klist.push_back(key);
        }
    }

    return klist;
}

template <>
void FieldProps::erase<int>(const std::string& keyword)
{
    this->int_data.erase(keyword);
}

template <>
void FieldProps::erase<double>(const std::string& keyword)
{
    this->double_data.erase(keyword);
}

template <>
std::vector<int> FieldProps::extract<int>(const std::string& keyword)
{
    auto field_iter = this->int_data.find(keyword);

    auto field = std::move(field_iter->second);
    std::vector<int> data = std::move(field.data);

    this->int_data.erase(field_iter);

    return data;
}

template <>
std::vector<double> FieldProps::extract<double>(const std::string& keyword)
{
    auto field_iter = this->double_data.find(keyword);

    auto field = std::move(field_iter->second);
    std::vector<double> data = std::move(field.data);

    this->double_data.erase(field_iter);

    return data;
}

double FieldProps::getSIValue(const std::string& keyword, const double raw_value) const
{
    if (this->tran.find(keyword) != this->tran.end()) {
        return this->unit_system.to_si(UnitSystem::measure::transmissibility, raw_value);
    }

    if (const auto& kw_info = Fieldprops::keywords::global_kw_info<double>(keyword);
        kw_info.unit)
    {
        const auto& dim = this->unit_system.parse(*kw_info.unit);

        return dim.convertRawToSi(raw_value);
    }

    return raw_value;
}

double FieldProps::getSIValue(const ScalarOperation op,
                              const std::string& keyword,
                              const double raw_value) const
{
    return (op == ScalarOperation::MUL)
        ? raw_value
        : this->getSIValue(keyword, raw_value);
}

void FieldProps::handle_int_keyword(const Fieldprops::keywords::keyword_info<int>& kw_info,
                                    const DeckKeyword& keyword,
                                    const Box& box)
{
    auto& field_data = this->init_get<int>(keyword.name());

    const auto& deck_data = keyword.getIntData();
    const auto& deck_status = keyword.getValueStatus();

    assign_deck(kw_info, keyword, field_data, deck_data, deck_status, box);
}

void FieldProps::handle_double_keyword(const Section section,
                                       const Fieldprops::keywords::keyword_info<double>& kw_info,
                                       const DeckKeyword& keyword,
                                       const std::string& keyword_name,
                                       const Box& box)
{
    // if second paramter is true then this will not be the actual keyword
    // but one prefixed with __MULT__ that will be used to construct the
    // multiplier for later application to the actual keyword.
    auto& field_data = this->init_get<double>
        (keyword_name, kw_info, (section == Section::EDIT) && kw_info.multiplier);

    const auto& deck_data = keyword.getSIDoubleData();
    const auto& deck_status = keyword.getValueStatus();

    if ((section == Section::SCHEDULE) && kw_info.multiplier) {
        // Apply all multipliers cumulatively
        multiply_deck(kw_info, keyword, field_data, deck_data, deck_status, box);
    }
    else {
        // Apply only latest multiplier (overwrite these previous one)
        assign_deck(kw_info, keyword, field_data, deck_data, deck_status, box);
    }

    if (section == Section::GRID) {
        if (field_data.valid()) {
            return;
        }

        if (kw_info.top) {
            this->distribute_toplayer(field_data, deck_data, box);
        }
    }
}

void FieldProps::handle_double_keyword(const Section section,
                                       const Fieldprops::keywords::keyword_info<double>& kw_info,
                                       const DeckKeyword& keyword,
                                       const Box& box)
{
    this->handle_double_keyword(section, kw_info, keyword, keyword.name(), box);
}

double FieldProps::get_alpha(const std::string& func_name,
                             const std::string& target_array,
                             const double       raw_alpha)
{
    return ((func_name == "ADDX") ||
            (func_name == "MAXLIM") ||
            (func_name == "MINLIM"))
        ? this->getSIValue(target_array, raw_alpha)
        : raw_alpha;
}

double FieldProps::get_beta(const std::string& func_name,
                            const std::string& target_array,
                            const double       raw_beta)
{
    return (func_name == "MULTA")
        ? this->getSIValue(target_array, raw_beta)
        : raw_beta;
}

template <typename T>
void FieldProps::operate(const DeckRecord&                   record,
                         Fieldprops::FieldData<T>&           target_data,
                         const Fieldprops::FieldData<T>&     src_data,
                         const std::vector<Box::cell_index>& index_list,
                         const bool                          global)
{
    const auto target_array = record.getItem("TARGET_ARRAY").getTrimmedString(0);
    if (this->tran.find(target_array) != this->tran.end()) {
        throw std::logic_error {
            "The OPERATE keyword cannot be used for "
            "manipulations of TRANX, TRANY or TRANZ"
        };
    }

    const auto func_name    = record.getItem("OPERATION").getTrimmedString(0);
    const auto check_target = (func_name == "MULTIPLY") || (func_name == "POLY");

    const auto alpha = this->get_alpha(func_name, target_array, record.getItem("PARAM1").get<double>(0));
    const auto beta  = this->get_beta(func_name, target_array, record.getItem("PARAM2").get<double>(0));
    const auto func  = Operate::get(func_name, alpha, beta);

    auto& to_data = global? *target_data.global_data : target_data.data;
    auto& to_status = global? *target_data.global_value_status : target_data.value_status;
    const auto& from_data = global? *src_data.global_data : src_data.data;
    auto& from_status = global? *src_data.global_value_status : src_data.value_status;

    for (const auto& cell_index : index_list) {
        // This is the global index if global is true and global storage is used.
        const auto ix = cell_index.active_index;

        if (value::has_value(from_status[ix])) {
            if (!check_target || value::has_value(to_status[ix])) {
                to_data[ix] = func(to_data[ix], from_data[ix]);
                to_status[ix] = from_status[ix];
            }
            else {
                throw std::invalid_argument {
                    "Tried to use unset property value "
                    "in OPERATE/OPERATER keyword"
                };
            }
        }
        else {
            throw std::invalid_argument {
                "Tried to use unset property value in "
                "OPERATE/OPERATER keyword"
            };
        }
    }
}

void FieldProps::handle_operateR(const DeckKeyword& keyword)
{
    // Special case handling for OPERATE*R*.  General keyword structure is
    //
    //   OPERATER
    //     ResArray  RegionID  Operation  SrcArray  a  b  RegionSet /
    //   -- ...
    //   /
    //
    // which applies the 'Operation' to all SrcArray elements within region
    // 'RegionID' of the specified 'RegionSet'.  The operation typically
    // incorporates at least one of the parameters 'a' and 'b'.  Resulting
    // values overwrite the corresponding elements of the result/target
    // array (ResArray).

    for (const auto& record : keyword) {
        const auto target_kw = Fieldprops::keywords::
            get_keyword_from_alias(record.getItem(0).getTrimmedString(0));

        if (! FieldProps::supported<double>(target_kw)) {
            continue;
        }

        if (this->tran.find(target_kw) != this->tran.end()) {
            throw std::logic_error {
                "The region operations cannot be used for "
                "manipulations of TRANX, TRANY or TRANZ"
            };
        }

        const int region_value = record.getItem("REGION_NUMBER").get<int>(0);

        auto& field_data = this->init_get<double>(target_kw);

        // For the OPERATER keyword we fetch the region name from the deck
        // record with no extra hoops.
        const auto reg_name = record.getItem("REGION_NAME").getTrimmedString(0);
        const auto src_kw = record.getItem("ARRAY_PARAMETER").getTrimmedString(0);

        const auto& [index_list, all_active] = this->region_index(reg_name, region_value);
        if (index_list.empty()) {
            log_empty_region(keyword, reg_name, region_value, src_kw);
            continue;
        }

        const auto& src_data = this->init_get<double>(src_kw);
        FieldProps::operate(record, field_data, src_data, index_list);

        // Supporting region operations on global storage arrays would
        // require global storage for the *NUM region set arrays (i.e.,
        // FLUXNUM, MULTNUM, OPERNUM).  In order to avoid global storage
        // for a significant portion of the property arrays we have
        // decided, as a project policy, to not support such operations
        // at this time.  There is, however, no technical reason for why
        // the current implementation could not be extended to support
        // region operations on global storage arrays.
        // For now regions do not 100% support global storage.
        // Make sure that the global storage at least reflects the local one.
        if (field_data.global_data) {
            if (!all_active) {
                const auto message =
                    fmt::format(R"(Region operation on 3D field {} with global storage will not update inactive cells.
Note that this might cause problems for PINCH option 4 or 5 being ALL.)", target_kw);

                OpmLog::warning(Log::fileMessage(keyword.location(), message));
            }
            update_global_from_local(field_data, index_list);
        }
    }
}

void FieldProps::handle_region_operation(const DeckKeyword& keyword)
{
    if (keyword.name() == ParserKeywords::OPERATER::keywordName) {
        // Special case handling for OPERATER.
        this->handle_operateR(keyword);
        return;
    }

    // If we get here, we're processing ADDREG, EQUALREG, or MULTIREG.
    // General keyword structure is as follows
    //
    //   OpKeywd
    //     Array Scalar RegionID RegionSet /
    //   -- ...
    //   /
    //
    // which applies the operation (ADDREG = addition of 'Scalar', EQUALREG
    // = assignment of 'Scalar', MULTIREG = multiplication by 'Scalar') to
    // each Array element for which the corresponding cell is in RegionID of
    // the specified RegionSet.

    const auto operation = fromString(keyword.name());

    for (const auto& record : keyword) {
        const auto target_kw = Fieldprops::keywords::
            get_keyword_from_alias(record.getItem(0).getTrimmedString(0));

        if (this->tran.find(target_kw) != this->tran.end()) {
            throw std::logic_error {
                "The region operations cannot be used for "
                "manipulations of TRANX, TRANY or TRANZ"
            };
        }

        const int region_value = record.getItem("REGION_NUMBER").get<int>(0);

        if (FieldProps::supported<double>(target_kw)) {
            auto& field_data = this->init_get<double>(target_kw);

            const auto reg_name = this->region_name(record.getItem("REGION_NAME"));
            const auto& [index_list, all_active] = this->region_index(reg_name, region_value);
            if (index_list.empty()) {
                log_empty_region(keyword, reg_name, region_value, target_kw);
                continue;
            }

            const auto scalar_value =
                this->getSIValue(operation, target_kw, record.getItem(1).get<double>(0));

            apply(operation, keyword.location(), target_kw,
                  field_data.data, field_data.value_status,
                  scalar_value, index_list);

            // Supporting region operations on global storage arrays would
            // require global storage for the *NUM region set arrays (i.e.,
            // FLUXNUM, MULTNUM, OPERNUM).  In order to avoid global storage
            // for a significant portion of the property arrays we have
            // decided, as a project policy, to not support such operations
            // at this time.  There is, however, no technical reason for why
            // the current implementation could not be extended to support
            // region operations on global storage arrays..
            // Make sure that the global storage at least reflects the local one.
            if (field_data.global_data) {
                if (!all_active) {
                    const auto message =
                        fmt::format(R"(Region operation on 3D field {} with global storage will not update inactive cells.
Note that this might cause problems for PINCH option 4 or 5 being ALL.)", target_kw);

                    OpmLog::warning(Log::fileMessage(keyword.location(), message));
                }
                update_global_from_local(field_data, index_list);
            }
            continue;
        }

        if (FieldProps::supported<int>(target_kw)) {
            continue;
        }
    }
}

void FieldProps::handle_OPERATE(const DeckKeyword& keyword, Box box)
{
    // Implementation of the OPERATE keyword.  General keyword structure is
    //
    //   OPERATE
    //     ResArray  Box  Operation  SrcArray  a  b /
    //   -- ...
    //   /
    //
    // which applies the 'Operation' to all SrcArray elements within the
    // 'Box', typically incorporating the parameters 'a' and 'b'.  Resulting
    // values overwrite the corresponding elements of the result/target
    // array (ResArray).

    for (const auto& record : keyword) {
        box.update(record);

        const auto target_kw = Fieldprops::keywords::
            get_keyword_from_alias(record.getItem(0).getTrimmedString(0));

        auto& field_data = this->init_get<double>(target_kw);

        const auto src_kw = record.getItem("ARRAY").getTrimmedString(0);
        const auto& src_data = this->init_get<double>(src_kw);

        FieldProps::operate(record, field_data, src_data, box.index_list());

        if (field_data.global_data)
        {
            if (!src_data.global_data) {
                throw std::logic_error {
                    "The OPERATE and OPERATER keywords are only "
                    "supported between keywords with same storage"
                };
            }

            FieldProps::operate(record, field_data, src_data, box.global_index_list(), true);
        }
    }
}

void FieldProps::handle_operation(const Section      section,
                                  const DeckKeyword& keyword,
                                  Box                box,
                                  const bool         onlyACTNUM)
{
    // Keyword handler for ADD, EQUALS, MAXVALUE, MINVALUE, and MULTIPLY.
    //
    // General keyword structure:
    //
    //   KWNAME -- e.g., ADD or EQUALS
    //     Array Scalar Box /
    //     Array Scalar Box /
    //   -- ...
    //   /
    //
    // For example
    //
    //    ADD
    //      PERMX  12.3  1 10   2 4   5 6 /
    //    /
    //
    // which adds 12.3 mD to the PERMX array in all cells within the box
    // {1..10, 2..4, 5..6}.  The array being operated on must already exist
    // and be well defined for all elements in the box unless the operation
    // is "EQUALS" or we're processing one of the TRAN* arrays or PORV in
    // the EDIT section.  Final TRAN* array processing is deferred to a
    // later stage and at this point we're only collecting descriptors of
    // what operations to apply when we get to that stage.

    const auto mustExist = keyword.name() != ParserKeywords::EQUALS::keywordName;
    const auto editSect  = section == Section::EDIT;
    const auto operation = fromString(keyword.name());

    std::unordered_map<std::string, std::string> tran_fields;

    for (const auto& record : keyword) {
        const auto target_kw = Fieldprops::keywords::
            get_keyword_from_alias(record.getItem(0).getTrimmedString(0));

        if (onlyACTNUM && target_kw != "ACTNUM") {
            continue;
        }

        box.update(record);

        if (auto tran_iter = this->tran.find(target_kw);
            FieldProps::supported<double>(target_kw) ||
            (tran_iter != this->tran.end()))
        {
            auto kw_info = Fieldprops::keywords::global_kw_info<double>
                (target_kw, /* allow_unsupported = */ tran_iter != this->tran.end());

            auto unique_name = target_kw;

            // Check if the target array is one of TRANX, TRANY, or TRANZ.
            if (tran_iter != this->tran.end()) {
                // The transmissibility calculations are applied to one
                // "work" array per direction and per array operation (i.e.,
                // keyword.name()).  If we have not already seen this
                // transmissibility direction while processing this keyword,
                // then we register a new transmissibility calculator
                // operation.
                auto tran_field_iter = tran_fields.find(target_kw);
                if (tran_field_iter == tran_fields.end()) {
                    unique_name = tran_iter->second.next_name();

                    tran_fields.emplace(target_kw, unique_name);
                    tran_iter->second.add_action(operation, unique_name);

                    kw_info = tran_iter->second.make_kw_info(operation, target_kw);
                }
                else {
                    unique_name = tran_field_iter->second;
                }
            }
            else if (mustExist && !kw_info.multiplier &&
                     !(editSect && (unique_name == ParserKeywords::PORV::keywordName)) &&
                     (this->double_data.find(unique_name) == this->double_data.end()))
            {
                // Note exceptions for the MULT* arrays (i.e., MULT[XYZ] and
                // MULT[XYZ]-).  We always support operating on defaulted
                // array values (i.e, all elements equal to one) in the case
                // of those arrays, even if the operation is not assignment.

                throw OpmInputError {
                    fmt::format("Target array {} must already "
                                "exist when operated upon in {}.",
                                target_kw, keyword.name()),
                    keyword.location()
                };
            }

            const auto scalar_value = this->
                getSIValue(operation, target_kw, record.getItem(1).get<double>(0));

            auto& field_data = this->init_get<double>
                (unique_name, kw_info, /* multiplier_in_edit =*/ editSect && kw_info.multiplier);

            apply(operation, keyword.location(), target_kw,
                  field_data.data, field_data.value_status,
                  scalar_value, box.index_list());

            if (field_data.global_data) {
                apply(operation, keyword.location(), target_kw,
                      *field_data.global_data,
                      *field_data.global_value_status,
                      scalar_value, box.global_index_list());
            }

            continue;
        }

        if (FieldProps::supported<int>(target_kw)) {
            if (mustExist && (this->int_data.find(target_kw) == this->int_data.end())) {
                throw OpmInputError {
                    fmt::format("Target array {} must already "
                                "exist when operated upon in {}.",
                                target_kw, keyword.name()),
                    keyword.location()
                };
            }

            const auto scalar_value = static_cast<int>(record.getItem(1).get<double>(0));

            auto& field_data = this->init_get<int>(target_kw);

            apply(operation, keyword.location(), target_kw,
                  field_data.data,
                  field_data.value_status,
                  scalar_value, box.index_list());

            continue;
        }

        throw OpmInputError {
            fmt::format("Target array {} is not supported in the "
                        "{} operation", target_kw, keyword.name()),
            keyword.location()
        };
    }
}

void FieldProps::handle_COPY(const DeckKeyword& keyword,
                             Box                box,
                             const bool         isRegionOperation)
{
    auto arrayName = [](const auto& item)
    {
        return Fieldprops::keywords::get_keyword_from_alias(item.getTrimmedString(0));
    };

    for (const auto& record : keyword) {
        const auto src_kw    = arrayName(record.getItem(0));
        const auto target_kw = arrayName(record.getItem(1));

        std::vector<Box::cell_index> index_list;
        auto srcDescr = std::string {};

        if (isRegionOperation) {
            using Kw = ParserKeywords::COPYREG;
            const auto  regionId   = record.getItem<Kw::REGION_NUMBER>().get<int>(0);
            const auto& regionName = this->region_name(record.getItem<Kw::REGION_NAME>());

            index_list = this->region_index(regionName, regionId).first;
            srcDescr = fmt::format("{} in region {} of region set {}",
                                   src_kw, regionId, regionName);
        }
        else {
            box.update(record);
            index_list = box.index_list();

            srcDescr = fmt::format("{} in BOX ({}-{}, {}-{}, {}-{})",
                                   src_kw,
                                   box.I1() + 1, box.I2() + 1,
                                   box.J1() + 1, box.J2() + 1,
                                   box.K1() + 1, box.K2() + 1);
        }

        if (FieldProps::supported<double>(src_kw)) {
            const auto& src_data = this->try_get<double>(src_kw, TryGetFlags::MustExist);
            src_data.verify_status(keyword.location(), "Source array", "COPY");

            auto& target_data = this->init_get<double>(target_kw);
            target_data.checkInitialisedCopy(src_data.field_data(), index_list,
                                             srcDescr, target_kw,
                                             keyword.location());
            if (target_data.global_data && !isRegionOperation) {
                if (!src_data.field_data().global_data) {
                    throw std::logic_error {
                        fmt::format
                        (R"(The copying is only supported between keywords with same storage.
 (COPY {} {})", src_kw, target_kw)
                    };
                }
                target_data.checkInitialisedCopy(src_data.field_data(), box.global_index_list(),
                                                 srcDescr, target_kw,
                                                 keyword.location(),
                                                 true);
            }
            continue;
        }

        if (FieldProps::supported<int>(src_kw)) {
            const auto& src_data = this->try_get<int>(src_kw, TryGetFlags::MustExist);
            src_data.verify_status(keyword.location(), "Source array", "COPY");

            auto& target_data = this->init_get<int>(target_kw);
            target_data.checkInitialisedCopy(src_data.field_data(), index_list,
                                             srcDescr, target_kw,
                                             keyword.location());
            continue;
        }
    }
}

void FieldProps::handle_keyword(const Section      section,
                                const DeckKeyword& keyword,
                                Box&               box,
                                const bool         onlyACTNUM)
{
    const auto& name = keyword.name();

    if (Fieldprops::keywords::oper_keywords.count(name) == 1) {
        this->handle_operation(section, keyword, box, onlyACTNUM);
    }

    else if (name == ParserKeywords::OPERATE::keywordName) {
        this->handle_OPERATE(keyword, box);
    }

    else if (Fieldprops::keywords::region_oper_keywords.count(name) == 1) {
        this->handle_region_operation(keyword);
    }

    else if (Fieldprops::keywords::box_keywords.count(name) == 1) {
        handle_box_keyword(keyword, box);
    }

    else if ((name == ParserKeywords::COPY::keywordName) ||
             (name == ParserKeywords::COPYREG::keywordName))
    {
        handle_COPY(keyword, box, name == ParserKeywords::COPYREG::keywordName);
    }
}

// ---------------------------------------------------------------------------

void FieldProps::init_tempi(Fieldprops::FieldData<double>& tempi)
{
    if (this->tables.hasTables("RTEMPVD")) {
        const auto& eqlnum = this->get<int>("EQLNUM");
        const auto& rtempvd = this->tables.getRtempvdTables();
        std::vector< double > tempi_values( this->active_size, 0 );

        for (size_t active_index = 0; active_index < this->active_size; active_index++) {
            const auto& table = rtempvd.getTable<RtempvdTable>(eqlnum[active_index] - 1);
            double depth = this->cell_depth[active_index];
            tempi_values[active_index] = table.evaluate("Temperature", depth);
        }

        tempi.default_update(tempi_values);
    } else
        tempi.default_assign(this->tables.rtemp());
}

void FieldProps::init_porv(Fieldprops::FieldData<double>& porv)
{
    auto& porv_data = porv.data;
    auto& porv_status = porv.value_status;

    const auto& poro = this->init_get<double>("PORO");
    const auto& poro_status = poro.value_status;
    const auto& poro_data = poro.data;

    for (std::size_t active_index = 0; active_index < this->active_size; active_index++) {
        if (value::has_value(poro_status[active_index])) {
            porv_data[active_index] = this->cell_volume[active_index] * poro_data[active_index];
            porv_status[active_index] = value::status::valid_default;
        }
    }

    if (this->has<double>("NTG")) {
        const auto& ntg = this->get<double>("NTG");
        for (std::size_t active_index = 0; active_index < this->active_size; active_index++)
            porv_data[active_index] *= ntg[active_index];
    }

    if (this->has<double>("MULTPV")) {
        const auto& multpv = this->get<double>("MULTPV");
        std::transform(porv_data.begin(), porv_data.end(), multpv.begin(), porv_data.begin(), std::multiplies<>());
    }

    for (const auto& mregp: this->multregp) {
        const auto index_list = this->region_index(mregp.region_name, mregp.region_value).first;
        for (const auto& cell_index : index_list)
            porv_data[cell_index.active_index] *= mregp.multiplier;
    }
}

std::string FieldProps::canonical_fipreg_name(const std::string& fipreg)
{
    constexpr auto num_unique_chars = std::string::size_type{6};
    const auto shortname = fipreg.substr(0, num_unique_chars);

    auto canonicalPos = std::unordered_map<std::string, std::string>::const_iterator{};
    auto inserted = false;

    std::tie(canonicalPos, inserted) = this->fipreg_shortname_translation
        .emplace(shortname, fipreg);

    if (inserted || (fipreg.size() <= num_unique_chars)) {
        // New FIP keyword (inserted == true) or we're looking up the
        // canonical name of an existing FIP array based on a unique prefix
        // string, such as "FIPUNIT" from "FIPUNI" (size <= ...).
        return canonicalPos->second;
    }

    // New FIP keyword with the same unique prefix as an existing FIP
    // keyword.  Override the translation table entry for this prefix,
    // because "last entry wins".
    canonicalPos = this->fipreg_shortname_translation
        .insert_or_assign(canonicalPos, shortname, fipreg);

    return canonicalPos->second;
}

const std::string&
FieldProps::canonical_fipreg_name(const std::string& fipreg) const
{
    auto canonicalPos = this->fipreg_shortname_translation.find(fipreg.substr(0, 6));

    return (canonicalPos != this->fipreg_shortname_translation.end())
        ? canonicalPos->second
        : fipreg;
}

// Generate a combined ACTNUM property array from three distinct data
// sources.
//
//    1. The property array stored internally in this FieldProps object.
//
//    2. Direct ACTNUM array operations of the form
//
//       EQUALS
//           ACTNUM 0 1 10 1 10 1 3 /
//       /
//
//    3. Cells with PORV == 0 will get ACTNUM = 0.
//
// Note that steps 2 and 3 generally forms an ACTNUM property which differs
// from the ACTNUM property stored internally in the FieldProps instance.
std::vector<int> FieldProps::actnum()
{
    auto actnum = this->m_actnum;

    // Avoid de-activating all cells if PORO has not yet been read (typically in tests)
    if (!this->has<double>("PORO")) {
        return actnum;
    }

    const auto& deck_actnum = this->init_get<int>("ACTNUM");

    std::vector<int> global_map(this->active_size);
    {
        std::size_t active_index = 0;
        for (std::size_t g = 0; g < this->global_size; g++) {
            if (this->m_actnum[g]) {
                global_map[active_index] = g;
                active_index++;
            }
        }
    }


    const auto& porv = this->init_get<double>("PORV");
    const auto& porv_data = porv.data;
    for (std::size_t active_index = 0; active_index < this->active_size; active_index++) {
        auto global_index = global_map[active_index];
        actnum[global_index] = deck_actnum.data[active_index];
        if (porv_data[active_index] == 0)
            actnum[global_index] = 0;
    }
    return actnum;
}

const std::vector<int>& FieldProps::actnumRaw() const
{
    return m_actnum;
}

void FieldProps::processMULTREGP(const Deck& deck)
{
    using Kw = ParserKeywords::MULTREGP;

    for (const auto& keyword : deck[Kw::keywordName]) {
        for (const auto& record : keyword) {
            const int region_value = record.getItem<Kw::REGION>().get<int>(0);
            if (region_value <= 0) {
                continue;
            }

            const auto reg_name =
                make_region_name(record.getItem<Kw::REGION_TYPE>().getTrimmedString(0));

            // Can't use getSIDouble(0) here as there's no defined dimension
            // for the multiplier item in keyword MULTREGP.
            const auto multiplier = record.getItem<Kw::MULTIPLIER>().get<double>(0);

            auto iter = std::find_if(this->multregp.begin(), this->multregp.end(),
                                     [region_value](const MultregpRecord& mregp)
                                     { return mregp.region_value == region_value; });

            // There is some weirdness if the same region value is entered
            // in several records, then only the last applies.
            if (iter != this->multregp.end()) {
                iter->region_name = reg_name;
                iter->multiplier = multiplier;
            }
            else {
                this->multregp.emplace_back(region_value, multiplier, reg_name);
            }
        }
    }
}

void FieldProps::scanGRIDSection(const GRIDSection& grid_section)
{
    auto box = makeGlobalGridBox(this->grid_ptr);

    for (const auto& keyword : grid_section) {
        if (auto kwPos = Fieldprops::keywords::GRID::double_keywords.find(keyword.name());
            kwPos != Fieldprops::keywords::GRID::double_keywords.end())
        {
            this->handle_double_keyword(Section::GRID, kwPos->second, keyword, box);
            continue;
        }

        if (auto kwPos = Fieldprops::keywords::GRID::int_keywords.find(keyword.name());
            kwPos != Fieldprops::keywords::GRID::int_keywords.end())
        {
            this->handle_int_keyword(kwPos->second, keyword, box);
            continue;
        }

        this->handle_keyword(Section::GRID, keyword, box);
    }
}

void FieldProps::scanGRIDSectionOnlyACTNUM(const GRIDSection& grid_section)
{
    Box box(*this->grid_ptr, [](const std::size_t) { return true; }, [](const std::size_t i) { return i; });

    for (const auto& keyword : grid_section) {
        const std::string& name = keyword.name();

        if (name == "ACTNUM") {
            this->handle_int_keyword(Fieldprops::keywords::GRID::int_keywords.at(name), keyword, box);
        }
        else if ((name == "EQUALS") || (Fieldprops::keywords::box_keywords.count(name) == 1)) {
            this->handle_keyword(Section::GRID, keyword, box, true);
        }
    }

    if (auto iter = this->int_data.find("ACTNUM");
        iter == this->int_data.end())
    {
        m_actnum.assign(this->grid_ptr->getCartesianSize(), 1);
    }
    else {
        m_actnum = iter->second.data;
    }
}

void FieldProps::scanEDITSection(const EDITSection& edit_section)
{
    auto box = makeGlobalGridBox(this->grid_ptr);

    for (const auto& keyword : edit_section) {
        const std::string& name = keyword.name();

        if (auto tran_iter = this->tran.find(name);
            tran_iter!= this->tran.end())
        {
            auto& tran_calc = tran_iter->second;
            auto unique_name = tran_calc.next_name();
            this->handle_double_keyword(Section::EDIT, {}, keyword, unique_name, box);
            tran_calc.add_action( Fieldprops::ScalarOperation::EQUAL, unique_name );
            continue;
        }

        if (auto kwPos = Fieldprops::keywords::EDIT::double_keywords.find(name);
            kwPos != Fieldprops::keywords::EDIT::double_keywords.end())
        {
            this->handle_double_keyword(Section::EDIT, kwPos->second, keyword, box);
            continue;
        }

        if (auto kwPos = Fieldprops::keywords::EDIT::int_keywords.find(name);
            kwPos != Fieldprops::keywords::EDIT::int_keywords.end())
        {
            this->handle_int_keyword(kwPos->second, keyword, box);
            continue;
        }

        this->handle_keyword(Section::EDIT, keyword, box);
    }
    // Multiplier will not have been applied yet to prevent EQUALS MULT* from overwriting values
    // and to only honor the last MULT* occurrence
    // apply recorded multipliers of section to existing ones
    apply_multipliers();
}

void FieldProps::init_satfunc(const std::string& keyword,
                              Fieldprops::FieldData<double>& satfunc)
{
    if (!this->m_rtep.has_value())
        this->m_rtep = satfunc::getRawTableEndpoints(this->tables, this->m_phases,
                                                     this->m_satfuncctrl.minimumRelpermMobilityThreshold());

    const auto& endnum = this->get<int>("ENDNUM");
    const auto& satreg = (keyword[0] == 'I')
        ? this->get<int>("IMBNUM")
        : this->get<int>("SATNUM");

    satfunc.default_update(satfunc::init(keyword, this->tables, this->m_phases, this->m_rtep.value(), this->cell_depth, satreg, endnum));
}

void FieldProps::scanPROPSSection(const PROPSSection& props_section)
{
    auto box = makeGlobalGridBox(this->grid_ptr);

    for (const auto& keyword : props_section) {
        const std::string& name = keyword.name();
        if (Fieldprops::keywords::PROPS::satfunc.count(name) == 1) {
            Fieldprops::keywords::keyword_info<double> sat_info{};
            this->handle_double_keyword(Section::PROPS, sat_info, keyword, box);
            continue;
        }

        if (auto kwPos = Fieldprops::keywords::PROPS::double_keywords.find(name);
            kwPos != Fieldprops::keywords::PROPS::double_keywords.end())
        {
            this->handle_double_keyword(Section::PROPS, kwPos->second, keyword, box);
            continue;
        }

        if (auto kwPos = Fieldprops::keywords::PROPS::int_keywords.find(name);
            kwPos != Fieldprops::keywords::PROPS::int_keywords.end())
        {
            this->handle_int_keyword(kwPos->second, keyword, box);
            continue;
        }

        this->handle_keyword(Section::PROPS, keyword, box);
    }
}

void FieldProps::scanREGIONSSection(const REGIONSSection& regions_section)
{
    auto box = makeGlobalGridBox(this->grid_ptr);

    for (const auto& keyword : regions_section) {
        const std::string& name = keyword.name();

        if (auto kwPos = Fieldprops::keywords::REGIONS::int_keywords.find(name);
            kwPos != Fieldprops::keywords::REGIONS::int_keywords.end())
        {
            this->handle_int_keyword(kwPos->second, keyword, box);
            continue;
        }

        if (Fieldprops::keywords::isFipxxx(name)) {
            auto kw_info = Fieldprops::keywords::keyword_info<int>{};
            kw_info.init(1);
            this->handle_int_keyword(kw_info, keyword, box);
            continue;
        }

        this->handle_keyword(Section::REGIONS, keyword, box);
    }
}

void FieldProps::scanSOLUTIONSection(const SOLUTIONSection& solution_section,
                                     const std::size_t      ncomps)
{
    auto box = makeGlobalGridBox(this->grid_ptr);

    for (const auto& keyword : solution_section) {
        const std::string& name = keyword.name();

        if (auto kwPos = Fieldprops::keywords::SOLUTION::double_keywords.find(name);
            kwPos != Fieldprops::keywords::SOLUTION::double_keywords.end())
        {
            this->handle_double_keyword(Section::SOLUTION, kwPos->second, keyword, box);
            continue;
        }

        if (auto kwPos = Fieldprops::keywords::SOLUTION::composition_keywords.find(name);
            kwPos != Fieldprops::keywords::SOLUTION::composition_keywords.end())
        {
            if (ncomps < 1) {
                throw std::logic_error {
                    fmt::format("With compostional keyword {} defined in SOLUTION, while the DATA file "
                                "does not appear to be a compostional case.", name)
                };
            }

            // TODO: maybe we should go to the function handle_keyword for more flexibility
            const auto& kw_info = kwPos->second.num_value_per_cell(ncomps);
            this->handle_double_keyword(Section::SOLUTION, kw_info, keyword, box);
            continue;
        }

        this->handle_keyword(Section::SOLUTION, keyword, box);
    }
}

void FieldProps::handle_schedule_keywords(const std::vector<DeckKeyword>& keywords)
{
    auto box = makeGlobalGridBox(this->grid_ptr, &this->m_actnum, &this->m_active_index);

    // When called in the SCHEDULE section the context is that the scaling factors
    // have already been applied. We set them to zero for reuse here.
    for (const auto& [kw, _] : Fieldprops::keywords::SCHEDULE::double_keywords) {
        (void)_;
        if (this->has<double>(kw)) {
            auto& field_data = this->init_get<double>(kw);
            field_data.default_assign(1.0);
        }
    }

    for (const auto& keyword : keywords) {
        const std::string& name = keyword.name();

        if (auto kwPos = Fieldprops::keywords::SCHEDULE::double_keywords.find(name);
            kwPos != Fieldprops::keywords::SCHEDULE::double_keywords.end())
        {
            this->handle_double_keyword(Section::SCHEDULE, kwPos->second, keyword, box);
            continue;
        }

        if (Fieldprops::keywords::box_keywords.count(name) == 1) {
            handle_box_keyword(keyword, box);
            continue;
        }
    }
}

const std::string& FieldProps::default_region() const
{
    return this->m_default_region;
}

void FieldProps::apply_tran(const std::string& keyword, std::vector<double>& data)
{
    ::Opm::apply_tran(this->tran, this->double_data, this->active_size, keyword, data);
}


void FieldProps::apply_tranz_global(const std::vector<size_t>& indices, std::vector<double>& data) const
{
    ::Opm::apply_tran(this->tran.at("TRANZ"), this->double_data, indices, data);
}

bool FieldProps::tran_active(const std::string& keyword) const
{
    auto calculator = this->tran.find(keyword);
    return calculator != this->tran.end() && calculator->second.size() > 0;
}

void FieldProps::apply_numerical_aquifers(const NumericalAquifers& numerical_aquifers)
{
    auto& porv_data = this->init_get<double>("PORV").data;
    auto& poro_data = this->init_get<double>("PORO").data;
    auto& satnum_data = this->init_get<int>("SATNUM").data;
    auto& pvtnum_data = this->init_get<int>("PVTNUM").data;

    auto& permx_data = this->init_get<double>("PERMX").data;
    auto& permy_data = this->init_get<double>("PERMY").data;
    auto& permz_data = this->init_get<double>("PERMZ").data;

    const auto& aqu_cell_props = numerical_aquifers.aquiferCellProps();
    for (const auto& [global_index, cellprop] : aqu_cell_props) {
        const size_t active_index = this->grid_ptr->activeIndex(global_index);
        this->cell_volume[active_index] = cellprop.volume;
        this->cell_depth[active_index] = cellprop.depth;

        porv_data[active_index] = cellprop.pore_volume;
        poro_data[active_index] = cellprop.porosity;
        satnum_data[active_index] = cellprop.satnum;
        pvtnum_data[active_index] = cellprop.pvtnum;

        // isolate the numerical aquifer cells by setting permeability to be zero
        permx_data[active_index] = 0.;
        permy_data[active_index] = 0.;
        permz_data[active_index] = 0.;
    }
}

std::vector<std::string> FieldProps::fip_regions() const
{
    constexpr auto maxchars = std::string::size_type{6};

    std::vector<std::string> result;
    for (const auto& [key, fd_value] : this->int_data) {
        if (fd_value.valid() && Fieldprops::keywords::isFipxxx(key)) {
            result.push_back(key.substr(0, maxchars));
        }
    }

    return result;
}

void FieldProps::set_active_indices(const std::vector<int>& indices)
{
    m_active_index.clear();
    std::size_t idx = 0;
    for (int index : indices) {
        m_active_index.emplace(index, idx++);
    }
}

template std::vector<bool> FieldProps::defaulted<int>(const std::string& keyword);
template std::vector<bool> FieldProps::defaulted<double>(const std::string& keyword);
}
