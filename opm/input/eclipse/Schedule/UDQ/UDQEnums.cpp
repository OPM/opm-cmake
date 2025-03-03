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

#include <opm/input/eclipse/Schedule/UDQ/UDQEnums.hpp>

#include <opm/common/utility/String.hpp>

#include <charconv>
#include <map>
#include <optional>
#include <regex>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include <fmt/format.h>

namespace {

    std::string_view dequote(std::string_view s)
    {
        const auto b = s.find('\'');
        if (b == std::string_view::npos) {
            return s;
        }

        const auto e = s.find('\'', b + 1);
        if (e == std::string_view::npos) {
            throw std::invalid_argument {
                "Invalid quoted string |" + std::string{s} + '|'
            };
        }

        return s.substr(b + 1, e - b - 1);
    }

    bool is_asterisk(std::string_view s)
    {
        const auto ast = std::regex { R"(\s*\*\s*)" };
        return std::regex_match(s.begin(), s.end(), ast);
    }

    bool is_no_mix(const Opm::UDQVarType t)
    {
        return (t == Opm::UDQVarType::CONNECTION_VAR)
            || (t == Opm::UDQVarType::REGION_VAR)
            || (t == Opm::UDQVarType::SEGMENT_VAR)
            || (t == Opm::UDQVarType::AQUIFER_VAR)
            || (t == Opm::UDQVarType::BLOCK_VAR)
            || (t == Opm::UDQVarType::WELL_VAR)
            || (t == Opm::UDQVarType::GROUP_VAR)
            ;
    }

    bool is_valid_vartype(const Opm::UDQVarType t)
    {
        return is_no_mix(t)
            || (t == Opm::UDQVarType::NONE)
            || (t == Opm::UDQVarType::SCALAR)
            || (t == Opm::UDQVarType::FIELD_VAR)
            || (t == Opm::UDQVarType::TABLE_LOOKUP)
            ;
    }

    bool is_wgname_set(const std::vector<std::string>& selector)
    {
        return selector.empty()
            || (selector.front().find("*") != std::string::npos);
    }

    bool is_single_segment_number(std::string_view segNum0)
    {
        auto segNum = dequote(segNum0);

        if (segNum.empty()) {
            // Not specified
            return false;
        }

        auto result = 0;
        auto [ptr, ec] { std::from_chars(segNum.data(), segNum.data() + segNum.size(), result) };

        if ((ec == std::errc{}) && (ptr == segNum.data() + segNum.size())) {
            // Segment number is "123", or "'-1'", or something similar.
            return result > 0;
        }
        else if ((ec == std::errc::invalid_argument) && is_asterisk(segNum)) {
            // Segment number is '*'.  Treat as all segments.
            return false;
        }
        else {
            // Segment number is some unrecognized number string other than '*'.
            throw std::invalid_argument {
                "Invalid segment number string |" + std::string{segNum0} + '|'
            };
        }
    }

    bool is_segment_set(const std::vector<std::string>& selector)
    {
        return (selector.size() < std::vector<std::string>::size_type{2})
            || ! is_single_segment_number(selector[1])
            || selector.front().empty()
            || (selector.front().find("*") != std::string::npos);
    }

    bool is_single_region(const std::vector<std::string>& selector)
    {
        if (selector.size() != std::vector<std::string>::size_type{1}) {
            return false;
        }

        auto regNumStr = dequote(selector.front());
        if (regNumStr.empty()) {
            // Not specified
            return false;
        }

        const auto end = regNumStr.data() + regNumStr.size();
        auto regNum = 0;
        auto [ptr, ec] { std::from_chars(regNumStr.data(), end, regNum) };

        if ((ec == std::errc{}) && (ptr == end)) {
            return regNum > 0;
        }
        else if ((ec == std::errc::invalid_argument) && is_asterisk(regNumStr)) {
            // Region number is '*'.  Treat as all regions.
            return false;
        }
        else {
            // Region number is some unrecognized number string other than '*'.
            throw std::invalid_argument {
                "Invalid region number string |" + selector.front() + '|'
            };
        }
    }

    const auto cmp_func = std::set<Opm::UDQTokenType> {
        Opm::UDQTokenType::binary_cmp_eq,
        Opm::UDQTokenType::binary_cmp_ne,
        Opm::UDQTokenType::binary_cmp_le,
        Opm::UDQTokenType::binary_cmp_ge,
        Opm::UDQTokenType::binary_cmp_lt,
        Opm::UDQTokenType::binary_cmp_gt,
    };

    const auto binary_func = std::set<Opm::UDQTokenType> {
        Opm::UDQTokenType::binary_op_add,
        Opm::UDQTokenType::binary_op_mul,
        Opm::UDQTokenType::binary_op_sub,
        Opm::UDQTokenType::binary_op_div,
        Opm::UDQTokenType::binary_op_pow,
        Opm::UDQTokenType::binary_op_uadd,
        Opm::UDQTokenType::binary_op_umul,
        Opm::UDQTokenType::binary_op_umin,
        Opm::UDQTokenType::binary_op_umax,
        Opm::UDQTokenType::binary_cmp_eq,
        Opm::UDQTokenType::binary_cmp_ne,
        Opm::UDQTokenType::binary_cmp_le,
        Opm::UDQTokenType::binary_cmp_ge,
        Opm::UDQTokenType::binary_cmp_lt,
        Opm::UDQTokenType::binary_cmp_gt,
    };

    const auto set_func = std::set<Opm::UDQTokenType> {
        Opm::UDQTokenType::binary_op_uadd,
        Opm::UDQTokenType::binary_op_umul,
        Opm::UDQTokenType::binary_op_umin,
        Opm::UDQTokenType::binary_op_umax,
    };

    const auto scalar_func = std::set<Opm::UDQTokenType> {
        Opm::UDQTokenType::scalar_func_sum,
        Opm::UDQTokenType::scalar_func_avea,
        Opm::UDQTokenType::scalar_func_aveg,
        Opm::UDQTokenType::scalar_func_aveh,
        Opm::UDQTokenType::scalar_func_max,
        Opm::UDQTokenType::scalar_func_min,
        Opm::UDQTokenType::scalar_func_norm1,
        Opm::UDQTokenType::scalar_func_norm2,
        Opm::UDQTokenType::scalar_func_normi,
        Opm::UDQTokenType::scalar_func_prod,
    };

    const auto unary_elemental_func = std::set<Opm::UDQTokenType> {
        Opm::UDQTokenType::elemental_func_randn,
        Opm::UDQTokenType::elemental_func_randu,
        Opm::UDQTokenType::elemental_func_rrandn,
        Opm::UDQTokenType::elemental_func_rrandu,
        Opm::UDQTokenType::elemental_func_abs,
        Opm::UDQTokenType::elemental_func_def,
        Opm::UDQTokenType::elemental_func_exp,
        Opm::UDQTokenType::elemental_func_idv,
        Opm::UDQTokenType::elemental_func_ln,
        Opm::UDQTokenType::elemental_func_log,
        Opm::UDQTokenType::elemental_func_nint,
        Opm::UDQTokenType::elemental_func_sorta,
        Opm::UDQTokenType::elemental_func_sortd,
        Opm::UDQTokenType::elemental_func_undef,
    };

    const auto func_type = std::unordered_map<std::string, Opm::UDQTokenType> {
        {"+"    , Opm::UDQTokenType::binary_op_add},
        {"-"    , Opm::UDQTokenType::binary_op_sub},
        {"/"    , Opm::UDQTokenType::binary_op_div},
        {"DIV"  , Opm::UDQTokenType::binary_op_div},
        {"*"    , Opm::UDQTokenType::binary_op_mul},
        {"^"    , Opm::UDQTokenType::binary_op_pow},
        {"UADD" , Opm::UDQTokenType::binary_op_uadd},
        {"UMUL" , Opm::UDQTokenType::binary_op_umul},
        {"UMIN" , Opm::UDQTokenType::binary_op_umin},
        {"UMAX" , Opm::UDQTokenType::binary_op_umax},
        {"=="   , Opm::UDQTokenType::binary_cmp_eq},
        {"!="   , Opm::UDQTokenType::binary_cmp_ne},
        {"<="   , Opm::UDQTokenType::binary_cmp_le},
        {">="   , Opm::UDQTokenType::binary_cmp_ge},
        {"<"    , Opm::UDQTokenType::binary_cmp_lt},
        {">"    , Opm::UDQTokenType::binary_cmp_gt},
        {"RANDN", Opm::UDQTokenType::elemental_func_randn},
        {"RANDU", Opm::UDQTokenType::elemental_func_randu},
        {"RRNDN", Opm::UDQTokenType::elemental_func_rrandn},
        {"RRNDU", Opm::UDQTokenType::elemental_func_rrandu},
        {"ABS"  , Opm::UDQTokenType::elemental_func_abs},
        {"DEF"  , Opm::UDQTokenType::elemental_func_def},
        {"EXP"  , Opm::UDQTokenType::elemental_func_exp},
        {"IDV"  , Opm::UDQTokenType::elemental_func_idv},
        {"LN"   , Opm::UDQTokenType::elemental_func_ln},
        {"LOG"  , Opm::UDQTokenType::elemental_func_log},
        {"NINT" , Opm::UDQTokenType::elemental_func_nint},
        {"SORTA", Opm::UDQTokenType::elemental_func_sorta},
        {"SORTD", Opm::UDQTokenType::elemental_func_sortd},
        {"UNDEF", Opm::UDQTokenType::elemental_func_undef},
        {"SUM"  , Opm::UDQTokenType::scalar_func_sum},
        {"AVEA" , Opm::UDQTokenType::scalar_func_avea},
        {"AVEG" , Opm::UDQTokenType::scalar_func_aveg},
        {"AVEH" , Opm::UDQTokenType::scalar_func_aveh},
        {"MAX"  , Opm::UDQTokenType::scalar_func_max},
        {"MIN"  , Opm::UDQTokenType::scalar_func_min},
        {"NORM1", Opm::UDQTokenType::scalar_func_norm1},
        {"NORM2", Opm::UDQTokenType::scalar_func_norm2},
        {"NORMI", Opm::UDQTokenType::scalar_func_normi},
        {"PROD" , Opm::UDQTokenType::scalar_func_prod},
    };

} // Anonymous namespace

namespace Opm { namespace UDQ {

UDQVarType targetType(const std::string&              keyword,
                      const std::vector<std::string>& selector)
{
    const auto var_type = targetType(keyword);

    if (var_type == UDQVarType::NONE || var_type == UDQVarType::FIELD_VAR) {
        // Field variables and "none" (e.g., TIME) are treated as scalar.
        return UDQVarType::SCALAR;
    }

    if ((var_type == UDQVarType::WELL_VAR || var_type == UDQVarType::GROUP_VAR)
        && !is_wgname_set(selector))
    {
        // Well/group variables that apply to a single well/group are
        // treated as scalar.
        return UDQVarType::SCALAR;
    }

    if (var_type == UDQVarType::SEGMENT_VAR && !is_segment_set(selector)) {
        // Segment variables that apply to a single segment in a single MS
        // well are treated as scalar.
        return UDQVarType::SCALAR;
    }

    if (var_type == UDQVarType::REGION_VAR && is_single_region(selector)) {
        // Region variables which apply to a single region are treated as
        // scalars.
        return UDQVarType::SCALAR;
    }

    // All other targets keep the type inferred from the keyword.
    return var_type;
}

UDQVarType targetType(const std::string& keyword)
{
    if (keyword.empty()) {
        return UDQVarType::NONE;
    }

    if (keyword.substr(0,3) == "TU_") {
        return UDQVarType::TABLE_LOOKUP;
    }

    switch (keyword.front()) {
    case 'C': return UDQVarType::CONNECTION_VAR;
    case 'R': return UDQVarType::REGION_VAR;
    case 'F': return UDQVarType::FIELD_VAR;
    case 'S': return UDQVarType::SEGMENT_VAR;
    case 'A': return UDQVarType::AQUIFER_VAR;
    case 'B': return UDQVarType::BLOCK_VAR;
    case 'W': return UDQVarType::WELL_VAR;
    case 'G': return UDQVarType::GROUP_VAR;

    default:
        if (const auto double_value = try_parse_double(keyword);
            double_value.has_value())
        {
            return UDQVarType::SCALAR;
        }

        return UDQVarType::NONE;
    }
}

UDQVarType varType(const std::string& keyword)
{
    if ((keyword.size() < std::string::size_type{2}) || (keyword[1] != 'U')) {
        throw std::invalid_argument("Keyword: '" + keyword + "' is not of UDQ type");
    }

    const char first_char =  keyword[0];
    switch (first_char) {
    case 'W': return UDQVarType::WELL_VAR;
    case 'G': return UDQVarType::GROUP_VAR;
    case 'C': return UDQVarType::CONNECTION_VAR;
    case 'R': return UDQVarType::REGION_VAR;
    case 'F': return UDQVarType::FIELD_VAR;
    case 'S': return UDQVarType::SEGMENT_VAR;
    case 'A': return UDQVarType::AQUIFER_VAR;
    case 'B': return UDQVarType::BLOCK_VAR;

    default:
        throw std::invalid_argument("Keyword: " + keyword + " is not of UDQ type");
    }
}

UDQAction actionType(const std::string& action_string)
{
    if (action_string == "ASSIGN") {
        return UDQAction::ASSIGN;
    }

    if (action_string == "DEFINE") {
        return UDQAction::DEFINE;
    }

    if (action_string == "UNITS") {
        return UDQAction::UNITS;
    }

    if (action_string == "UPDATE") {
        return UDQAction::UPDATE;
    }

    throw std::invalid_argument("Invalid action string " + action_string);
}

UDQUpdate updateType(const std::string& update_string)
{
    if (update_string == "ON") {
        return UDQUpdate::ON;
    }

    if (update_string == "OFF") {
        return UDQUpdate::OFF;
    }

    if (update_string == "NEXT") {
        return UDQUpdate::NEXT;
    }

    throw std::invalid_argument("Invalid status update string " + update_string);
}

UDQUpdate updateType(const int int_value)
{
    switch (int_value) {
    case 0: return UDQUpdate::OFF;
    case 1: return UDQUpdate::NEXT;
    case 2: return UDQUpdate::ON;

    default:
        throw std::logic_error("Invalid integer for UDQUpdate type");
    }
}

bool binaryFunc(const UDQTokenType token_type)
{
    return binary_func.find(token_type) != binary_func.end();
}

bool scalarFunc(const UDQTokenType token_type)
{
    return scalar_func.find(token_type) != scalar_func.end();
}

bool elementalUnaryFunc(const UDQTokenType token_type)
{
    return unary_elemental_func.find(token_type) != unary_elemental_func.end();
}

bool cmpFunc(const UDQTokenType token_type)
{
    return cmp_func.find(token_type) != cmp_func.end();
}

bool setFunc(const UDQTokenType token_type)
{
    return set_func.find(token_type) != set_func.end();
}

UDQTokenType funcType(const std::string& func_name)
{
    {
        auto func = func_type.find(func_name);
        if (func != func_type.end()) {
            return func->second;
        }
    }

    if (func_name.substr(0, 2) == "TU") {
        return UDQTokenType::table_lookup;
    }

    return UDQTokenType::error;
}

UDQTokenType tokenType(const std::string& token)
{
    if (const auto token_type = funcType(token);
        token_type != UDQTokenType::error)
    {
        return token_type;
    }

    if (token == "(") {
        return UDQTokenType::open_paren;
    }

    if (token == "[") {
        return UDQTokenType::table_lookup_start;
    }

    if (token == ")") {
        return UDQTokenType::close_paren;
    }

    if (token == "]") {
        return UDQTokenType::table_lookup_end;
    }

    if (const auto number = try_parse_double(token);
        number.has_value())
    {
        return UDQTokenType::number;
    }

    return UDQTokenType::ecl_expr;
}

UDQVarType coerce(const UDQVarType t1, const UDQVarType t2)
{
    if (! (is_valid_vartype(t1) && is_valid_vartype(t2))) {
        // Note: Can't use typeName() here since that would throw another
        // exception.
        throw std::logic_error {
            fmt::format("Cannot coerce between {} and {}.",
                        static_cast<std::underlying_type_t<UDQVarType>>(t1),
                        static_cast<std::underlying_type_t<UDQVarType>>(t2))
        };
    }

    if (t1 == t2) {
        return t1;
    }

    const auto is_restricted_t1 = is_no_mix(t1);
    const auto is_restricted_t2 = is_no_mix(t2);

    if (is_restricted_t1 && is_restricted_t2) {
        // t1 != t2, but neither can be coerced into the other.
        throw std::logic_error {
            fmt::format("Cannot coerce between {} and {}",
                        typeName(t1), typeName(t2))
        };
    }

    if (is_restricted_t1) {
        return t1;
    }

    if (is_restricted_t2) {
        return t2;
    }

    if (t1 == UDQVarType::NONE) {
        return t2;
    }

    if (t2 == UDQVarType::NONE) {
        return t1;
    }

    return t1;
}

std::string typeName(const UDQVarType var_type)
{
    switch (var_type) {
    case UDQVarType::NONE:
        return "NONE";

    case UDQVarType::SCALAR:
        return "SCALAR";

    case UDQVarType::WELL_VAR:
        return "WELL_VAR";

    case UDQVarType::CONNECTION_VAR:
        return "CONNECTION_VAR";

    case UDQVarType::FIELD_VAR:
        return "FIELD_VAR";

    case UDQVarType::GROUP_VAR:
        return "GROUP_VAR";

    case UDQVarType::REGION_VAR:
        return "REGION_VAR";

    case UDQVarType::SEGMENT_VAR:
        return "SEGMENT_VAR";

    case UDQVarType::AQUIFER_VAR:
        return "AQUIFER_VAR";

    case UDQVarType::BLOCK_VAR:
        return "BLOCK_VAR";

    default:
        throw std::runtime_error {
            fmt::format("Should not be here: {}",
                        static_cast<std::underlying_type_t<UDQVarType>>(var_type))
        };
    }
}

bool trailingSpace(const UDQTokenType token_type)
{
    return binaryFunc(token_type)
        || cmpFunc(token_type);
}

bool leadingSpace(const UDQTokenType token_type)
{
    return binaryFunc(token_type)
        || cmpFunc(token_type);
}

namespace {
    template <typename Value>
    Value lookup_control_map_value(const std::map<UDAControl, Value>& map,
                                   const UDAControl                   control)
    {
        auto pos = map.find(control);
        if (pos == map.end()) {
            throw std::logic_error {
                fmt::format("Unrecognized enum type ({}) - internal error",
                            static_cast<std::underlying_type_t<UDAControl>>(control))
            };
        }

        return pos->second;
    }
} // Anonymous

UDAKeyword keyword(const UDAControl control)
{
    static const auto c2k = std::map<UDAControl, UDAKeyword> {
        {UDAControl::WCONPROD_ORAT, UDAKeyword::WCONPROD},
        {UDAControl::WCONPROD_WRAT, UDAKeyword::WCONPROD},
        {UDAControl::WCONPROD_GRAT, UDAKeyword::WCONPROD},
        {UDAControl::WCONPROD_LRAT, UDAKeyword::WCONPROD},
        {UDAControl::WCONPROD_RESV, UDAKeyword::WCONPROD},
        {UDAControl::WCONPROD_BHP,  UDAKeyword::WCONPROD},
        {UDAControl::WCONPROD_THP,  UDAKeyword::WCONPROD},
        {UDAControl::WCONPROD_LIFT, UDAKeyword::WCONPROD},

        // --------------------------------------------------------------
        {UDAControl::WCONINJE_RATE, UDAKeyword::WCONINJE},
        {UDAControl::WCONINJE_RESV, UDAKeyword::WCONINJE},
        {UDAControl::WCONINJE_BHP,  UDAKeyword::WCONINJE},
        {UDAControl::WCONINJE_THP,  UDAKeyword::WCONINJE},

        // --------------------------------------------------------------
        {UDAControl::WELTARG_ORAT, UDAKeyword::WELTARG},
        {UDAControl::WELTARG_WRAT, UDAKeyword::WELTARG},
        {UDAControl::WELTARG_GRAT, UDAKeyword::WELTARG},
        {UDAControl::WELTARG_LRAT, UDAKeyword::WELTARG},
        {UDAControl::WELTARG_RESV, UDAKeyword::WELTARG},
        {UDAControl::WELTARG_BHP,  UDAKeyword::WELTARG},
        {UDAControl::WELTARG_THP,  UDAKeyword::WELTARG},
        {UDAControl::WELTARG_LIFT, UDAKeyword::WELTARG},

        // --------------------------------------------------------------
        {UDAControl::GCONPROD_OIL_TARGET,    UDAKeyword::GCONPROD},
        {UDAControl::GCONPROD_WATER_TARGET,  UDAKeyword::GCONPROD},
        {UDAControl::GCONPROD_GAS_TARGET,    UDAKeyword::GCONPROD},
        {UDAControl::GCONPROD_LIQUID_TARGET, UDAKeyword::GCONPROD},

        // --------------------------------------------------------------
        {UDAControl::GCONINJE_SURFACE_MAX_RATE,      UDAKeyword::GCONINJE},
        {UDAControl::GCONINJE_RESV_MAX_RATE,         UDAKeyword::GCONINJE},
        {UDAControl::GCONINJE_TARGET_REINJ_FRACTION, UDAKeyword::GCONINJE},
        {UDAControl::GCONINJE_TARGET_VOID_FRACTION,  UDAKeyword::GCONINJE},
    };

    return lookup_control_map_value(c2k, control);
}

bool group_control(const UDAControl control)
{
    try {
        const auto kw = keyword(control);
        return (kw == UDAKeyword::GCONPROD)
            || (kw == UDAKeyword::GCONINJE);
    }
    catch (const std::logic_error&) {
        return false;
    }
}

bool well_control(const UDAControl control)
{
    try {
        const auto kw = keyword(control);
        return (kw == UDAKeyword::WCONPROD)
            || (kw == UDAKeyword::WCONINJE)
            || (kw == UDAKeyword::WELTARG);
    }
    catch (const std::logic_error&) {
        return false;
    }
}

bool is_well_injection_control(const UDAControl control,
                               const bool       isInjector)
{
    try {
        const auto kw = keyword(control);
        return (kw == UDAKeyword::WCONINJE)
            || (isInjector && (kw == UDAKeyword::WELTARG));
    }
    catch (const std::logic_error&) {
        return false;
    }
}

bool is_well_production_control(const UDAControl control,
                                const bool       isProducer)
{
    try {
        const auto kw = keyword(control);
        return (kw == UDAKeyword::WCONPROD)
            || (isProducer && (kw == UDAKeyword::WELTARG));
    }
    catch (const std::logic_error&) {
        return false;
    }
}

bool is_group_injection_control(const UDAControl control)
{
    try {
        return keyword(control) == UDAKeyword::GCONINJE;
    }
    catch (const std::logic_error&) {
        return false;
    }
}

bool is_group_production_control(const UDAControl control)
{
    try {
        return keyword(control) == UDAKeyword::GCONPROD;
    }
    catch (const std::logic_error&) {
        return false;
    }
}

int udaCode(const UDAControl control)
{
    static const auto c2uda = std::map<UDAControl, int> {
        {UDAControl::WCONPROD_ORAT,   300'004},
        {UDAControl::WCONPROD_WRAT,   400'004},
        {UDAControl::WCONPROD_GRAT,   500'004},
        {UDAControl::WCONPROD_LRAT,   600'004},
        {UDAControl::WCONPROD_RESV,   700'004},
        {UDAControl::WCONPROD_BHP,    800'004},
        {UDAControl::WCONPROD_THP,    900'004},
        {UDAControl::WCONPROD_LIFT, 1'100'004},

        // --------------------------------------------------------------
        {UDAControl::WCONINJE_RATE, 400'003},
        {UDAControl::WCONINJE_RESV, 500'003},
        {UDAControl::WCONINJE_BHP,  600'003},
        {UDAControl::WCONINJE_THP,  700'003},

        // --------------------------------------------------------------
        {UDAControl::GCONPROD_OIL_TARGET,    200'020},
        {UDAControl::GCONPROD_WATER_TARGET,  300'020},
        {UDAControl::GCONPROD_GAS_TARGET,    400'020},
        {UDAControl::GCONPROD_LIQUID_TARGET, 500'020},

        // --------------------------------------------------------------
        {UDAControl::GCONINJE_SURFACE_MAX_RATE,      300'018}, // Surface injection rate
        {UDAControl::GCONINJE_RESV_MAX_RATE,         400'018}, // Reservoir volume injection rate
        {UDAControl::GCONINJE_TARGET_REINJ_FRACTION, 500'018}, // Reinjection fraction
        {UDAControl::GCONINJE_TARGET_VOID_FRACTION,  600'018}, // Voidage replacement fraction

        // --------------------------------------------------------------
        {UDAControl::WELTARG_ORAT,        16},
        {UDAControl::WELTARG_WRAT,   100'016},
        {UDAControl::WELTARG_GRAT,   200'016},
        {UDAControl::WELTARG_LRAT,   300'016},
        {UDAControl::WELTARG_RESV,   400'016},
        {UDAControl::WELTARG_BHP,    500'016},
        {UDAControl::WELTARG_THP,    600'016},
        {UDAControl::WELTARG_LIFT, 1'000'016},
    };

    return lookup_control_map_value(c2uda, control);
}

namespace { namespace UDAVersionCompat {

    namespace v0 {
        UDAControl udaControl(const int uda_code)
        {
            switch (uda_code) {
            case   300'004: return UDAControl::WCONPROD_ORAT;
            case   400'004: return UDAControl::WCONPROD_WRAT;
            case   500'004: return UDAControl::WCONPROD_GRAT;
            case   600'004: return UDAControl::WCONPROD_LRAT;
            case   700'004: return UDAControl::WCONPROD_RESV;
            case   800'004: return UDAControl::WCONPROD_BHP;
            case   900'004: return UDAControl::WCONPROD_THP;
            case 1'100'004: return UDAControl::WCONPROD_LIFT;

            case   400'003: return UDAControl::WCONINJE_RATE;
            case   500'003: return UDAControl::WCONINJE_RESV;
            case   600'003: return UDAControl::WCONINJE_BHP;
            case   700'003: return UDAControl::WCONINJE_THP;

            case   200'019: return UDAControl::GCONPROD_OIL_TARGET;
            case   300'019: return UDAControl::GCONPROD_WATER_TARGET;
            case   400'019: return UDAControl::GCONPROD_GAS_TARGET;
            case   500'019: return UDAControl::GCONPROD_LIQUID_TARGET;

            case   300'017: return UDAControl::GCONINJE_SURFACE_MAX_RATE;
            case   400'017: return UDAControl::GCONINJE_RESV_MAX_RATE;
            case   500'017: return UDAControl::GCONINJE_TARGET_REINJ_FRACTION;
            case   600'017: return UDAControl::GCONINJE_TARGET_VOID_FRACTION;

            case        16: return UDAControl::WELTARG_ORAT;
            case   100'016: return UDAControl::WELTARG_WRAT;
            case   200'016: return UDAControl::WELTARG_GRAT;
            case   300'016: return UDAControl::WELTARG_LRAT;
            case   400'016: return UDAControl::WELTARG_RESV;
            case   500'016: return UDAControl::WELTARG_BHP;
            case   600'016: return UDAControl::WELTARG_THP;
            case 1'000'016: return UDAControl::WELTARG_LIFT;
            }

            throw std::logic_error {
                fmt::format("Unknown UDA integer control code {}", uda_code)
            };
        }
    } // namespace v0

    namespace vCurr {
        UDAControl udaControl(const int uda_code)
        {
            switch (uda_code) {
            case   300'004: return UDAControl::WCONPROD_ORAT;
            case   400'004: return UDAControl::WCONPROD_WRAT;
            case   500'004: return UDAControl::WCONPROD_GRAT;
            case   600'004: return UDAControl::WCONPROD_LRAT;
            case   700'004: return UDAControl::WCONPROD_RESV;
            case   800'004: return UDAControl::WCONPROD_BHP;
            case   900'004: return UDAControl::WCONPROD_THP;
            case 1'100'004: return UDAControl::WCONPROD_LIFT;

            case   400'003: return UDAControl::WCONINJE_RATE;
            case   500'003: return UDAControl::WCONINJE_RESV;
            case   600'003: return UDAControl::WCONINJE_BHP;
            case   700'003: return UDAControl::WCONINJE_THP;

            case   200'020: return UDAControl::GCONPROD_OIL_TARGET;
            case   300'020: return UDAControl::GCONPROD_WATER_TARGET;
            case   400'020: return UDAControl::GCONPROD_GAS_TARGET;
            case   500'020: return UDAControl::GCONPROD_LIQUID_TARGET;

            case   300'018: return UDAControl::GCONINJE_SURFACE_MAX_RATE;
            case   400'018: return UDAControl::GCONINJE_RESV_MAX_RATE;
            case   500'018: return UDAControl::GCONINJE_TARGET_REINJ_FRACTION;
            case   600'018: return UDAControl::GCONINJE_TARGET_VOID_FRACTION;

            case        16: return UDAControl::WELTARG_ORAT;
            case   100'016: return UDAControl::WELTARG_WRAT;
            case   200'016: return UDAControl::WELTARG_GRAT;
            case   300'016: return UDAControl::WELTARG_LRAT;
            case   400'016: return UDAControl::WELTARG_RESV;
            case   500'016: return UDAControl::WELTARG_BHP;
            case   600'016: return UDAControl::WELTARG_THP;
            case 1'000'016: return UDAControl::WELTARG_LIFT;
            }

            throw std::logic_error {
                fmt::format("Unknown UDA integer control code {}", uda_code)
            };
        }
    } // namespace vCurr

}} // namespace (anonymous)::UDAVersionCompat

RstFileUDAVersion rstFileUDAVersion(const int rstFileVersion)
{
    if (rstFileVersion == 201802) {
        return RstFileUDAVersion::v0;
    }
    else if (rstFileVersion == 202204) {
        return RstFileUDAVersion::vCurr;
    }

    return RstFileUDAVersion::vUnsupp;
}

UDAControl udaControl(const int uda_code, const RstFileUDAVersion udaVersion)
{
    switch (udaVersion) {
    case RstFileUDAVersion::v0:
        return UDAVersionCompat::v0::udaControl(uda_code);

    case RstFileUDAVersion::vCurr:
        return UDAVersionCompat::vCurr::udaControl(uda_code);

    default:
        throw std::invalid_argument {
            fmt::format("Restart file UDA version {} is unsupported",
                        static_cast<std::underlying_type_t<RstFileUDAVersion>>(udaVersion))
        };
    }
}

std::string controlName(const UDAControl control)
{
    switch (control) {
    case UDAControl::GCONPROD_OIL_TARGET:
        return "GCONPROD_ORAT";

    case UDAControl::GCONPROD_WATER_TARGET:
        return "GCONPROD_WRAT";

    case UDAControl::GCONPROD_GAS_TARGET:
        return "GCONPROD_GRAT";

    case UDAControl::GCONPROD_LIQUID_TARGET:
        return "GCONPROD_LRAT";

    case UDAControl::GCONINJE_SURFACE_MAX_RATE:
        return "GCONINJE_SURFACE_RATE";

    case UDAControl::GCONINJE_RESV_MAX_RATE:
        return "GCONINJE_RESERVOIR_RATE";

    case UDAControl::GCONINJE_TARGET_REINJ_FRACTION:
        return "GCONINJE_REINJ_FRACTION";

    case UDAControl::GCONINJE_TARGET_VOID_FRACTION:
        return "GCONINJE_VOID_FRACTION";

    case UDAControl::WCONPROD_ORAT:
        return "WCONPROD_ORAT";

    case UDAControl::WCONPROD_GRAT:
        return "WCONPROD_GRAT";

    case UDAControl::WCONPROD_WRAT:
        return "WCONPROD_WRAT";

    case UDAControl::WCONPROD_LRAT:
        return "WCONPROD_LRAT";

    case UDAControl::WCONPROD_RESV:
        return "WCONPROD_RESV";

    case UDAControl::WCONPROD_BHP:
        return "WCONPROD_BHP";

    case UDAControl::WCONPROD_THP:
        return "WCONPROD_THP";

    case UDAControl::WCONPROD_LIFT:
        return "WCONPROD_LIFT";

    case UDAControl::WCONINJE_RATE:
        return "WCONINJE_RATE";

    case UDAControl::WCONINJE_RESV:
        return "WCONINJE_RESV";

    case UDAControl::WCONINJE_BHP:
        return "WCONINJE_BHP";

    case UDAControl::WCONINJE_THP:
        return "WCONINJE_THP";

    case UDAControl::WELTARG_ORAT: return "WELTARG_ORAT";
    case UDAControl::WELTARG_WRAT: return "WELTARG_WRAT";
    case UDAControl::WELTARG_GRAT: return "WELTARG_GRAT";
    case UDAControl::WELTARG_LRAT: return "WELTARG_LRAT";
    case UDAControl::WELTARG_RESV: return "WELTARG_RESV";
    case UDAControl::WELTARG_BHP:  return "WELTARG_BHP";
    case UDAControl::WELTARG_THP:  return "WELTARG_THP";
    case UDAControl::WELTARG_LIFT: return "WELTARG_LIFT";
    }

    throw std::logic_error {
        fmt::format("Unknown UDA control keyword '{}'",
                    static_cast<std::underlying_type_t<UDAControl>>(control))
    };
}

}} // namespace Opm::UDQ
