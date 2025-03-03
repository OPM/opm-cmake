/*
  Copyright 2018 Statoil ASA

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

#include <opm/output/eclipse/AggregateUDQData.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/output/eclipse/AggregateGroupData.hpp>
#include <opm/output/eclipse/InteHEAD.hpp>
#include <opm/output/eclipse/UDQDims.hpp>
#include <opm/output/eclipse/VectorItems/intehead.hpp>
#include <opm/output/eclipse/VectorItems/udq.hpp>
#include <opm/output/eclipse/WriteRestartHelpers.hpp>

#include <opm/input/eclipse/Schedule/Group/Group.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/ScheduleState.hpp>

#include <opm/input/eclipse/Schedule/UDQ/UDQActive.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQAssign.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQConfig.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQDefine.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQEnums.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQFunctionTable.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQInput.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQParams.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQState.hpp>

#include <opm/input/eclipse/Schedule/Well/Well.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstring>
#include <iterator>
#include <map>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <fmt/format.h>

// ###########################################################################
// Class Opm::RestartIO::Helpers::AggregateUDQData
// ---------------------------------------------------------------------------

namespace VI = ::Opm::RestartIO::Helpers::VectorItems;

namespace {

    // maximum number of groups
    std::size_t ngmaxz(const std::vector<int>& inteHead)
    {
        return inteHead[20];
    }

    // maximum number of wells
    std::size_t nwmaxz(const std::vector<int>& inteHead)
    {
        return inteHead[163];
    }

    // function to return true if token is a function
    bool isTokenTypeFunc(const Opm::UDQTokenType token)
    {
        return Opm::UDQ::scalarFunc(token)
            || Opm::UDQ::elementalUnaryFunc(token)
            || (token == Opm::UDQTokenType::table_lookup);
    }

    // function to return true if token is a binary operator: type power
    // (exponentiation)
    bool isTokenTypeBinaryPowOp(const Opm::UDQTokenType token)
    {
        return token == Opm::UDQTokenType::binary_op_pow;
    }

    // function to return true if token is a binary operator: type multiply
    // or divide
    bool isTokenTypeBinaryMulDivOp(const Opm::UDQTokenType token)
    {
        const auto type_1 = std::array {
            Opm::UDQTokenType::binary_op_div,
            Opm::UDQTokenType::binary_op_mul,
        };

        return std::find(type_1.begin(), type_1.end(), token) != type_1.end();
    }

    // function to return true if token is a binary operator: type add or
    // subtract
    bool isTokenTypeBinaryAddSubOp(const Opm::UDQTokenType token)
    {
        const auto type_1 = std::array {
            Opm::UDQTokenType::binary_op_add,
            Opm::UDQTokenType::binary_op_sub,
        };

        return std::find(type_1.begin(), type_1.end(), token) != type_1.end();
    }

    // function to return true if token is a binary union operator
    bool isTokenTypeBinaryUnionOp(const Opm::UDQTokenType token)
    {
        const auto type_1 = std::array {
            Opm::UDQTokenType::binary_op_uadd,
            Opm::UDQTokenType::binary_op_umul,
            Opm::UDQTokenType::binary_op_umin,
            Opm::UDQTokenType::binary_op_umax,
        };

        return std::find(type_1.begin(), type_1.end(), token) != type_1.end();
    }

    // function to return true if token is an open or close parenthesis token
    bool isTokenTypeParen(const Opm::UDQTokenType token)
    {
        const auto type_1 = std::array {
            Opm::UDQTokenType::open_paren,
            Opm::UDQTokenType::close_paren,
        };

        return std::find(type_1.begin(), type_1.end(), token) != type_1.end();
    }

    // A function to return true if the token is an operator
    bool isOperatorToken(const Opm::UDQTokenType token)
    {
        return Opm::UDQ::scalarFunc(token)
            || Opm::UDQ::elementalUnaryFunc(token)
            || Opm::UDQ::binaryFunc(token)
            || Opm::UDQ::setFunc(token);
    }

    // function to return index number of last binary token not inside
    // bracket that is ending the expression
    int numOperators(const std::vector<Opm::UDQToken>& modTokens)
    {
        return std::count_if(modTokens.begin(), modTokens.end(),
                             [](const auto& token)
                             {
                                 return isOperatorToken(token.type())
                                     || isTokenTypeParen(token.type());
                             });
    }

    // function to return the precedence of the current operator/function
    int opFuncPrec(const Opm::UDQTokenType token)
    {
        int prec = 0;
        if (isTokenTypeFunc(token)) prec = 6;
        if (Opm::UDQ::cmpFunc(token)) prec = 5;
        if (isTokenTypeBinaryPowOp(token)) prec = 4;
        if (isTokenTypeBinaryMulDivOp(token)) prec = 3;
        if (isTokenTypeBinaryAddSubOp(token)) prec = 2;
        if (isTokenTypeBinaryUnionOp(token)) prec = 1;
        return prec;
    }

    struct substOuterParentheses
    {
        std::vector<Opm::UDQToken> highestLevOperators;
        std::map<std::size_t, std::vector<Opm::UDQToken>> substitutedTokens;
        int noleadingOpenPar;
        bool leadChangeSign;
    };

    // function to return
    //      a vector of functions and operators at the highest level,
    //      a map of substituted tokens,
    //      the number of leading open_paren that bracket the whole expression,
    //      a logical flag indicating whether there is a leading change of sign in the expression

    substOuterParentheses
    substitute_outer_parenthesis(const std::vector<Opm::UDQToken>& modTokens,
                                 int                               noLeadOpenPar,
                                 bool                              leadChgSgn)
    {
        std::map <std::size_t, std::vector<Opm::UDQToken>> substTok;
        std::vector<Opm::UDQToken> highLevOp;
        std::vector<std::size_t> startParen;
        std::vector<std::size_t> endParen;
        std::size_t level = 0;
        std::size_t search_pos = 0;

        while (search_pos < modTokens.size()) {
            if (modTokens[search_pos].type() == Opm::UDQTokenType::open_paren  && level == 0) {
                startParen.push_back(search_pos);
                ++level;
            }
            else if (modTokens[search_pos].type() == Opm::UDQTokenType::open_paren) {
                ++level;
            }
            else if (modTokens[search_pos].type() == Opm::UDQTokenType::close_paren && level == 1) {
                endParen.push_back(search_pos);
                --level;
            }
            else if (modTokens[search_pos].type() == Opm::UDQTokenType::close_paren) {
                --level;
            }

            ++search_pos;
        }


        //
        //Store all the operators at the highest level
        //include the ecl_expr tokens and replace the content of parentheses with a comp_expr
        if (startParen.size() >= 1) {
            if (startParen[0] > 0) {
                //First store all tokens before the first start_paren
                for (std::size_t i = 0; i < startParen[0]; ++i) {
                    highLevOp.emplace_back(modTokens[i]);
                }
            }

            //
            // Replace content of all parentheses at the highest level by an comp_expr
            // store all tokens including () for all tokens inside a pair of ()
            // also store the tokens between sets of () and at the end of an expression

            for (std::size_t ind = 0; ind < startParen.size();  ++ind) {
                std::vector<Opm::UDQToken> substringToken;
                for (std::size_t i = startParen[ind]; i < endParen[ind]+1; ++i) {
                    substringToken.emplace_back(modTokens[i]);
                }

                // store the content inside the parenthesis
                substTok.emplace(ind, std::move(substringToken));

                //
                // make the vector of high level tokens
                //
                //first add ecl_expr instead of content of (...)

                highLevOp.emplace_back(std::to_string(ind), Opm::UDQTokenType::comp_expr);
                //
                // store all tokens between end_paren before and start_paren after current ()
                std::size_t subS_max = (ind == startParen.size()-1) ? modTokens.size() : startParen[ind+1];

                if ((endParen[ind] + 1) < subS_max) {
                    for (std::size_t i = endParen[ind] + 1; i < subS_max; ++i) {
                        highLevOp.emplace_back(modTokens[i]);
                    }
                }
            }
        }
        else {
            //
            // treat the case with no ()
            for (std::size_t i = 0; i < modTokens.size(); ++i) {
                highLevOp.emplace_back(modTokens[i]);
            }
        }

        //
        // check if there is a leading minus-sign (change sign)
        if ((modTokens[0].type() == Opm::UDQTokenType::binary_op_sub)) {
            if (startParen.size() > 0) {
                // if followed by start_paren linked to end_paren before end of data
                // set flag and remove from operator list because it is considered as a highest precedence operator
                // unless () go from token 2 two the end of expression
                if ((startParen[0] == 1)  && (endParen[0] < modTokens.size()-1)) {
                    leadChgSgn = true;
                }
            }
            else {
                // set flag and remove from operator list
                leadChgSgn = true;
            }

            if (leadChgSgn) {
                // remove from operator list because it is considered as a highest precedence operator and is
                // therefore not a normal "binary_op_sub" operator
                std::vector<Opm::UDQToken> temp_high_lev_op(highLevOp.begin()+1, highLevOp.end());
                highLevOp = temp_high_lev_op;
            }
        }
        else if (startParen.size() >= 1) {
            //
            // check for leading start_paren combined with end_paren at end of data
            if ((startParen[0] == 0) && (endParen[0] == modTokens.size()-1)) {
                //
                // remove leading and trailing ()
                const std::vector<Opm::UDQToken> modTokens_red(modTokens.begin()+1, modTokens.end()-1);
                noLeadOpenPar += 1;
                //
                // recursive call to itself to re-interpret the token-input

                substOuterParentheses substOpPar =
                    substitute_outer_parenthesis(modTokens_red, noLeadOpenPar, leadChgSgn);

                highLevOp = substOpPar.highestLevOperators;
                substTok = substOpPar.substitutedTokens;
                noLeadOpenPar = substOpPar.noleadingOpenPar;
                leadChgSgn = substOpPar.leadChangeSign;
            }
        }
        // interpretation of token input is completed return resulting object

        return {
            highLevOp,
            substTok,
            noLeadOpenPar,
            leadChgSgn
        };
    }

    // Categorize function in terms of which token-types are used in formula
    //
    // The define_type is (-) the location among a set of tokens of the
    // "top" of the parse tree (AST - abstract syntax tree) i.e. the
    // location of the lowest precedence operator relative to the total set
    // of operators, functions and open-/close - parenthesis
    int define_type(const std::vector<Opm::UDQToken>& tokens)
    {
        int def_type = 0;
        int noLeadOpenPar = 0;
        bool leadChgSgn = false;

        //
        // analyse the expression

        substOuterParentheses expr = substitute_outer_parenthesis(tokens, noLeadOpenPar, leadChgSgn);

        //
        // loop over high level operators to find operator with lowest precedence and highest index

        int curPrec  = 100;
        std::size_t indLowestPrecOper = 0;
        for (std::size_t ind = 0; ind < expr.highestLevOperators.size(); ++ind) {
            if ((expr.highestLevOperators[ind].type() != Opm::UDQTokenType::ecl_expr) &&
                (expr.highestLevOperators[ind].type() != Opm::UDQTokenType::comp_expr) &&
                (expr.highestLevOperators[ind].type() != Opm::UDQTokenType::number))
            {
                const int tmpPrec = opFuncPrec(expr.highestLevOperators[ind].type());
                if (tmpPrec <= curPrec) {
                    curPrec = tmpPrec;
                    indLowestPrecOper = ind;
                }
            }
        }

        //
        // if lowest precedence operator is the first token (and not equal to change sign)
        // NOTE: also for the case with outer () removed
        if ((!expr.leadChangeSign) && (indLowestPrecOper == 0)) {
            // test if operator is a function (precedence = 6)
            if ((curPrec == 6) || (expr.highestLevOperators[indLowestPrecOper].type() == Opm::UDQTokenType::binary_op_sub) ) {
                def_type = -1;
                def_type -= expr.noleadingOpenPar;
            } else {
                // def type is 1 when for all other situations (ecl-expression or number)
                def_type = 1;
            }
        } else {
            //
            // treat cases which start either with (ecl_experessions, open-parenthes or leadChangeSign
            def_type = (expr.leadChangeSign) ?  -1 : 0;
            def_type -= expr.noleadingOpenPar;
            // calculate position of lowest precedence operator
            // account for leading change sign operator
            for (std::size_t ind = 0; ind <= indLowestPrecOper; ++ind) {
                //
                //count operators, including functions and parentheses (not original ecl_experessions)
                if (isOperatorToken(expr.highestLevOperators[ind].type())) {
                    // single operator - subtract one
                    --def_type;
                } else if (expr.highestLevOperators[ind].type() == Opm::UDQTokenType::comp_expr) {
                    // expression in parentheses -  add all operators
                    std::size_t ind_ce = static_cast<std::size_t>(std::stoi(expr.highestLevOperators[ind].str()));
                    auto indSubstTok = expr.substitutedTokens.find(ind_ce);
                    if (indSubstTok != expr.substitutedTokens.end()) {
                        // count the number of operators & parenthes in this sub-expression
                        def_type -= numOperators(indSubstTok->second);
                    } else {
                        const auto msg = fmt::format("Invalid compound expression index {}", ind_ce);
                        Opm::OpmLog::error(msg);
                        throw std::invalid_argument { msg };
                    }
                }
                else if ((expr.highestLevOperators[ind].type() != Opm::UDQTokenType::ecl_expr) &&
                         (expr.highestLevOperators[ind].type() != Opm::UDQTokenType::number))
                {
                    // unknown token - write warning
                    Opm::OpmLog::warning(fmt::format("Unknown tokenType '{}' in define_type()",
                                                     expr.highestLevOperators[ind].str()));
                }
            }
        }

        return def_type;
    }

    template <typename T>
    std::pair<bool, int>
    findInVector(const std::vector<T>& vecOfElements, const T& element)
    {
        std::pair<bool, int> result;

        // Find given element in vector
        auto it = std::find(vecOfElements.begin(), vecOfElements.end(), element);

        if (it != vecOfElements.end()) {
            result.second = std::distance(vecOfElements.begin(), it);
            result.first = true;
        }
        else {
            result.first = false;
            result.second = -1;
        }

        return result;
    }

    namespace iUdq {

        Opm::RestartIO::Helpers::WindowedArray<int>
        allocate(const Opm::UDQDims& udqDims)
        {
            using WV = Opm::RestartIO::Helpers::WindowedArray<int>;

            return WV {
                WV::NumWindows{ std::max(udqDims.totalNumUDQs(), std::size_t{1}) },
                WV::WindowSize{ Opm::UDQDims::entriesPerIUDQ() }
            };
        }

        template <class IUDQArray>
        void staticContrib(const Opm::UDQInput& udq_input, IUDQArray& iUdq)
        {
            if (udq_input.is<Opm::UDQDefine>()) {
                const auto& udq_define = udq_input.get<Opm::UDQDefine>();
                const auto& update_status =  udq_define.status();
                const auto& tokens = udq_define.tokens();

                iUdq[0] = (update_status.first == Opm::UDQUpdate::ON)
                    ? 2 : 0;

                iUdq[1] = define_type(tokens);
            }
            else {
                iUdq[0] = iUdq[1] = 0;
            }

            iUdq[2] = udq_input.index.typed_insert_index;
        }

    } // iUdq

    namespace iUad {

        template <class IUADArray>
        int staticContrib(const Opm::UDQActive::OutputRecord& iuad_record,
                          const bool                          is_field_uda,
                          const int                           iuap_offset,
                          IUADArray&                          iUad)
        {
            namespace VI = Opm::RestartIO::Helpers::VectorItems::IUad;

            using Ix   = VI::index;;
            using Kind = VI::Value::IuapElems;

            iUad[Ix::UDACode] = iuad_record.uda_code;

            // +1 for one-based indices.
            iUad[Ix::UDQIndex] = iuad_record.input_index + 1;

            iUad[Ix::NumIuapElm] = is_field_uda
                ? Kind::Field : Kind::Regular;

            iUad[Ix::UseCount] = iuad_record.use_count;

            // +1 for one-based indices.
            iUad[Ix::Offset] = iuap_offset + 1;

            return iUad[Ix::UseCount] * iUad[Ix::NumIuapElm];
        }

    } // iUad

    namespace zUdn {

        Opm::RestartIO::Helpers::WindowedArray<
            Opm::EclIO::PaddedOutputString<8>
        >
        allocate(const Opm::UDQDims& udqDims)
        {
            using WV = Opm::RestartIO::Helpers::WindowedArray<
                Opm::EclIO::PaddedOutputString<8>
                >;

            const auto nwin = std::max(udqDims.totalNumUDQs(), std::size_t{1});

            return WV {
                WV::NumWindows{ nwin },
                WV::WindowSize{ Opm::UDQDims::entriesPerZUDN() }
            };
        }

        template <class zUdnArray>
        void staticContrib(const Opm::UDQInput& udq_input, zUdnArray& zUdn)
        {
            using Ix = Opm::RestartIO::Helpers::VectorItems::ZUdn::index;;

            zUdn[Ix::Keyword] = udq_input.keyword();
            zUdn[Ix::Unit] = udq_input.unit();
        }

    } // zUdn

    namespace zUdl {

        Opm::RestartIO::Helpers::WindowedArray<
            Opm::EclIO::PaddedOutputString<8>
        >
        allocate(const Opm::UDQDims& udqDims)
        {
            using WV = Opm::RestartIO::Helpers::WindowedArray<
                Opm::EclIO::PaddedOutputString<8>
                >;

            const auto nwin = std::max(udqDims.totalNumUDQs(), std::size_t{1});

            return WV {
                WV::NumWindows{ nwin },
                WV::WindowSize{ Opm::UDQDims::entriesPerZUDL() }
            };
        }

        template <class zUdlArray>
        void staticContrib(const Opm::UDQInput& input, zUdlArray& zUdl)
        {
            // Write out the input formula if key is a DEFINE udq
            if (! input.is<Opm::UDQDefine>()) {
                return;
            }

            const auto l_sstr    = std::string::size_type {8};
            const auto max_l_str = Opm::UDQDims::entriesPerZUDL() * l_sstr;

            const auto& udq_define = input.get<Opm::UDQDefine>();
            const auto& z_data = udq_define.input_string();

            if (z_data.size() > max_l_str) {
                const auto msg =
                    fmt::format(R"(DEFINE expression for UDQ {} is too long.
  Number of characters {} exceeds upper limit of {}.
  Expression: {})",
                                udq_define.keyword(),
                                z_data.size(), max_l_str,
                                z_data);

                throw std::invalid_argument { msg };
            }

            const auto n_sstr = z_data.size() / l_sstr;
            for (auto i = 0*n_sstr; i < n_sstr; ++i) {
                if (i == 0) {
                    auto temp_str = z_data.substr(i * l_sstr, l_sstr);

                    // If first character is a minus sign, change to ~
                    if (temp_str.compare(0, 1, "-") == 0) {
                        temp_str.replace(0, 1, "~");
                    }

                    zUdl[i] = temp_str;
                }
                else {
                    zUdl[i] = z_data.substr(i * l_sstr, l_sstr);
                }
            }

            // Add remainder of last non-zero string
            if ((z_data.size() % l_sstr) > 0) {
                zUdl[n_sstr] = z_data.substr(n_sstr * l_sstr);
            }
        }

    } // zUdl

    namespace iGph {

        std::vector<int>
        phaseVector(const Opm::Schedule&    sched,
                    const std::size_t       simStep,
                    const std::vector<int>& inteHead)
        {
            auto inj_phase = std::vector<int>(ngmaxz(inteHead), 0);

            auto update_phase = [](const int phase, const int new_phase) {
                if (phase == 0) {
                    return new_phase;
                }

                throw std::logic_error {
                    "Cannot write restart files with UDA "
                    "control on multiple phases in same group"
                };
            };

            for (const auto* group : sched.restart_groups(simStep)) {
                if ((group == nullptr) || !group->isInjectionGroup()) {
                    continue;
                }

                auto& int_phase = (group->name() == "FIELD")
                    ? inj_phase.back()
                    : inj_phase[group->insert_index() - 1];

                int_phase = 0;
                for (const auto& [phase, int_value] : std::array {
                        std::pair {Opm::Phase::OIL,   1},
                        std::pair {Opm::Phase::WATER, 2},
                        std::pair {Opm::Phase::GAS,   3},
                    })
                {
                    if (! group->hasInjectionControl(phase)) {
                        continue;
                    }

                    if (group->injectionProperties(phase).uda_phase()) {
                        int_phase = update_phase(int_phase, int_value);
                    }
                }
            }

            return inj_phase;
        }

    } // iGph

    namespace iUap {

        std::vector<int>
        data(const Opm::ScheduleState&                       sched,
             const std::vector<Opm::UDQActive::InputRecord>& iuap)
        {
            // Construct the current list of well or group sequence numbers
            // to output the IUAP array.
            auto wg_no = std::vector<int>{};

            for (const auto& udaRecord : iuap) {
                switch (Opm::UDQ::keyword(udaRecord.control)) {
                case Opm::UDAKeyword::WCONPROD:
                case Opm::UDAKeyword::WCONINJE:
                case Opm::UDAKeyword::WELTARG:
                    // Well level control.  Use well's insertion index as
                    // the IUAP entry (+1 for one-based indices).
                    wg_no.push_back(sched.wells(udaRecord.wgname).seqIndex() + 1);
                    break;

                case Opm::UDAKeyword::GCONPROD:
                case Opm::UDAKeyword::GCONINJE: {
                    // Group level control.  Need to distinguish between the
                    // FIELD and the non-FIELD cases.

                    if (const auto& gname = udaRecord.wgname; gname != "FIELD") {
                        // The Schedule object inserts 'FIELD' at index
                        // zero.  The group's insert_index() is therefore,
                        // serendipitously, already suitably adjusted to
                        // one-based indices for output purposes.
                        wg_no.push_back(sched.groups(gname).insert_index());
                    }
                    else {
                        // IUAP for field level UDAs is represented by two
                        // copies of the numeric ID '1'.
                        wg_no.insert(wg_no.end(), 2, 1);
                    }
                }
                    break;

                default: {
                    const auto msg =
                        fmt::format("Invalid control keyword {} for UDQ {}",
                                    static_cast<std::underlying_type_t<Opm::UDAControl>>(udaRecord.control),
                                    udaRecord.udq);

                    Opm::OpmLog::error(msg);

                    throw std::invalid_argument { msg };
                }
                    break;
                }
            }

            return wg_no;
        }

    } // iUap

    namespace dUdf {

        std::optional<Opm::RestartIO::Helpers::WindowedArray<double>>
        allocate(const Opm::UDQDims& udqDims)
        {
            using WV = Opm::RestartIO::Helpers::WindowedArray<double>;

            auto dudf = std::optional<WV>{};

            if (udqDims.numFieldUDQs() > 0) {
                dudf.emplace(WV::NumWindows { udqDims.numFieldUDQs() },
                             WV::WindowSize { std::size_t{1} });
            }

            return dudf;
        }

        template <class DUDFArray>
        void staticContrib(const Opm::UDQState& udq_state,
                           const std::string&   udq,
                           DUDFArray&           dUdf)
        {
            // Set value for group name "FIELD"
            dUdf[0] = udq_state.has(udq)
                ? udq_state.get(udq)
                : Opm::UDQ::restart_default;
        }

    } // dUdf -- Field level UDQ values (DUDF restart array)

    namespace dUdg {

        std::optional<Opm::RestartIO::Helpers::WindowedArray<double>>
        allocate(const Opm::UDQDims& udqDims)
        {
            using WV = Opm::RestartIO::Helpers::WindowedArray<double>;

            auto dudg = std::optional<WV>{};

            if (udqDims.numGroupUDQs() > 0) {
                dudg.emplace(WV::NumWindows{ udqDims.numGroupUDQs() },
                             WV::WindowSize{ udqDims.maxNumGroups() });
            }

            return dudg;
        }

        template <class DUDGArray>
        void staticContrib(const Opm::UDQState&                  udq_state,
                           const std::vector<const Opm::Group*>& groups,
                           const std::string&                    udq,
                           const std::size_t                     ngmaxz,
                           DUDGArray&                            dUdg)
        {
            for (std::size_t ind = 0; ind < groups.size(); ++ind) {
                const auto* group = groups[ind];

                const auto useDflt = (group == nullptr)
                    || (ind == ngmaxz - 1)
                    || ! udq_state.has_group_var(group->name(), udq);

                dUdg[ind] = useDflt ? Opm::UDQ::restart_default
                    : udq_state.get_group_var(group->name(), udq);
            }
        }

    } // dUdg -- Group level UDQ values (DUDG restart array)

    namespace dUds {

        std::vector<std::string>
        allMsWells(const Opm::ScheduleState&       scheduleBlock,
                   const std::vector<std::string>& allWells)
        {
            auto msWells = std::vector<std::string>{};
            msWells.reserve(allWells.size());

            std::copy_if(allWells.begin(), allWells.end(),
                         std::back_inserter(msWells),
                         [&scheduleBlock](const std::string& wname)
                         {
                             auto wptr = scheduleBlock.wells.get_ptr(wname);
                             return (wptr != nullptr) && wptr->isMultiSegment();
                         });

            return msWells;
        }

        std::optional<Opm::RestartIO::Helpers::WindowedMatrix<double>>
        allocate(const Opm::UDQDims& udqDims)
        {
            using WM = Opm::RestartIO::Helpers::WindowedMatrix<double>;

            auto duds = std::optional<WM>{};

            if (udqDims.numSegmentUDQs() > 0) {
                // maxNumSegments() for each of
                //    maxNumMsWells() for each of
                //       numSegmentUDQs().
                //
                // Initial value restart_default simplifies collection logic.
                duds.emplace(WM::NumRows    { udqDims.numSegmentUDQs() },
                             WM::NumCols    { udqDims.maxNumMsWells()  },
                             WM::WindowSize { udqDims.maxNumSegments() },
                             Opm::UDQ::restart_default);
            }

            return duds;
        }

    } // dUds -- Segment level UDQ values (DUDS restart array)

    namespace dUdw {

        std::optional<Opm::RestartIO::Helpers::WindowedArray<double>>
        allocate(const Opm::UDQDims& udqDims)
        {
            using WV = Opm::RestartIO::Helpers::WindowedArray<double>;

            auto dudw = std::optional<WV>{};

            if (udqDims.numWellUDQs() > 0) {
                const auto numWells =
                    std::max(udqDims.maxNumWells(), std::size_t{1});

                dudw.emplace(WV::NumWindows{ udqDims.numWellUDQs() },
                             WV::WindowSize{ numWells });
            }

            return dudw;
        }

        template <class DUDWArray>
        void staticContrib(const Opm::UDQState&            udq_state,
                           const std::vector<std::string>& wells,
                           const std::string&              udq,
                           const std::size_t               nwmaxz,
                           DUDWArray&                      dUdw)
        {
            // Initialize array to the default value for the array
            std::fill_n(dUdw.begin(), nwmaxz, Opm::UDQ::restart_default);

            for (std::size_t ind = 0; ind < wells.size(); ++ind) {
                const auto& wname = wells[ind];

                if (udq_state.has_well_var(wname, udq)) {
                    dUdw[ind] = udq_state.get_well_var(wname, udq);
                }
            }
        }
    } // dUdw -- Well level UDQ values (DUDW restart array)
}

// ===========================================================================

Opm::RestartIO::Helpers::AggregateUDQData::
AggregateUDQData(const UDQDims& udqDims)
    : iUDQ_ { iUdq::allocate(udqDims) }
    , zUDN_ { zUdn::allocate(udqDims) }
    , zUDL_ { zUdl::allocate(udqDims) }
      // ------------------------------------------------------------
    , dUDF_ { dUdf::allocate(udqDims) }
    , dUDG_ { dUdg::allocate(udqDims) }
    , dUDS_ { dUds::allocate(udqDims) }
    , dUDW_ { dUdw::allocate(udqDims) }
{}

// ---------------------------------------------------------------------------

void
Opm::RestartIO::Helpers::AggregateUDQData::
captureDeclaredUDQData(const Schedule&         sched,
                       const std::size_t       simStep,
                       const UDQState&         udq_state,
                       const std::vector<int>& inteHead)
{
    const auto udqInput = sched.getUDQConfig(simStep).input();

    const auto allWells = this->dUDW_.has_value() || this->dUDS_.has_value()
        ? sched.wellNames(simStep)
        : std::vector<std::string>{};

    this->collectUserDefinedQuantities(udqInput, inteHead);

    this->collectUserDefinedArguments(sched, simStep, inteHead);

    if (this->dUDF_.has_value()) {
        this->collectFieldUDQValues(udqInput, udq_state,
                                    inteHead[VI::intehead::NO_FIELD_UDQS]);
    }

    if (this->dUDG_.has_value()) {
        this->collectGroupUDQValues(udqInput, udq_state, ngmaxz(inteHead),
                                    sched.restart_groups(simStep),
                                    inteHead[VI::intehead::NO_GROUP_UDQS]);
    }

    if (this->dUDS_.has_value()) {
        const auto msWells = dUds::allMsWells(sched[simStep], allWells);

        if (! msWells.empty()) {
            this->collectSegmentUDQValues(udqInput, udq_state, msWells);
        }
    }

    if (this->dUDW_.has_value()) {
        this->collectWellUDQValues(udqInput, udq_state, nwmaxz(inteHead),
                                   allWells,
                                   inteHead[VI::intehead::NO_WELL_UDQS]);
    }
}

// ---------------------------------------------------------------------------

void
Opm::RestartIO::Helpers::AggregateUDQData::
collectUserDefinedQuantities(const std::vector<UDQInput>& udqInput,
                             const std::vector<int>&      inteHead)
{
    const auto expectNumUDQ = inteHead[VI::intehead::NO_WELL_UDQS]
        + inteHead[VI::intehead::NO_GROUP_UDQS]
        + inteHead[VI::intehead::NO_FIELD_UDQS]
        + inteHead[VI::intehead::NO_SEG_UDQS];

    int cnt = 0;
    for (const auto& udq_input : udqInput) {
        const auto udq_index = udq_input.index.insert_index;

        auto iudq = this->iUDQ_[udq_index];
        iUdq::staticContrib(udq_input, iudq);

        auto zudn = this->zUDN_[udq_index];
        zUdn::staticContrib(udq_input, zudn);

        auto zudl = this->zUDL_[udq_index];
        zUdl::staticContrib(udq_input, zudl);

        ++cnt;
    }

    if (cnt != expectNumUDQ) {
        OpmLog::error(fmt::format("Inconsistent total number of UDQs: {}, "
                                  "and sum of field, group, segment, "
                                  "and well UDQs: {}", cnt, expectNumUDQ));
    }
}

// ---------------------------------------------------------------------------

void
Opm::RestartIO::Helpers::AggregateUDQData::
collectUserDefinedArguments(const Schedule&         sched,
                            const std::size_t       simStep,
                            const std::vector<int>& inteHead)
{
    const auto& udq_active = sched[simStep].udq_active();
    if (! udq_active) {
        // No UDAs at this report step.  Nothing to do.
        return;
    }

    // Collect UDAs into the IUAD, IUAP, and IGPH restart vectors.

    assert (inteHead[VI::intehead::NO_IUADS] > 0);

    this->collectIUAD(udq_active, inteHead[VI::intehead::NO_IUADS]);

    // 2. Form IUAP.
    this->collectIUAP(iUap::data(sched[simStep], udq_active.iuap()),
                      inteHead[VI::intehead::NO_IUAPS]);

    // 3. Form IGPH.
    this->collectIGPH(iGph::phaseVector(sched, simStep, inteHead),
                      inteHead[VI::intehead::NGMAXZ]);
}

// ---------------------------------------------------------------------------

void
Opm::RestartIO::Helpers::AggregateUDQData::
collectFieldUDQValues(const std::vector<UDQInput>& udqInput,
                      const UDQState&              udq_state,
                      const int                    expectNumFieldUDQs)
{
    auto ix = std::size_t {0};

    int cnt = 0;
    for (const auto& udq_input : udqInput) {
        if (udq_input.var_type() == UDQVarType::FIELD_VAR) {
            auto dudf = (*this->dUDF_)[ix];

            dUdf::staticContrib(udq_state, udq_input.keyword(), dudf);

            ++ix;
            ++cnt;
        }
    }

    if (cnt != expectNumFieldUDQs) {
        OpmLog::error(fmt::format("Inconsistent number of DUDF elements: {}, "
                                  "expected number of DUDF elements {}.",
                                  cnt, expectNumFieldUDQs));
    }
}

// ---------------------------------------------------------------------------

void
Opm::RestartIO::Helpers::AggregateUDQData::
collectGroupUDQValues(const std::vector<UDQInput>&     udqInput,
                      const UDQState&                  udqState,
                      const std::size_t                ngmax,
                      const std::vector<const Group*>& groups,
                      const int                        expectedNumGroupUDQs)
{
    auto ix = std::size_t{0};

    int cnt = 0;
    for (const auto& udq_input : udqInput) {
        if (udq_input.var_type() == UDQVarType::GROUP_VAR) {
            auto dudg = (*this->dUDG_)[ix];

            dUdg::staticContrib(udqState, groups,
                                udq_input.keyword(),
                                ngmax, dudg);

            ++ix;
            ++cnt;
        }
    }

    if (cnt != expectedNumGroupUDQs) {
        OpmLog::error(fmt::format("Inconsistent number of DUDG elements: {}, "
                                  "expected number of DUDG elements {}.",
                                  cnt, expectedNumGroupUDQs));
    }
}

// ---------------------------------------------------------------------------

void
Opm::RestartIO::Helpers::AggregateUDQData::
collectSegmentUDQValues(const std::vector<UDQInput>&    udqInput,
                        const UDQState&                 udqState,
                        const std::vector<std::string>& msWells)
{
    assert (msWells.size() <= this->dUDS_->numCols());

    auto udqIdx = WindowedMatrix<double>::Idx{0};
    for (const auto& udq_input : udqInput) {
        if (udq_input.var_type() != UDQVarType::SEGMENT_VAR) {
            continue;
        }

        if (udqIdx >= this->dUDS_->numRows()) {
            throw std::logic_error {
                fmt::format("UDQ variable index {} exceeds number "
                            "of declared segment level UDQs {}",
                            udqIdx, this->dUDS_->numRows())
            };
        }

        for (auto mswIdx = 0*msWells.size(); mswIdx < msWells.size(); ++mswIdx) {
            auto duds = (*this->dUDS_)(udqIdx, mswIdx);
            udqState.exportSegmentUDQ(udq_input.keyword(), msWells[mswIdx], duds);
        }

        ++udqIdx;
    }
}

// ---------------------------------------------------------------------------

void
Opm::RestartIO::Helpers::AggregateUDQData::
collectWellUDQValues(const std::vector<UDQInput>&    udqInput,
                     const UDQState&                 udqState,
                     const std::size_t               nwmax,
                     const std::vector<std::string>& wells,
                     const int                       expectedNumWellUDQs)
{
    auto ix = std::size_t {0};

    int cnt = 0;
    for (const auto& udq_input : udqInput) {
        if (udq_input.var_type() == UDQVarType::WELL_VAR) {
            auto dudw = (*this->dUDW_)[ix];

            dUdw::staticContrib(udqState, wells,
                                udq_input.keyword(),
                                nwmax, dudw);

            ++ix;
            ++cnt;
        }
    }

    if (cnt != expectedNumWellUDQs) {
        OpmLog::error(fmt::format("Inconsistent number of DUDW elements: {}, "
                                  "expected number of DUDW elements {}.",
                                  cnt, expectedNumWellUDQs));
    }
}

// ---------------------------------------------------------------------------

void
Opm::RestartIO::Helpers::AggregateUDQData::
collectIUAD(const UDQActive& udqActive, const std::size_t expectNumIUAD)
{
    const auto& iuad_records = udqActive.iuad();
    if (iuad_records.size() != expectNumIUAD) {
        OpmLog::error(fmt::format("Number of actual IUADs ({}) incommensurate "
                                  "with expected number of IUADs from INTEHEAD ({}).",
                                  iuad_records.size(), expectNumIUAD));

        return;
    }

    using WV = Opm::RestartIO::Helpers::WindowedArray<int>;

    this->iUAD_.emplace(WV::NumWindows{ expectNumIUAD },
                        WV::WindowSize{ Opm::UDQDims::entriesPerIUAD() });

    auto iuap_offset = 0;
    auto index = std::size_t{0};
    for (const auto& iuad_record : iuad_records) {
        const auto kw = UDQ::keyword(iuad_record.control);
        const auto is_field_uda =
            ((kw == UDAKeyword::GCONPROD) ||
             (kw == UDAKeyword::GCONINJE))
            && (iuad_record.wg_name() == "FIELD");

        auto iuad = (*this->iUAD_)[index];

        iuap_offset +=
            iUad::staticContrib(iuad_record, is_field_uda, iuap_offset, iuad);

        ++index;
    }
}

void
Opm::RestartIO::Helpers::AggregateUDQData::
collectIUAP(const std::vector<int>& wgIndex,
            const std::size_t       expectNumIUAP)
{
    if (wgIndex.size() != expectNumIUAP) {
        OpmLog::error(fmt::format("Number of actual IUAPs ({}) incommensurate "
                                  "with expected number of IUAPs from INTEHEAD ({}).",
                                  wgIndex.size(), expectNumIUAP));

        return;
    }

    using WV = Opm::RestartIO::Helpers::WindowedArray<int>;

    this->iUAP_.emplace(WV::NumWindows{ 1 }, WV::WindowSize{ expectNumIUAP });

    std::copy(wgIndex.begin(), wgIndex.end(), (*this->iUAP_)[0].begin());
}

void
Opm::RestartIO::Helpers::AggregateUDQData::
collectIGPH(const std::vector<int>& phase_vector,
            const std::size_t       expectNumIGPH)
{
    if (phase_vector.size() != expectNumIGPH) {
        OpmLog::error(fmt::format("Number of actual IGPHs ({}) incommensurate "
                                  "with expected number of IGPHs from INTEHEAD ({}).",
                                  phase_vector.size(), expectNumIGPH));

        return;
    }

    using WV = Opm::RestartIO::Helpers::WindowedArray<int>;

    this->iGPH_.emplace(WV::NumWindows{ 1 }, WV::WindowSize{ expectNumIGPH });

    std::copy(phase_vector.begin(), phase_vector.end(), (*this->iGPH_)[0].begin());
}
