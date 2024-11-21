/*
  Copyright 2016 Statoil ASA.

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

#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/utility/OpmInputError.hpp>
#include <opm/common/utility/shmatch.hpp>

#include <opm/input/eclipse/EclipseState/Aquifer/AquiferConfig.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FieldPropsManager.hpp>
#include <opm/input/eclipse/EclipseState/Grid/GridDims.hpp>

#include <opm/input/eclipse/Schedule/Group/Group.hpp>
#include <opm/input/eclipse/Schedule/Network/ExtNetwork.hpp>
#include <opm/input/eclipse/Schedule/MSW/WellSegments.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQConfig.hpp>
#include <opm/input/eclipse/Schedule/Well/Connection.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>

#include <opm/input/eclipse/Parser/ErrorGuard.hpp>
#include <opm/input/eclipse/Parser/InputErrorAction.hpp>
#include <opm/input/eclipse/Parser/ParseContext.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Deck/DeckItem.hpp>
#include <opm/input/eclipse/Deck/DeckKeyword.hpp>
#include <opm/input/eclipse/Deck/DeckRecord.hpp>
#include <opm/input/eclipse/Deck/DeckSection.hpp>

#include <opm/io/eclipse/EclUtil.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <iterator>
#include <map>
#include <regex>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <fmt/format.h>

namespace Opm {

namespace {

struct SummaryConfigContext {
    std::unordered_map<std::string, std::set<int>> regions;
};


    const std::vector<std::string> ALL_keywords = {
        "FAQR",  "FAQRG", "FAQT", "FAQTG", "FGIP", "FGIPG", "FGIPL",
        "FGIR",  "FGIT",  "FGOR", "FGPR",  "FGPT", "FOIP",  "FOIPG",
        "FOIPL", "FOIR",  "FOIT", "FOPR",  "FOPT", "FPR",   "FVIR",
        "FVIT",  "FVPR",  "FVPT", "FWCT",  "FWGR", "FWIP",  "FWIR",
        "FWIT",  "FWPR",  "FWPT",
        "GGIR",  "GGIT",  "GGOR", "GGPR",  "GGPT", "GOIR",  "GOIT",
        "GOPR",  "GOPT",  "GVIR", "GVIT",  "GVPR", "GVPT",  "GWCT",
        "GWGR",  "GWIR",  "GWIT", "GWPR",  "GWPT",
        "WBHP",  "WGIR",  "WGIT", "WGOR",  "WGPR", "WGPT",  "WOIR",
        "WOIT",  "WOPR",  "WOPT", "WPI",   "WTHP", "WVIR",  "WVIT",
        "WVPR",  "WVPT",  "WWCT", "WWGR",  "WWIR", "WWIT",  "WWPR",
        "WWPT",  "WGLIR",
        // ALL will not expand to these keywords yet
        // Analytical aquifer keywords
        "AAQR",  "AAQRG", "AAQT", "AAQTG"
    };

    const std::vector<std::string> GMWSET_keywords = {
        "GMWPT", "GMWPR", "GMWPA", "GMWPU", "GMWPG", "GMWPO", "GMWPS",
        "GMWPV", "GMWPP", "GMWPL", "GMWIT", "GMWIN", "GMWIA", "GMWIU", "GMWIG",
        "GMWIS", "GMWIV", "GMWIP", "GMWDR", "GMWDT", "GMWWO", "GMWWT"
    };

    const std::vector<std::string> FMWSET_keywords = {
        "FMCTF", "FMWPT", "FMWPR", "FMWPA", "FMWPU", "FMWPF", "FMWPO", "FMWPS",
        "FMWPV", "FMWPP", "FMWPL", "FMWIT", "FMWIN", "FMWIA", "FMWIU", "FMWIF",
        "FMWIS", "FMWIV", "FMWIP", "FMWDR", "FMWDT", "FMWWO", "FMWWT"
    };


    const std::vector<std::string> PERFORMA_keywords = {
        "TCPU", "ELAPSED","NEWTON","NLINEARS","NLINSMIN", "NLINSMAX","MLINEARS",
        "MSUMLINS","MSUMNEWT","TIMESTEP","TCPUTS","TCPUDAY","STEPTYPE","TELAPLIN"
    };

    const std::vector<std::string> NMESSAGE_keywords = {
        "MSUMBUG", "MSUMCOMM", "MSUMERR", "MSUMMESS", "MSUMPROB", "MSUMWARN"
    };

    const std::vector<std::string> DATE_keywords = {
         "DAY", "MONTH", "YEAR"
    };

    const std::map<std::string, std::vector<std::string>> meta_keywords = {
        {"PERFORMA", PERFORMA_keywords},
        {"NMESSAGE", NMESSAGE_keywords},
        {"DATE",     DATE_keywords},
        {"ALL",      ALL_keywords},
        {"FMWSET",   FMWSET_keywords},
        {"GMWSET",   GMWSET_keywords},
    };

    // This is a hardcoded mapping between 3D field keywords,
    // e.g. 'PRESSURE' and 'SWAT' and summary keywords like 'RPR' and
    // 'BPR'. The purpose of this mapping is to maintain an overview of
    // which 3D field keywords are needed by the Summary calculation
    // machinery, based on which summary keywords are requested.
    const std::map<std::string, std::set<std::string>> required_fields = {
         {"PRESSURE", {"FPR", "RPR*", "BPR"}},
         {"RPV",      {"FRPV", "RRPV*"}},
         {"OIP",      {"ROIP*", "FOIP", "FOE"}},
         {"OIPR",     {"FOIPR"}},
         {"OIPL",     {"ROIPL*", "FOIPL"}},
         {"OIPG",     {"ROIPG*", "FOIPG"}},
         {"GIP",      {"RGIP*", "FGIP" }},
         {"GIPR",     {"FGIPR"}},
         {"GIPL",     {"RGIPL*", "FGIPL"}},
         {"GIPG",     {"RGIPG*", "FGIPG"}},
         {"WIP",      {"RWIP*", "FWIP" }},
         {"WIPR",     {"FWIPR"}},
         {"WIPL",     {"RWIPL*", "FWIPL"}},
         {"WIPG",     {"RWIPG*", "FWIPG"}},
         {"WCD",      {"RWCD", "FWCD" }},
         {"GCDI",     {"RGCDI", "FGCDI"}},
         {"GCDM",     {"RGCDM", "FGCDM"}},
         {"GKDI",     {"RGKDI", "FGKDI"}},
         {"GKDM",     {"RGKDM", "FGKDM"}},
         {"SWAT",     {"BSWAT"}},
         {"SGAS",     {"BSGAS"}},
         {"SALT",     {"FSIP"}},
         {"TEMP",     {"BTCNFHEA"}},
         {"GMIP",     {"RGMIP", "FGMIP"}},
         {"GMGP",     {"RGMGP", "FGMGP"}},
         {"GMDS",     {"RGMDS", "FGMDS"}},
         {"GMTR",     {"RGMTR", "FGMTR"}},
         {"GMST",     {"RGMST", "FGMST"}},
         {"GMMO",     {"RGMMO", "FGMMO"}},
         {"GMUS",     {"RGMUS", "FGMUS"}},
         {"GKTR",     {"RGKTR", "FGKTR"}},
         {"GKMO",     {"RGKMO", "FGKMO"}},
    };

    using keyword_set = std::unordered_set<std::string>;

    inline bool is_in_set(const keyword_set& set, const std::string& keyword) {
        return set.find(keyword) != set.end();
    }

    bool is_special(const std::string& keyword) {
        static const keyword_set specialkw {
            "ELAPSED",
            "MAXDPR",
            "MAXDSG",
            "MAXDSO",
            "MAXDSW",
            "NAIMFRAC",
            "NEWTON",
            "NLINEARS",
            "NLINSMAX",
            "NLINSMIN",
            "STEPTYPE",
            "WNEWTON",
        };

        return is_in_set(specialkw, keyword);
    }

    bool is_udq_blacklist(const std::string& keyword) {
        static const keyword_set udq_blacklistkw {
            "SUMTHIN",
        };

        return is_in_set(udq_blacklistkw, keyword);
    }

    bool is_processing_instruction(const std::string& keyword) {
        static const keyword_set processing_instructionkw {
            "NARROW",
            "RPTONLY",
            "RUNSUM",
            "SEPARATE",
            "SUMMARY",
        };

        return is_in_set(processing_instructionkw, keyword);
    }

    bool is_udq(const std::string& keyword) {
        // Does 'keyword' match one of the patterns
        //   AU*, BU*, CU*, FU*, GU*, RU*, SU*, or WU*?
        using sz_t = std::string::size_type;
        return (keyword.size() > sz_t{1})
            && (keyword[1] == 'U')
            && !is_udq_blacklist(keyword)
            && (keyword.find_first_of("WGFCRBSA") == sz_t{0});
    }

    bool is_pressure(const std::string& keyword) {
        static const keyword_set presskw {
            "BHP", "BHPH", "THP", "THPH", "PR",
            "PRD", "PRDH", "PRDF", "PRDA",
            "AQP", "NQP",
        };

        return is_in_set(presskw, keyword.substr(1));
    }

    bool is_rate(const std::string& keyword) {
        static const keyword_set ratekw {
            "OPR", "GPR", "WPR", "GLIR", "LPR", "NPR", "CPR", "VPR", "TPR", "TPC",

            "OFR", "OFRF", "OFRS", "OFR+", "OFR-", "TFR",
            "GFR", "GFRF", "GFRS", "GFR+", "GFR-",
            "WFR", "WFR+", "WFR-",

            "OPGR", "GPGR", "WPGR", "VPGR",
            "OPRH", "GPRH", "WPRH", "LPRH",
            "OVPR", "GVPR", "WVPR",
            "OPRS", "GPRS", "OPRF", "GPRF",

            "OIR", "GIR", "WIR", "LIR", "NIR", "CIR", "VIR", "TIR", "TIC"
            "OIGR", "GIGR", "WIGR",
            "OIRH", "GIRH", "WIRH",
            "OVIR", "GVIR", "WVIR",
            "GMIR",

            "OPI", "OPP", "GPI", "GPP", "WPI", "WPP",

            "AQR", "AQRG", "NQR",
        };

        return is_in_set(ratekw, keyword.substr(1))
            || ((keyword.length() > 4) &&
                is_in_set({ "TPR", "TPC", "TIR", "TIC", "TFR" },
                          keyword.substr(1, 3)));
    }

    bool is_ratio(const std::string& keyword) {
        static const keyword_set ratiokw {
            "GLR" , "GOR" , "OGR", "WCT" , "WGR",
            "GLRH", "GORH",        "WCTH",
        };

        return is_in_set(ratiokw, keyword.substr(1));
    }

    bool is_total(const std::string& keyword) {
        static const keyword_set totalkw {
            "OPT", "GPT", "WPT", "LPT", "NPT", "CPT",
            "VPT", "TPT", "OVPT", "GVPT", "WVPT",
            "WPTH", "OPTH", "GPTH", "LPTH",
            "GPTS", "OPTS", "GPTF", "OPTF",

            "OFT", "OFT+", "OFT-", "OFTL", "OFTG",
            "GFT", "GFT+", "GFT-", "GFTL", "GFTG",
            "WFT", "WFT+", "WFT-",

            "WIT", "OIT", "GIT", "LIT", "NIT", "CIT", "VIT", "TIT",
            "WITH", "OITH", "GITH", "WVIT", "OVIT", "GVIT", "GMIT",

            "AQT", "AQTG", "NQT",
        };

        return is_in_set(totalkw, keyword.substr(1))
            || ((keyword.length() > 4) &&
                is_in_set({ "TPT", "TIT" },
                          keyword.substr(1, 3)));
    }

    bool is_count(const std::string& keyword) {
        static const keyword_set countkw {
            "MWIN", "MWIT", "MWPR", "MWPT"
        };

        return is_in_set(countkw, keyword);
    }

    bool is_control_mode(const std::string& keyword) {
        static const keyword_set countkw {
            "MCTP", "MCTW", "MCTG"
        };

        return (keyword == "WMCTL")
            || is_in_set(countkw, keyword.substr(1));
    }

    bool is_prod_index(const std::string& keyword) {
        static const keyword_set countkw {
            "PI", "PI1", "PI4", "PI5", "PI9",
            "PIO", "PIG", "PIW", "PIL",
        };

        return !keyword.empty()
            && ((keyword[0] == 'W') || (keyword[0] == 'C'))
            && is_in_set(countkw, keyword.substr(1));
    }

    bool is_supported_region_to_region(const std::string& keyword)
    {
        static const auto supported_kw = std::regex {
            R"~~(R[OGW]F[RT][-+GL_]?([A-Z0-9_]{3})?)~~"
        };

        // R[OGW]F[RT][-+GL]? (e.g., "ROFTG", "RGFR+", or "RWFT")
        return std::regex_match(keyword, supported_kw);
    }

    bool is_unsupported_region_to_region(const std::string& keyword)
    {
        static const auto unsupported_kw = std::regex {
            R"~~(R([EK]|NL)F[RT][-+_]?([A-Z0-9_]{3})?)~~"
        };

        // R[EK]F[RT][-+]? (e.g., "REFT" or "RKFR+")
        // RNLF[RT][-+]? (e.g., "RNLFR-" or "RNLFT")
        return std::regex_match(keyword, unsupported_kw);
    }

    bool is_region_to_region(const std::string& keyword)
    {
        return is_supported_region_to_region  (keyword)
            || is_unsupported_region_to_region(keyword);
    }

    bool is_aquifer(const std::string& keyword)
    {
        static const auto aqukw = keyword_set {
            "AQP", "AQR", "AQRG", "AQT", "AQTG",
            "LQR", "LQT", "LQRG", "LQTG",
            "NQP", "NQR", "NQT",
            "AQTD", "AQPD",
        };

        return (keyword.size() >= std::string::size_type{4})
            && (keyword[0] == 'A')
            && is_in_set(aqukw, keyword.substr(1));
    }

    bool is_numeric_aquifer(const std::string& keyword)
    {
        // ANQP, ANQR, ANQT
        return is_aquifer(keyword) && (keyword[1] == 'N');
    }

    bool is_connection_completion(const std::string& keyword)
    {
        static const auto conn_compl_kw = std::regex {
            R"(C[OGW][IP][RT]L)"
        };

        return std::regex_match(keyword, conn_compl_kw);
    }

    bool is_well_completion(const std::string& keyword)
    {
        static const auto well_compl_kw = std::regex {
            R"(W[OGWLV][PIGOLCF][RT]L([0-9_]{2}[0-9])?)"
        };

        // True, e.g., for WOPRL, WOPRL__8, WOPRL123, but not WOPRL___ or
        // WKITL.
        return std::regex_match(keyword, well_compl_kw);
    }

    bool is_well_comp(const std::string& keyword)
    {
        static const auto well_comp_kw = keyword_set {
             "WAMF", "WXMF", "WYMF", "WZMF", "WCGMR", "WCOMR"
        };

        return is_in_set(well_comp_kw, keyword);
    }

    bool is_node_keyword(const std::string& keyword)
    {
        static const auto nodekw = keyword_set {
            "GPR", "GPRG", "GPRW", "NPR", "GNETPR"
        };

        return is_in_set(nodekw, keyword);
    }

    bool need_node_names(const SUMMARYSection& sect)
    {
        // We need the the node names if there is any node-related summary
        // keywords in the input deck's SUMMARY section.  The reason is that
        // we need to be able to fill out all node names in the case of a
        // keyword that does not specify any nodes (e.g., "GPR /"), and to
        // check for missing nodes if a keyword is erroneously specified.

        return std::any_of(sect.begin(), sect.end(),
            [](const DeckKeyword& keyword)
        {
            return is_node_keyword(keyword.name());
        });
    }

    std::vector<std::string> collect_node_names(const Schedule& sched)
    {
        auto node_names = std::vector<std::string>{};
        auto names = std::unordered_set<std::string>{};

        const auto nstep = sched.size() - 1;
        for (auto step = 0*nstep; step < nstep; ++step) {
            const auto& nodes = sched[step].network.get().node_names();
            names.insert(nodes.begin(), nodes.end());
            // Also insert wells belonging to groups in the network to be able to report network-computed THPs
            for (const auto& node : nodes) {
                if (!sched.hasGroup(node, step)) continue;
                const auto& group = sched.getGroup(node, step);
                for (const std::string& wellname : group.wells()) {
                    names.insert(wellname);
                }
            }            
        }

        node_names.assign(names.begin(), names.end());
        std::sort(node_names.begin(), node_names.end());

        return node_names;
    }

    SummaryConfigNode::Category
    distinguish_group_from_node(const std::string& keyword)
    {
        return is_node_keyword(keyword)
            ? SummaryConfigNode::Category::Node
            : SummaryConfigNode::Category::Group;
    }

    SummaryConfigNode::Category
    distinguish_connection_from_completion(const std::string& keyword)
    {
        return is_connection_completion(keyword)
            ? SummaryConfigNode::Category::Completion
            : SummaryConfigNode::Category::Connection;
    }

    SummaryConfigNode::Category
    distinguish_well_from_completion(const std::string& keyword)
    {
        return is_well_completion(keyword)
            ? SummaryConfigNode::Category::Completion
            : SummaryConfigNode::Category::Well;
    }

void handleMissingWell( const ParseContext& parseContext, ErrorGuard& errors, const KeywordLocation& location, const std::string& well) {
    std::string msg_fmt = fmt::format("Request for missing well {} in {{keyword}}\n"
                                      "In {{file}} line {{line}}", well);
    parseContext.handleError( ParseContext::SUMMARY_UNKNOWN_WELL , msg_fmt, location, errors );
}


void handleMissingGroup( const ParseContext& parseContext , ErrorGuard& errors, const KeywordLocation& location, const std::string& group) {
    std::string msg_fmt = fmt::format("Request for missing group {} in {{keyword}}\n"
                                      "In {{file}} line {{line}}", group);
    parseContext.handleError( ParseContext::SUMMARY_UNKNOWN_GROUP , msg_fmt, location, errors );
}

void handleMissingNode( const ParseContext& parseContext, ErrorGuard& errors, const KeywordLocation& location, const std::string& node_name )
{
    std::string msg_fmt = fmt::format("Request for missing network node {} in {{keyword}}\n"
                                      "In {{file}} line {{line}}", node_name);
    parseContext.handleError( ParseContext::SUMMARY_UNKNOWN_NODE, msg_fmt, location, errors );
}

void handleMissingAquifer( const ParseContext& parseContext,
                           ErrorGuard& errors,
                           const KeywordLocation& location,
                           const int id,
                           const bool is_numeric)
{
    std::string msg_fmt = fmt::format("Request for missing {} aquifer {} in {{keyword}}\n"
                                      "In {{file}} line {{line}}",
                                      is_numeric ? "numeric" : "anlytic", id);
    parseContext.handleError(ParseContext::SUMMARY_UNKNOWN_AQUIFER, msg_fmt, location, errors);
}

inline void keywordW( SummaryConfig::keyword_list& list,
                      const std::vector<std::string>& well_names,
                      SummaryConfigNode baseWellParam)
{
    std::transform(well_names.begin(), well_names.end(),
                   std::back_inserter(list),
                   [&baseWellParam](const auto& wname)
                   {
                       return baseWellParam.namedEntity(wname);
                   });
}

inline void keywordAquifer( SummaryConfig::keyword_list& list,
                            const std::vector<int>& aquiferIDs,
                            SummaryConfigNode baseAquiferParam)
{
    std::transform(aquiferIDs.begin(), aquiferIDs.end(),
                   std::back_inserter(list),
                   [&baseAquiferParam](const auto& id)
                   {
                       return baseAquiferParam.number(id);
                   });
}

// later check whether parseContext and errors are required
// maybe loc will be needed
void keywordAquifer( SummaryConfig::keyword_list& list,
                     const std::vector<int>& analyticAquiferIDs,
                     const std::vector<int>& numericAquiferIDs,
                     const ParseContext& parseContext,
                     ErrorGuard& errors,
                     const DeckKeyword& keyword)
{
    /*
      The keywords starting with AL take as arguments a list of Aquiferlists -
      this is not supported at all.
    */
    if (keyword.name().find("AL") == std::string::size_type{0}) {
        Opm::OpmLog::warning(Opm::OpmInputError::format("Unhandled summary keyword {keyword}\n"
                                                        "In {file} line {line}", keyword.location()));
        return;
    }

    auto param = SummaryConfigNode {
        keyword.name(), SummaryConfigNode::Category::Aquifer, keyword.location()
    }
    .parameterType( parseKeywordType(keyword.name()) )
    .isUserDefined( is_udq(keyword.name()) );

    const auto is_numeric = is_numeric_aquifer(keyword.name());
    const auto& pertinentIDs = is_numeric
        ? numericAquiferIDs : analyticAquiferIDs;

    if (keyword.empty() ||
        ! keyword.getDataRecord().getDataItem().hasValue(0))
    {
        keywordAquifer(list, pertinentIDs, param);
    }
    else {
        auto ids = std::vector<int>{};
        auto end = pertinentIDs.end();

        for (const int id : keyword.getIntData()) {
            // Note: std::find() could be std::lower_bound() here, but we
            // typically expect the number of pertinent aquifer IDs to be
            // small (< 10) so there's no big gain from a log(N) algorithm
            // in the common case.
            if (std::find(pertinentIDs.begin(), end, id) == end) {
                handleMissingAquifer(parseContext, errors, keyword.location(), id, is_numeric);
                continue;
            }

            ids.push_back(id);
        }

        keywordAquifer(list, ids, param);
    }
}


inline std::array< int, 3 > getijk( const DeckRecord& record ) {
    return {{
        record.getItem( "I" ).get< int >( 0 ) - 1,
        record.getItem( "J" ).get< int >( 0 ) - 1,
        record.getItem( "K" ).get< int >( 0 ) - 1
    }};
}


inline void keywordCL(SummaryConfig::keyword_list& list,
                      const ParseContext& parseContext,
                      ErrorGuard& errors,
                      const DeckKeyword& keyword,
                      const Schedule& schedule,
                      const GridDims& dims)
{
    auto node = SummaryConfigNode {
        keyword.name(), SummaryConfigNode::Category::Completion, keyword.location()
    }
    .parameterType(parseKeywordType(keyword.name()))
    .isUserDefined(is_udq(keyword.name()));

    for (const auto& record : keyword) {
        const auto& pattern = record.getItem(0).get<std::string>(0);
        auto well_names = schedule.wellNames(pattern, schedule.size() - 1);

        if (well_names.empty()) {
            handleMissingWell(parseContext, errors, keyword.location(), pattern);
        }

        const auto ijk_defaulted = record.getItem(1).defaultApplied(0);
        for (const auto& wname : well_names) {
            const auto& well = schedule.getWellatEnd(wname);
            const auto& all_connections = well.getConnections();

            node.namedEntity(wname);
            if (ijk_defaulted) {
                std::transform(all_connections.begin(), all_connections.end(),
                               std::back_inserter(list),
                               [&node](const auto& conn)
                               {
                                    return node.number(1 + conn.global_index());
                               });
            } else {
                const auto& ijk = getijk(record);
                auto global_index = dims.getGlobalIndex(ijk[0], ijk[1], ijk[2]);

                if (all_connections.hasGlobalIndex(global_index)) {
                    const auto& conn = all_connections.getFromGlobalIndex(global_index);
                    list.push_back( node.number( 1 + conn.global_index()));
                } else {
                    std::string msg = fmt::format("Problem with keyword {{keyword}}\n"
                                                  "In {{file}} line {{line}}\n"
                                                  "Connection ({},{},{}) not defined for well {}",
                                                  ijk[0] + 1, ijk[1] + 1, ijk[2] + 1, wname);
                    parseContext.handleError( ParseContext::SUMMARY_UNHANDLED_KEYWORD, msg, keyword.location(), errors);
                }
            }
        }
    }
}

inline void keywordWL(SummaryConfig::keyword_list& list,
                      const ParseContext& parseContext,
                      ErrorGuard& errors,
                      const DeckKeyword& keyword,
                      const Schedule& schedule)
{
    for (const auto& record : keyword) {
        const auto& pattern = record.getItem(0).get<std::string>(0);
        const auto well_names = schedule.wellNames(pattern, schedule.size() - 1);

        if (well_names.empty()) {
            handleMissingWell(parseContext, errors, keyword.location(), pattern);
            continue;
        }

        const auto completion = record.getItem(1).get<int>(0);

        // Use an amended KEYWORDS entry incorporating the completion ID,
        // e.g. "WOPRL_12", for the W*L summary vectors.  This is special
        // case treatment for compatibility reasons as the more common entry
        // here would be to just use "keyword.name()".
        auto node = SummaryConfigNode {
            fmt::format("{}{:_>3}", keyword.name(), completion),
            SummaryConfigNode::Category::Completion, keyword.location()
        }
        .parameterType(parseKeywordType(keyword.name()))
        .isUserDefined(is_udq(keyword.name()))
        .number(completion);

        for (const auto& wname : well_names) {
            if (schedule.getWellatEnd(wname).hasCompletion(completion)) {
                list.push_back(node.namedEntity(wname));
            }
            else {
                const auto msg = fmt::format("Problem with keyword {{keyword}}\n"
                                             "In {{file}} line {{line}}\n"
                                             "Completion number {} not defined for well {}",
                                             completion, wname);
                parseContext.handleError(ParseContext::SUMMARY_UNHANDLED_KEYWORD,
                                         msg, keyword.location(), errors);
            }
        }
    }
}

inline void keywordW( SummaryConfig::keyword_list& list,
                      const std::string& keyword,
                      KeywordLocation loc,
                      const Schedule& schedule) {
    auto param = SummaryConfigNode {
        keyword, SummaryConfigNode::Category::Well , std::move(loc)
    }
    .parameterType( parseKeywordType(keyword) )
    .isUserDefined( is_udq(keyword) );

    keywordW( list, schedule.wellNames(), param );
}


inline void keywordW( SummaryConfig::keyword_list& list,
                      const ParseContext& parseContext,
                      ErrorGuard& errors,
                      const DeckKeyword& keyword,
                      const Schedule& schedule ) {
    if (is_well_completion(keyword.name()))
        return keywordWL(list, parseContext, errors, keyword, schedule);

    auto param = SummaryConfigNode {
        keyword.name(), SummaryConfigNode::Category::Well, keyword.location()
    }
    .parameterType( parseKeywordType(keyword.name()) )
    .isUserDefined( is_udq(keyword.name()) );

    if (!keyword.empty() && keyword.getDataRecord().getDataItem().hasValue(0)) {
        for( const std::string& pattern : keyword.getStringData()) {
          auto well_names = schedule.wellNames( pattern, schedule.size() - 1 );

            if( well_names.empty() )
                handleMissingWell( parseContext, errors, keyword.location(), pattern );

            keywordW( list, well_names, param );
        }
    } else
        keywordW( list, schedule.wellNames(), param );
}

void keywordG(SummaryConfig::keyword_list& list,
              const std::string&           keyword,
              const KeywordLocation&       loc,
              const Schedule&              schedule)
{
    auto param = SummaryConfigNode {
        keyword, SummaryConfigNode::Category::Group, loc
    }
    .parameterType(parseKeywordType(keyword))
    .isUserDefined(is_udq(keyword));

    for (const auto& group : schedule.groupNames() ) {
        if (group == "FIELD") { continue; }

        list.push_back(param.namedEntity(group));
    }
}

void keywordG(SummaryConfig::keyword_list& list,
              const ParseContext&          parseContext,
              ErrorGuard&                  errors,
              const DeckKeyword&           keyword,
              const Schedule&              schedule,
              const bool excludeFieldFromGroupKw = true)
{
    if (keyword.name() == "GMWSET") {
        return;
    }

    auto param = SummaryConfigNode {
        keyword.name(), SummaryConfigNode::Category::Group, keyword.location()
    }
    .parameterType(parseKeywordType(keyword.name()))
    .isUserDefined(is_udq(keyword.name()));

    if (keyword.empty() ||
        ! keyword.getDataRecord().getDataItem().hasValue(0))
    {
        for (const auto& group : schedule.groupNames()) {
            if (excludeFieldFromGroupKw && (group == "FIELD")) {
                continue;
            }

            list.push_back(param.namedEntity(group));
        }

        return;
    }

    const auto& item = keyword.getDataRecord().getDataItem();

    for (const auto& group : item.getData<std::string>()) {
        if (schedule.back().groups.has(group)) {
            list.push_back( param.namedEntity(group));
        }
        else {
            handleMissingGroup(parseContext, errors, keyword.location(), group);
        }
    }
}

void keyword_node( SummaryConfig::keyword_list& list,
                   const std::vector<std::string>& node_names,
                   const ParseContext& parseContext,
                   ErrorGuard& errors,
                   const DeckKeyword& keyword)
{
    if (node_names.empty()) {
        const auto& location = keyword.location();
        std::string msg = "The network node keyword {keyword} is not supported in runs without networks\n"
                          "In {file} line {line}";
        parseContext.handleError( ParseContext::SUMMARY_UNHANDLED_KEYWORD, msg, location, errors);
        return;
    }

    auto param = SummaryConfigNode {
        keyword.name(), SummaryConfigNode::Category::Node, keyword.location()
    }
    .parameterType( parseKeywordType(keyword.name()) )
    .isUserDefined( is_udq(keyword.name()) );

    if( keyword.empty() ||
        !keyword.getDataRecord().getDataItem().hasValue( 0 ) )
    {
        std::transform(node_names.begin(), node_names.end(),
                       std::back_inserter(list),
                       [&param](const auto& node_name)
                       {
                           return param.namedEntity(node_name);
                       });
        return;
    }

    const auto& item = keyword.getDataRecord().getDataItem();

    for (const auto& node_name : item.getData<std::string>()) {
        auto pos = std::find(node_names.begin(),
                             node_names.end(), node_name);
        if (pos != node_names.end())
            list.push_back( param.namedEntity(node_name) );
        else
            handleMissingNode( parseContext, errors, keyword.location(), node_name );
    }
}

inline void keywordF( SummaryConfig::keyword_list& list,
                      const std::string& keyword,
                      KeywordLocation loc) {
    auto param = SummaryConfigNode {
        keyword, SummaryConfigNode::Category::Field, std::move(loc)
    }
    .parameterType( parseKeywordType(keyword) )
    .isUserDefined( is_udq(keyword) );

    list.push_back( std::move(param) );
}

inline void keywordAquifer( SummaryConfig::keyword_list& list,
                            const std::string& keyword,
                            const std::vector<int>& analyticAquiferIDs,
                            const std::vector<int>& numericAquiferIDs,
                            KeywordLocation loc)
{
    auto param = SummaryConfigNode {
        keyword, SummaryConfigNode::Category::Aquifer, std::move(loc)
    }
    .parameterType( parseKeywordType(keyword) )
    .isUserDefined( is_udq(keyword) );

    const auto& pertinentIDs = is_numeric_aquifer(keyword)
        ? numericAquiferIDs
        : analyticAquiferIDs;

    keywordAquifer(list, pertinentIDs, param);
}

inline void keywordF( SummaryConfig::keyword_list& list,
                      const DeckKeyword& keyword ) {
    if( keyword.name() == "FMWSET" ) return;
    keywordF( list, keyword.name(), keyword.location() );
}


inline std::array< int, 3 > getijk( const Connection& completion ) {
    return { { completion.getI(), completion.getJ(), completion.getK() }};
}


inline void keywordB( SummaryConfig::keyword_list& list,
                      const DeckKeyword& keyword,
                      const GridDims& dims) {
    auto param = SummaryConfigNode {
        keyword.name(), SummaryConfigNode::Category::Block, keyword.location()
    }
    .parameterType( parseKeywordType(keyword.name()) )
    .isUserDefined( is_udq(keyword.name()) );

    auto isValid = [&dims](const std::array<int,3>& ijk)
    {
        return (static_cast<std::size_t>(ijk[0]) < dims.getNX())
            && (static_cast<std::size_t>(ijk[1]) < dims.getNY())
            && (static_cast<std::size_t>(ijk[2]) < dims.getNZ());
    };

  for( const auto& record : keyword ) {
      auto ijk = getijk( record );

      if (! isValid(ijk)) {
          const auto msg_fmt =
              fmt::format("Block level summary keyword "
                          "{{keyword}}\n"
                          "In {{file}} line {{line}}\n"
                          "References invalid cell "
                          "{},{},{} in grid of dimensions "
                          "{},{},{}.\nThis block summary "
                          "vector request is ignored.",
                          ijk[0] + 1  , ijk[1] + 1  , ijk[2] + 1,
                          dims.getNX(), dims.getNY(), dims.getNZ());

          OpmLog::warning(OpmInputError::format(msg_fmt, keyword.location()));

          continue;
      }

      int global_index = 1 + dims.getGlobalIndex(ijk[0], ijk[1], ijk[2]);
      list.push_back( param.number(global_index) );
  }
}

inline std::optional<std::string>
establishRegionContext(const DeckKeyword&       keyword,
                       const FieldPropsManager& field_props,
                       const ParseContext&      parseContext,
                       ErrorGuard&              errors,
                       SummaryConfigContext&    context)
{
    auto region_name = std::string { "FIPNUM" };

    if (keyword.name().size() > 5) {
        region_name = "FIP" + keyword.name().substr(5, 3);

        if (! field_props.has_int(region_name)) {
            const auto msg_fmt =
                fmt::format("Problem with summary keyword {{keyword}}\n"
                            "In {{file}} line {{line}}\n"
                            "FIP region {} not defined in "
                            "REGIONS section - {{keyword}} ignored", region_name);

            parseContext.handleError(ParseContext::SUMMARY_INVALID_FIPNUM,
                                     msg_fmt, keyword.location(), errors);
            return std::nullopt;
        }
    }

    if (context.regions.count(region_name) == 0) {
        const auto& fipnum = field_props.get_int(region_name);
        context.regions.emplace(std::piecewise_construct,
                                std::forward_as_tuple(region_name),
                                std::forward_as_tuple(fipnum.begin(), fipnum.end()));
    }

    return { region_name };
}

inline void keywordR2R_unsupported(const DeckKeyword&  keyword,
                                   const ParseContext& parseContext,
                                   ErrorGuard&         errors)
{
    const auto msg_fmt = std::string {
        "Region to region summary keyword {keyword} is ignored\n"
        "In {file} line {line}"
    };

    parseContext.handleError(ParseContext::SUMMARY_UNHANDLED_KEYWORD,
                             msg_fmt, keyword.location(), errors);
}

inline void keywordR2R(const DeckKeyword&           keyword,
                       const FieldPropsManager&     field_props,
                       const ParseContext&          parseContext,
                       ErrorGuard&                  errors,
                       SummaryConfigContext&        context,
                       SummaryConfig::keyword_list& list)
{
    if (is_unsupported_region_to_region(keyword.name())) {
        keywordR2R_unsupported(keyword, parseContext, errors);
    }

    if (is_udq(keyword.name())) {
        throw std::invalid_argument {
            "Inter-Region quantity '"
           + keyword.name() + "' "
           + "cannot be a user-defined quantity"
        };
    }

    const auto region_name = establishRegionContext(keyword, field_props,
                                                    parseContext, errors,
                                                    context);

    if (! region_name.has_value()) {
        return;
    }

    auto param = SummaryConfigNode {
        keyword.name(), SummaryConfigNode::Category::Region, keyword.location()
    }
    .parameterType(parseKeywordType(keyword.name()))
    .fip_region(region_name.value())
    .isUserDefined(false);

    // Expected format:
    //
    //   ROFT
    //     1 2 /
    //     1 4 /
    //   /
    for (const auto& record : keyword) {
        // We *intentionally* record/use one-based region IDs here.
        const auto r1 = record.getItem("REGION1").get<int>(0);
        const auto r2 = record.getItem("REGION2").get<int>(0);

        list.push_back(param.number(EclIO::combineSummaryNumbers(r1, r2)));
    }
}


void keywordR(SummaryConfig::keyword_list& list,
              SummaryConfigContext&        context,
              const DeckKeyword&           deck_keyword,
              const Schedule&              schedule,
              const FieldPropsManager&     field_props,
              const ParseContext&          parseContext,
              ErrorGuard&                  errors)
{
    const auto keyword = deck_keyword.name();
    if (is_region_to_region(keyword)) {
        keywordR2R(deck_keyword, field_props, parseContext, errors, context, list);
        return;
    }

    const auto region_name =
        establishRegionContext(deck_keyword, field_props,
                               parseContext, errors,
                               context);

    if (! region_name.has_value()) {
        return;
    }

    auto regions = std::vector<int>{};

    // Assume that the FIPNUM array contains the values {1,2,4}; i.e. the
    // maximum value is 4 and the value 3 is missing.  Values which are too
    // large, i.e., > 4 in this case - and values which are missing in the
    // range are treated differently:
    //
    //    region_id >= 5: The requested region results are completely ignored.
    //
    //    region_id == 3: The summary file will contain a vector Rxxx:3 with the
    //    value 0.
    //
    // These behaviors are closely tied to the implementation in
    // opm-simulators which actually performs the region summation; and that
    // is also the main reason to treat these quite similar error conditions
    // differently.

    if (! deck_keyword.empty() &&
        (deck_keyword.getDataRecord().getDataItem().data_size() > 0))
    {
        const auto& item = deck_keyword.getDataRecord().getDataItem();

        for (const auto& region_id : item.getData<int>()) {
            const auto& region_set = context.regions.at(region_name.value());
            auto max_iter = region_set.rbegin();
            if (region_id > *max_iter) {
                std::string msg_fmt = fmt::format("Problem with summary keyword {{keyword}}\n"
                                                  "In {{file}} line {{line}}\n"
                                                  "FIP region {} not present in {} - ignored", region_id, region_name.value());
                parseContext.handleError(ParseContext::SUMMARY_REGION_TOO_LARGE, msg_fmt, deck_keyword.location(), errors);
                continue;
            }

            if (region_set.count(region_id) == 0) {
                std::string msg_fmt = fmt::format("Problem with summary keyword {{keyword}}\n"
                                                  "In {{file}} line {{line}}\n"
                                                  "FIP region {} not present in {} - will use 0", region_id, region_name.value());
                parseContext.handleError(ParseContext::SUMMARY_EMPTY_REGION, msg_fmt, deck_keyword.location(), errors);
            }

            regions.push_back(region_id);
        }
    }
    else {
        const auto cregions = context.regions.at(region_name.value());
        regions.assign(cregions.begin(), cregions.end());
    }

    // See comment on function roew() in Summary.cpp for this weirdness.
    if (keyword.rfind("ROEW", 0) == 0) {
        auto copt_node = SummaryConfigNode("COPT", SummaryConfigNode::Category::Connection, {});
        for (const auto& wname : schedule.wellNames()) {
            copt_node.namedEntity(wname);

            const auto& well = schedule.getWellatEnd(wname);
            for( const auto& connection : well.getConnections() ) {
                copt_node.number( connection.global_index() + 1 );
                list.push_back( copt_node );
            }
        }
    }

    auto param = SummaryConfigNode {
        keyword, SummaryConfigNode::Category::Region, deck_keyword.location()
    }
    .parameterType(parseKeywordType(EclIO::SummaryNode::normalise_region_keyword(keyword)))
    .fip_region   (region_name.value())
    .isUserDefined(is_udq(keyword));

    std::transform(regions.begin(), regions.end(),
                   std::back_inserter(list),
                   [&param](const auto& region)
                   {
                       return param.number(region);
                   });
}


inline void keywordMISC( SummaryConfig::keyword_list& list,
                         const std::string& keyword,
                         KeywordLocation loc)
{
    if (meta_keywords.find(keyword) == meta_keywords.end())
        list.emplace_back( keyword, SummaryConfigNode::Category::Miscellaneous , std::move(loc));
}

inline void keywordMISC( SummaryConfig::keyword_list& list,
                         const DeckKeyword& keyword)
{
    keywordMISC(list, keyword.name(), keyword.location());
}

  inline void keywordC( SummaryConfig::keyword_list& list,
                        const ParseContext& parseContext,
                        ErrorGuard& errors,
                        const DeckKeyword& keyword,
                        const Schedule& schedule,
                        const GridDims& dims) {

    if (is_connection_completion(keyword.name()))
        return keywordCL(list, parseContext, errors, keyword, schedule, dims);

    auto param = SummaryConfigNode {
        keyword.name(), SummaryConfigNode::Category::Connection, keyword.location()
    }
    .parameterType( parseKeywordType( keyword.name()) )
    .isUserDefined( is_udq(keyword.name()) );

    for( const auto& record : keyword ) {

        const auto& wellitem = record.getItem( 0 );

        const auto well_names = wellitem.defaultApplied( 0 )
                              ? schedule.wellNames()
                              : schedule.wellNames( wellitem.getTrimmedString( 0 ) );
        const auto ijk_defaulted = record.getItem( 1 ).defaultApplied( 0 );

        if( well_names.empty() )
            handleMissingWell( parseContext, errors, keyword.location(), wellitem.getTrimmedString( 0 ) );

        for(const auto& name : well_names) {
            param.namedEntity(name);
            const auto& well = schedule.getWellatEnd(name);
            /*
             * we don't want to add connections that don't exist, so we iterate
             * over a well's connections regardless of the desired block is
             * defaulted or not
             */
            for( const auto& connection : well.getConnections() ) {
                auto cijk = getijk( connection );

                if( ijk_defaulted || ( cijk == getijk(record) ) ) {
                    const int global_index = 1 + dims.getGlobalIndex(cijk[0], cijk[1], cijk[2]);
                    list.push_back( param.number(global_index) );
                }
            }
        }
    }
}

    bool isKnownSegmentKeyword(const DeckKeyword& keyword)
    {
        const auto& kw = keyword.name();

        if ((kw == "SUMMARY") || (kw == "SUMTHIN")) {
            return false;
        }

        if (kw[1] == 'U') {
            // User-defined quantity at segment level.  Unbounded set, so
            // assume this is well defined.
            return true;
        }

        const auto kw_whitelist = std::vector<const char*> {
            "SOFR" , "SOFRF", "SOFRS", "SOFT", "SOFV", "SOHF", "SOVIS",
            "SGFR" , "SGFRF", "SGFRS", "SGFT", "SGFV", "SGHF", "SGVIS",
            "SWFR" ,                   "SWFT", "SWFV", "SWHF", "SWVIS",
            "SGOR" , "SOGR" , "SWCT" , "SWGR" ,
            "SODEN", "SGDEN", "SWDEN", "SMDEN", "SDENM",
            "SPR"  , "SPRD" , "SPRDH", "SPRDF", "SPRDA",
        };

        return std::any_of(kw_whitelist.begin(), kw_whitelist.end(),
            [&kw](const char* known)
        {
            return kw == known;
        }) 
        || is_in_set({ "STFR", "STFC" }, kw.substr(0, 4));
    }


    int maxNumWellSegments(const Well& well)
    {
        return well.isMultiSegment()
            ? well.getSegments().size() : 0;
    }

    void makeSegmentNodes(const int                    segID,
                          const DeckKeyword&           keyword,
                          const Well&                  well,
                          SummaryConfig::keyword_list& list)
    {
        if (!well.isMultiSegment()) {
            // Not an MSW.  Don't create summary vectors for segments.
            return;
        }

        auto param = SummaryConfigNode {
            keyword.name(), SummaryConfigNode::Category::Segment, keyword.location()
        }
        .namedEntity(well.name())
        .parameterType(parseKeywordType(keyword.name()))
        .isUserDefined(is_udq(keyword.name()));

        if (segID < 1) {
            // Segment number defaulted.  Allocate a summary vector for each
            // segment.
            const auto nSeg = maxNumWellSegments(well);

            for (auto segNumber = 0*nSeg; segNumber < nSeg; ++segNumber) {
                list.push_back(param.number(segNumber + 1));
            }
        }
        else {
            // Segment number specified.  Allocate single summary vector for
            // that segment number.
            list.push_back(param.number(segID));
        }
    }

    void keywordSNoRecords(const DeckKeyword&           keyword,
                           const Schedule&              schedule,
                           SummaryConfig::keyword_list& list)
    {
        // No keyword records.  Allocate summary vectors for all
        // segments in all wells at all times.
        //
        // Expected format:
        //
        //   SGFR
        //   / -- All segments in all MS wells at all times.

        const auto segID = -1;

        for (const auto& well : schedule.getWellsatEnd()) {
            makeSegmentNodes(segID, keyword, well, list);
        }
    }

    void keywordSWithRecords(const ParseContext&          parseContext,
                             ErrorGuard&                  errors,
                             const DeckKeyword&           keyword,
                             const Schedule&              schedule,
                             SummaryConfig::keyword_list& list)
    {
        // Keyword has explicit records.  Process those and create
        // segment-related summary vectors for those wells/segments
        // that match the description.
        //
        // Expected formats:
        //
        //   SOFR
        //     'W1'   1 /
        //     'W1'  10 /
        //     'W3'     / -- All segments
        //   /
        //
        //   SPR
        //     1*   2 / -- Segment 2 in all multi-segmented wells
        //   /

        for (const auto& record : keyword) {
            const auto& wellitem = record.getItem(0);
            const auto& well_names = wellitem.defaultApplied(0)
                ? schedule.wellNames()
                : schedule.wellNames(wellitem.getTrimmedString(0));

            if (well_names.empty()) {
                handleMissingWell(parseContext, errors, keyword.location(),
                                  wellitem.getTrimmedString(0));
            }

            // Negative 1 (< 0) if segment ID defaulted.  Defaulted
            // segment number in record implies all segments.
            const auto segID = record.getItem(1).defaultApplied(0)
                ? -1 : record.getItem(1).get<int>(0);

            for (const auto& well_name : well_names) {
                makeSegmentNodes(segID, keyword, schedule.getWellatEnd(well_name), list);
            }
        }
    }

    inline void keywordS(SummaryConfig::keyword_list& list,
                         const ParseContext&          parseContext,
                         ErrorGuard&                  errors,
                         const DeckKeyword&           keyword,
                         const Schedule&              schedule)
    {
        // Generate SMSPEC nodes for SUMMARY keywords of the form
        //
        //   SOFR
        //     'W1'   1 /
        //     'W1'  10 /
        //     'W3'     / -- All segments
        //   /
        //
        //   SPR
        //     1*   2 / -- Segment 2 in all multi-segmented wells
        //   /
        //
        //   SGFR
        //   / -- All segments in all MS wells at all times.

        if (! isKnownSegmentKeyword(keyword)) {
            // Ignore keywords that have not been explicitly white-listed
            // for treatment as segment summary vectors.
            return;
        }

        if (! keyword.empty()) {
            // Keyword with explicit records.  Handle as alternatives SOFR
            // and SPR above
            keywordSWithRecords(parseContext, errors,
                                keyword, schedule, list);
        }
        else {
            // Keyword with no explicit records.  Handle as alternative SGFR
            // above.
            keywordSNoRecords(keyword, schedule, list);
        }
    }

    std::string to_string(const SummaryConfigNode::Category cat)
    {
        switch (cat) {
        case SummaryConfigNode::Category::Aquifer:       return "Aquifer";
        case SummaryConfigNode::Category::Well:          return "Well";
        case SummaryConfigNode::Category::Group:         return "Group";
        case SummaryConfigNode::Category::Field:         return "Field";
        case SummaryConfigNode::Category::Region:        return "Region";
        case SummaryConfigNode::Category::Block:         return "Block";
        case SummaryConfigNode::Category::Connection:    return "Connection";
        case SummaryConfigNode::Category::Completion:    return "Completion";
        case SummaryConfigNode::Category::Segment:       return "Segment";
        case SummaryConfigNode::Category::Node:          return "Node";
        case SummaryConfigNode::Category::Miscellaneous: return "Miscellaneous";
        }

        throw std::invalid_argument {
            "Unhandled Summary Parameter Category '"
            + std::to_string(static_cast<int>(cat)) + '\''
        };
    }

    void check_udq( const KeywordLocation& location,
                    const Schedule& schedule,
                    const ParseContext& parseContext,
                    ErrorGuard& errors ) {
        if (! is_udq(location.keyword))
            // Nothing to do
            return;

        const auto& udq = schedule.getUDQConfig(schedule.size() - 1);

        if (!udq.has_keyword(location.keyword)) {
            std::string msg = "Summary output requested for UDQ {keyword}\n"
                              "In {file} line {line}\n"
                              "No definition for this UDQ found in the SCHEDULE section";
            parseContext.handleError(ParseContext::SUMMARY_UNDEFINED_UDQ, msg, location, errors);
            return;
        }

        if (!udq.has_unit(location.keyword)) {
            std::string msg = "Summary output requested for UDQ {keyword}\n"
                              "In {file} line {line}\n"
                              "No unit defined in the SCHEDULE section for {keyword}";
            parseContext.handleError(ParseContext::SUMMARY_UDQ_MISSING_UNIT, msg, location, errors);
        }
    }

void handleKW(SummaryConfig::keyword_list&    list,
              SummaryConfigContext&           context,
              const std::vector<std::string>& node_names,
              const std::vector<int>&         analyticAquiferIDs,
              const std::vector<int>&         numericAquiferIDs,
              const DeckKeyword&              keyword,
              const Schedule&                 schedule,
              const FieldPropsManager&        field_props,
              const ParseContext&             parseContext,
              ErrorGuard&                     errors,
              const GridDims&                 dims,
              const bool excludeFieldFromGroupKw = true)
{
    using Cat = SummaryConfigNode::Category;

    const auto& name = keyword.name();
    check_udq(keyword.location(), schedule, parseContext, errors);

    const auto cat = parseKeywordCategory(name);
    switch (cat) {
    case Cat::Well:
        if (is_well_comp(keyword.name())) {
            Opm::OpmLog::warning(Opm::OpmInputError::format("Unhandled summary keyword {keyword}\n"
                                                            "In {file} line {line}", keyword.location()));
            return;
        }

        keywordW(list, parseContext, errors, keyword, schedule);
        break;

    case Cat::Group:
        keywordG(list, parseContext, errors, keyword, schedule,
                 excludeFieldFromGroupKw);
        break;

    case Cat::Field:
        keywordF(list, keyword);
        break;

    case Cat::Block:
        keywordB(list, keyword, dims);
        break;

    case Cat::Region:
        keywordR(list, context, keyword, schedule, field_props, parseContext, errors);
        break;

    case Cat::Connection:
        keywordC(list, parseContext, errors, keyword, schedule, dims);
        break;

    case Cat::Completion:
        if (is_well_completion(name)) {
            keywordWL(list, parseContext, errors, keyword, schedule);
        }
        else {
            keywordCL(list, parseContext, errors, keyword, schedule, dims);
        }
        break;

    case Cat::Segment:
        keywordS(list, parseContext, errors, keyword, schedule);
        break;

    case Cat::Node:
        keyword_node(list, node_names, parseContext, errors, keyword);
        break;

    case Cat::Aquifer:
        keywordAquifer(list, analyticAquiferIDs, numericAquiferIDs, parseContext, errors, keyword);
        break;

    case Cat::Miscellaneous:
        keywordMISC(list, keyword);
        break;

    default: {
        const auto msg_fmt = fmt::format("Summary output keyword {{keyword}} of "
                                         "type {} is not supported\n"
                                         "In {{file}} line {{line}}",
                                         to_string(cat));

        parseContext.handleError(ParseContext::SUMMARY_UNHANDLED_KEYWORD,
                                 msg_fmt, keyword.location(), errors);
    }
        break;
    }
}

inline void handleKW( SummaryConfig::keyword_list& list,
                      const std::string& keyword,
                      const std::vector<int>& analyticAquiferIDs,
                      const std::vector<int>& numericAquiferIDs,
                      const KeywordLocation& location,
                      const Schedule& schedule,
                      const ParseContext& /* parseContext */,
                      ErrorGuard& /* errors */)
{
    if (is_udq(keyword)) {
        throw std::logic_error {
            "UDQ keywords not handleded when expanding alias list"
        };
    }

    using Cat = SummaryConfigNode::Category;
    const auto cat = parseKeywordCategory(keyword);

    switch (cat) {
    case Cat::Well:
        keywordW(list, keyword, location, schedule);
        break;

    case Cat::Group:
        keywordG(list, keyword, location, schedule);
        break;

    case Cat::Field:
        keywordF(list, keyword, location);
        break;

    case Cat::Aquifer:
        keywordAquifer(list, keyword, analyticAquiferIDs,
                       numericAquiferIDs, location);
        break;

    case Cat::Miscellaneous:
        keywordMISC(list, keyword, location);
        break;

    default:
        throw std::logic_error {
            fmt::format("Keyword type {} is not supported in alias "
                        "lists.  Internal error handling keyword {}",
                        to_string(cat), keyword)
        };
    }
}


  inline void uniq( SummaryConfig::keyword_list& vec ) {
      std::sort( vec.begin(), vec.end());
      auto logical_end = std::unique( vec.begin(), vec.end() );
      vec.erase( logical_end, vec.end() );
      if (vec.empty())
          return;

      /*
        This is a desperate hack to ensure that the ROEW keywords come after
        WOPT keywords, to ensure that the WOPT keywords have been fully
        evaluated in the SummaryState when we evaluate the ROEW keywords.
      */
      std::size_t tail_index = vec.size() - 1;
      std::size_t item_index = 0;
      while (true) {
          if (item_index >= tail_index)
              break;

          auto& node = vec[item_index];
          if (node.keyword().rfind("ROEW", 0) == 0) {
              std::swap( node, vec[tail_index] );
              tail_index--;
          }
          item_index++;
      }
  }

} // Anonymous namespace

// =====================================================================

SummaryConfigNode::Type parseKeywordType(std::string keyword)
{
    if (parseKeywordCategory(keyword) == SummaryConfigNode::Category::Region) {
        keyword = EclIO::SummaryNode::normalise_region_keyword(keyword);
    }

    if (is_well_completion(keyword)) {
        keyword.pop_back();
    }

    if (is_connection_completion(keyword)) {
        keyword.pop_back();
    }

    if (is_rate(keyword)) { return SummaryConfigNode::Type::Rate; }
    if (is_total(keyword)) { return SummaryConfigNode::Type::Total; }
    if (is_ratio(keyword)) { return SummaryConfigNode::Type::Ratio; }
    if (is_pressure(keyword)) { return SummaryConfigNode::Type::Pressure; }
    if (is_count(keyword)) { return SummaryConfigNode::Type::Count; }
    if (is_control_mode(keyword)) { return SummaryConfigNode::Type::Mode; }
    if (is_prod_index(keyword)) { return SummaryConfigNode::Type::ProdIndex; }

    return SummaryConfigNode::Type::Undefined;
}

SummaryConfigNode::Category parseKeywordCategory(const std::string& keyword)
{
    using Cat = SummaryConfigNode::Category;

    if (is_special(keyword)) { return Cat::Miscellaneous; }

    switch (keyword[0]) {
    case 'A': if (is_aquifer(keyword)) return Cat::Aquifer; break;
    case 'W': return distinguish_well_from_completion(keyword);
    case 'G': return distinguish_group_from_node(keyword);
    case 'F': return Cat::Field;
    case 'C': return distinguish_connection_from_completion(keyword);
    case 'R': return Cat::Region;
    case 'B': return Cat::Block;
    case 'S': return Cat::Segment;
    case 'N': return Cat::Node;
    }

    // TCPU, MLINEARS, NEWTON, &c
    return Cat::Miscellaneous;
}


SummaryConfigNode::SummaryConfigNode(std::string keyword, const Category cat, KeywordLocation loc_arg)
    : keyword_ (std::move(keyword))
    , category_(cat)
    , loc      (std::move(loc_arg))
{}

SummaryConfigNode SummaryConfigNode::serializationTestObject()
{
    SummaryConfigNode result;
    result.keyword_ = "test1";
    result.category_ = Category::Region;
    result.loc = KeywordLocation::serializationTestObject();
    result.type_ = Type::Pressure;
    result.name_ = "test2";
    result.number_ = 2;
    result.userDefined_ = true;

    return result;
}

SummaryConfigNode& SummaryConfigNode::fip_region(const std::string& fip_region)
{
    this->fip_region_ = fip_region;
    return *this;
}


SummaryConfigNode& SummaryConfigNode::parameterType(const Type type)
{
    this->type_ = type;
    return *this;
}

SummaryConfigNode& SummaryConfigNode::namedEntity(std::string name)
{
    this->name_ = std::move(name);
    return *this;
}

SummaryConfigNode& SummaryConfigNode::number(const int num)
{
    this->number_ = num;
    return *this;
}

SummaryConfigNode& SummaryConfigNode::isUserDefined(const bool userDefined)
{
    this->userDefined_ = userDefined;
    return *this;
}

std::string SummaryConfigNode::uniqueNodeKey() const
{
    switch (this->category()) {
    case SummaryConfigNode::Category::Well: [[fallthrough]];
    case SummaryConfigNode::Category::Node: [[fallthrough]];
    case SummaryConfigNode::Category::Group:
        return this->keyword() + ':' + this->namedEntity();

    case SummaryConfigNode::Category::Field: [[fallthrough]];
    case SummaryConfigNode::Category::Miscellaneous:
        return this->keyword();

    case SummaryConfigNode::Category::Aquifer: [[fallthrough]];
    case SummaryConfigNode::Category::Region: [[fallthrough]];
    case SummaryConfigNode::Category::Block:
        return this->keyword() + ':' + std::to_string(this->number());

    case SummaryConfigNode::Category::Connection: [[fallthrough]];
    case SummaryConfigNode::Category::Completion: [[fallthrough]];
    case SummaryConfigNode::Category::Segment:
        return this->keyword() + ':' + this->namedEntity() + ':' + std::to_string(this->number());
    }

    throw std::invalid_argument {
        "Unhandled Summary Parameter Category '"
        + to_string(this->category()) + '\''
    };
}

bool operator==(const SummaryConfigNode& lhs, const SummaryConfigNode& rhs)
{
    if (lhs.keyword() != rhs.keyword()) return false;

    assert (lhs.category() == rhs.category());

    switch( lhs.category() ) {
        case SummaryConfigNode::Category::Field: [[fallthrough]];
        case SummaryConfigNode::Category::Miscellaneous:
            // Fully identified by keyword
            return true;

        case SummaryConfigNode::Category::Well: [[fallthrough]];
        case SummaryConfigNode::Category::Node: [[fallthrough]];
        case SummaryConfigNode::Category::Group:
            // Equal if associated to same named entity
            return lhs.namedEntity() == rhs.namedEntity();

        case SummaryConfigNode::Category::Aquifer: [[fallthrough]];
        case SummaryConfigNode::Category::Region: [[fallthrough]];
        case SummaryConfigNode::Category::Block:
            // Equal if associated to same numeric entity
            return lhs.number() == rhs.number();

        case SummaryConfigNode::Category::Connection: [[fallthrough]];
        case SummaryConfigNode::Category::Completion: [[fallthrough]];
        case SummaryConfigNode::Category::Segment:
            // Equal if associated to same numeric
            // sub-entity of same named entity
            return (lhs.namedEntity() == rhs.namedEntity())
                && (lhs.number()      == rhs.number());
    }

    return false;
}

bool operator<(const SummaryConfigNode& lhs, const SummaryConfigNode& rhs)
{
    if (lhs.keyword() < rhs.keyword()) return true;
    if (rhs.keyword() < lhs.keyword()) return false;

    // If we get here, the keyword are equal.

    switch( lhs.category() ) {
        case SummaryConfigNode::Category::Field: [[fallthrough]];
        case SummaryConfigNode::Category::Miscellaneous:
            // Fully identified by keyword.
            // Return false for equal keywords.
            return false;

        case SummaryConfigNode::Category::Well: [[fallthrough]];
        case SummaryConfigNode::Category::Node: [[fallthrough]];
        case SummaryConfigNode::Category::Group:
            // Ordering determined by namedEntityd entity
            return lhs.namedEntity() < rhs.namedEntity();

        case SummaryConfigNode::Category::Aquifer: [[fallthrough]];
        case SummaryConfigNode::Category::Region: [[fallthrough]];
        case SummaryConfigNode::Category::Block:
            // Ordering determined by numeric entity
            return lhs.number() < rhs.number();

        case SummaryConfigNode::Category::Connection: [[fallthrough]];
        case SummaryConfigNode::Category::Completion: [[fallthrough]];
        case SummaryConfigNode::Category::Segment:
        {
            // Ordering determined by pair of named entity and numeric ID.
            //
            // Would ideally implement this in terms of operator< for
            // std::tuple<std::string,int>, with objects generated by std::tie(),
            // but `namedEntity()` does not return an lvalue.
            const auto& lnm = lhs.namedEntity();
            const auto& rnm = rhs.namedEntity();

            return ( lnm <  rnm)
                || ((lnm == rnm) && (lhs.number() < rhs.number()));
        }
    }

    throw std::invalid_argument {
        "Unhandled Summary Parameter Category '" + to_string(lhs.category()) + '\''
    };
}

// =====================================================================

SummaryConfig::SummaryConfig( const Deck& deck,
                              const Schedule& schedule,
                              const FieldPropsManager& field_props,
                              const AquiferConfig& aquiferConfig,
                              const ParseContext& parseContext,
                              ErrorGuard& errors,
                              const GridDims& dims) {
    try {
        SUMMARYSection section( deck );
        SummaryConfigContext context;

        const auto node_names = need_node_names(section)
            ? collect_node_names(schedule)
            : std::vector<std::string> {};

        const auto analyticAquifers = analyticAquiferIDs(aquiferConfig);
        const auto numericAquifers = numericAquiferIDs(aquiferConfig);

        for (const auto& kw : section) {
            if (is_processing_instruction(kw.name())) {
                handleProcessingInstruction(kw.name());
            } else {
                handleKW(this->m_keywords, context,
                         node_names, analyticAquifers, numericAquifers,
                         kw, schedule, field_props, parseContext, errors, dims);
            }
        }

        for (const auto& meta_pair : meta_keywords) {
            if (section.hasKeyword(meta_pair.first)) {
                const auto& deck_keyword = section.getKeyword(meta_pair.first);
                for (const auto& kw : meta_pair.second) {
                    if (!this->hasKeyword(kw)) {
                        KeywordLocation location = deck_keyword.location();
                        location.keyword = fmt::format("{}/{}", meta_pair.first, kw);

                        handleKW(this->m_keywords, kw,
                                 analyticAquifers, numericAquifers,
                                 location, schedule, parseContext, errors);
                    }
                }
            }
        }

        uniq(this->m_keywords);
        for (const auto& kw : this->m_keywords) {
            this->short_keywords.insert(kw.keyword());
            this->summary_keywords.insert(kw.uniqueNodeKey());
        }
    }
    catch (const OpmInputError& opm_error) {
        throw;
    }
    catch (const std::exception& std_error) {
        OpmLog::error(fmt::format("An error occurred while configuring the summary properties\n"
                                  "Internal error: {}", std_error.what()));
        throw;
    }
}


SummaryConfig::SummaryConfig( const Deck& deck,
                              const Schedule& schedule,
                              const FieldPropsManager& field_props,
                              const AquiferConfig& aquiferConfig,
                              const ParseContext& parseContext,
                              ErrorGuard& errors) :
    SummaryConfig( deck , schedule, field_props, aquiferConfig, parseContext, errors, GridDims( deck ))
{ }


template <typename T>
SummaryConfig::SummaryConfig( const Deck& deck,
                              const Schedule& schedule,
                              const FieldPropsManager& field_props,
                              const AquiferConfig& aquiferConfig,
                              const ParseContext& parseContext,
                              T&& errors) :
    SummaryConfig(deck, schedule, field_props, aquiferConfig, parseContext, errors)
{}


SummaryConfig::SummaryConfig( const Deck& deck,
                              const Schedule& schedule,
                              const FieldPropsManager& field_props,
                              const AquiferConfig& aquiferConfig) :
    SummaryConfig(deck, schedule, field_props, aquiferConfig, ParseContext(), ErrorGuard())
{}


SummaryConfig::SummaryConfig(const keyword_list& kwds,
                             const std::set<std::string>& shortKwds,
                             const std::set<std::string>& smryKwds) :
    m_keywords(kwds), short_keywords(shortKwds), summary_keywords(smryKwds)
{}

SummaryConfig SummaryConfig::serializationTestObject()
{
    SummaryConfig result;
    result.m_keywords = {SummaryConfigNode::serializationTestObject()};
    result.short_keywords = {"test1"};
    result.summary_keywords = {"test2"};

    return result;
}

SummaryConfig::const_iterator SummaryConfig::begin() const {
    return this->m_keywords.cbegin();
}

SummaryConfig::const_iterator SummaryConfig::end() const {
    return this->m_keywords.cend();
}

SummaryConfig& SummaryConfig::merge( const SummaryConfig& other ) {
    this->m_keywords.insert( this->m_keywords.end(),
                             other.m_keywords.begin(),
                             other.m_keywords.end() );

    uniq( this->m_keywords );
    return *this;
}

SummaryConfig& SummaryConfig::merge( SummaryConfig&& other ) {
    auto fst = std::make_move_iterator( other.m_keywords.begin() );
    auto lst = std::make_move_iterator( other.m_keywords.end() );
    this->m_keywords.insert( this->m_keywords.end(), fst, lst );
    other.m_keywords.clear();

    uniq( this->m_keywords );
    return *this;
}

SummaryConfig::keyword_list
SummaryConfig::registerRequisiteUDQorActionSummaryKeys(const std::vector<std::string>& extraKeys,
                                                       const EclipseState&             es,
                                                       const Schedule&                 sched)
{
    auto summaryNodes = keyword_list{};

    if (extraKeys.empty()) {
        return summaryNodes;
    }

    auto candidateSummaryNodes = keyword_list{};

    // Note: When handling UDQs or, especially, ACTIONX condition blocks, it
    // is permissible to treat 'FIELD' as a regular group.  In particular,
    // an ACTIONX condition block might use conditions such as
    //
    //    GGOR 'FIELD' > 123.4 AND /
    //    GOPR 'FIELD' < 654.3 /
    //
    // and we need to be prepared to handle those.  We therefore bypass the
    // check for group == "FIELD" in keywordG() in this particular context.
    {
        const auto excludeFieldFromGroupKw = false;

        const auto node_names =
            std::any_of(extraKeys.begin(), extraKeys.end(), &is_node_keyword)
            ? collect_node_names(sched)
            : std::vector<std::string>{};

        const auto analyticAquifers = analyticAquiferIDs(es.aquifer());
        const auto numericAquifers = numericAquiferIDs(es.aquifer());

        const auto parseCtx = ParseContext { InputErrorAction::IGNORE };
        auto errors = ErrorGuard {};
        auto ctxt   = SummaryConfigContext {};

        for (const auto& vector_name : extraKeys) {
            handleKW(candidateSummaryNodes, ctxt, node_names,
                     analyticAquifers, numericAquifers,
                     DeckKeyword { KeywordLocation{}, vector_name },
                     sched, es.globalFieldProps(),
                     parseCtx, errors, es.gridDims(),
                     excludeFieldFromGroupKw);
        }
    }

    std::sort(candidateSummaryNodes.begin(), candidateSummaryNodes.end());

    summaryNodes.reserve(candidateSummaryNodes.size());
    std::set_difference(candidateSummaryNodes.begin(), candidateSummaryNodes.end(),
                        this->m_keywords.begin(), this->m_keywords.end(),
                        std::back_inserter(summaryNodes));

    if (summaryNodes.empty()) {
        // No new summary keywords encountered.
        return summaryNodes;
    }

    for (const auto& newKw : summaryNodes) {
        this->short_keywords.insert(newKw.keyword());
        this->summary_keywords.insert(newKw.uniqueNodeKey());
    }

    const auto numOrig = this->m_keywords.size();
    this->m_keywords.insert(this->m_keywords.end(),
                            summaryNodes.begin(),
                            summaryNodes.end());

    std::inplace_merge(this->m_keywords.begin(),
                       this->m_keywords.begin() + numOrig,
                       this->m_keywords.end());

    return summaryNodes;
}

bool SummaryConfig::hasKeyword( const std::string& keyword ) const {
    return short_keywords.find(keyword) != short_keywords.end();
}


bool SummaryConfig::hasSummaryKey(const std::string& keyword ) const {
    return summary_keywords.find(keyword) != summary_keywords.end();
}

const SummaryConfigNode& SummaryConfig::operator[](std::size_t index) const {
    return this->m_keywords[index];
}


bool SummaryConfig::match(const std::string& keywordPattern) const
{
    return std::any_of(this->short_keywords.begin(), this->short_keywords.end(),
                        [&keywordPattern](const auto& keyword)
                        {
                            return shmatch(keywordPattern, keyword);
                        });
}

SummaryConfig::keyword_list
SummaryConfig::keywords(const std::string& keywordPattern) const
{
    auto kw_list = keyword_list{};

    std::copy_if(this->m_keywords.begin(), this->m_keywords.end(),
                 std::back_inserter(kw_list),
                 [&keywordPattern](const auto& kw)
                 { return shmatch(keywordPattern, kw.keyword()); });

    return kw_list;
}


size_t SummaryConfig::size() const {
    return this->m_keywords.size();
}

// Can be used to query if a certain 3D field, e.g. PRESSURE, is required to
// calculate a summary variable.
bool SummaryConfig::require3DField(const std::string& keyword) const
{
    auto iter = required_fields.find(keyword);
    if (iter == required_fields.end()) {
        return false;
    }

    return std::any_of(iter->second.begin(), iter->second.end(),
                       [this](const std::string& smryKw)
                       { return this->match(smryKw); });
}


std::unordered_set<std::string> SummaryConfig::wbp_wells() const {
    std::unordered_set<std::string> wells;
    for (const auto& node : this->keywords("WBP*"))
        wells.insert( node.namedEntity() );
    return wells;
}


std::set<std::string> SummaryConfig::fip_regions() const {
    std::set<std::string> reg_set;
    for (const auto& node : this->m_keywords) {
        if (node.category() == EclIO::SummaryNode::Category::Region)
            reg_set.insert( node.fip_region() );
    }
    return reg_set;
}

std::set<std::string> SummaryConfig::fip_regions_interreg_flow() const
{
    using Category = EclIO::SummaryNode::Category;

    auto reg_set = std::set<std::string>{};

    for (const auto& node : this->m_keywords) {
        if ((node.category() == Category::Region) &&
            is_region_to_region(node.keyword()))
        {
            reg_set.insert(node.fip_region());
        }
    }

    return reg_set;
}

bool SummaryConfig::operator==(const Opm::SummaryConfig& data) const {
    return this->m_keywords == data.m_keywords &&
           this->short_keywords == data.short_keywords &&
           this->summary_keywords == data.summary_keywords;
}

void SummaryConfig::handleProcessingInstruction(const std::string& keyword) {
    if (keyword == "RUNSUM") {
        runSummaryConfig.create = true;
    } else if (keyword == "NARROW") {
        runSummaryConfig.narrow = true;
    } else if (keyword == "SEPARATE") {
        runSummaryConfig.separate = true;
    }
}

}
