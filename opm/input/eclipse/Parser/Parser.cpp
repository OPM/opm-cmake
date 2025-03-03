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

#include <opm/input/eclipse/Parser/Parser.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/OpmLog/LogUtil.hpp>
#include <opm/common/utility/OpmInputError.hpp>

#include <opm/input/eclipse/Parser/ErrorGuard.hpp>
#include <opm/input/eclipse/Parser/ParseContext.hpp>
#include <opm/input/eclipse/Parser/ParserItem.hpp>
#include <opm/input/eclipse/Parser/ParserKeyword.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/I.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/P.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/R.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/S.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/T.hpp>
#include <opm/input/eclipse/Parser/ParserRecord.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Deck/DeckItem.hpp>
#include <opm/input/eclipse/Deck/DeckKeyword.hpp>
#include <opm/input/eclipse/Deck/DeckRecord.hpp>
#include <opm/input/eclipse/Deck/DeckSection.hpp>
#include <opm/input/eclipse/Deck/ImportContainer.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FieldProps.hpp>

#include <opm/input/eclipse/Python/Python.hpp>

#include <opm/json/JsonObject.hpp>

#include <opm/common/utility/String.hpp>

#include "raw/RawConsts.hpp"
#include "raw/RawEnums.hpp"
#include "raw/RawKeyword.hpp"
#include "raw/RawRecord.hpp"
#include "raw/StarToken.hpp"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <iterator>
#include <optional>
#include <stack>
#include <stdexcept>
#include <string>
#include <string_view>
#include <regex>
#include <utility>
#include <vector>

#include <fmt/format.h>

namespace {
    /// \brief Whether a keyword is a global keyword.
    ///
    /// Those are allowed before RUNSPEC.
    /// \return whether keyword is global (ECHO, NOECHO, INCLUDE,
    ///         COLUMNS, FORMFEED, SKIP, SKIP100, SKIP300, END)
    bool isGlobalKeyword(const Opm::DeckKeyword& keyword )
    {
        const auto& name = keyword.name();
        using namespace std::string_view_literals;
        const auto kw_list = std::array {
            "ECHO"sv,
            "NOECHO"sv,
            "INCLUDE"sv,
            "COLUMNS"sv,
            "FORMFEED"sv,
            "SKIP"sv,
            "ENDSKIP"sv,
            "SKIP100"sv,
            "SKIP300"sv
        };
        return std::find(kw_list.begin(), kw_list.end(), name) != kw_list.end();
    }

    // If ROCKOPTS does NOT exist, then the number of records is NTPVT (= TABDIMS(2))
    //
    // Otherwise, the number of records depends on ROCKOPTS(3), (= TABLE_TYPE)
    //
    //  1) ROCKOPTS(3) == "SATNUM"  => NTSFUN (= TABDIMS( 1))
    //  2) ROCKOPTS(3) == "ROCKNUM" => NTROCC (= TABDIMS(13))
    //  3) ROCKOPTS(3) == "PVTNUM"  => NTPVT  (= TABDIMS( 2)) (default setting)
    //
    // Finally, if ROCKOPTS(3) == "ROCKNUM", but NTROCC is defaulted, then
    // the number of records is NTPVT.
    std::size_t targetSizeRockFromTabdims(const Opm::Deck& deck)
    {
        const auto& tabd = deck.get<Opm::ParserKeywords::TABDIMS>().back().getRecord(0);
        const auto ntPvt = tabd.getItem<Opm::ParserKeywords::TABDIMS::NTPVT>().get<int>(0);

        if (! deck.hasKeyword<Opm::ParserKeywords::ROCKOPTS>()) {
            return ntPvt;
        }

        const auto& tableTypeItem = deck.get<Opm::ParserKeywords::ROCKOPTS>().back()
            .getRecord(0).getItem<Opm::ParserKeywords::ROCKOPTS::TABLE_TYPE>();

        if (tableTypeItem.defaultApplied(0)) {
            return ntPvt;
        }

        const auto tableType = tableTypeItem.getTrimmedString(0);

        if (tableType == Opm::ParserKeywords::PVTNUM::keywordName) {
            return ntPvt;
        }

        if (tableType == Opm::ParserKeywords::SATNUM::keywordName) {
            return tabd.getItem<Opm::ParserKeywords::TABDIMS::NTSFUN>().get<int>(0);
        }

        if (tableType == Opm::ParserKeywords::ROCKNUM::keywordName) {
            const auto& ntrocc = tabd.getItem<Opm::ParserKeywords::TABDIMS::NTROCC>();

            return ntrocc.defaultApplied(0) ? ntPvt : ntrocc.get<int>(0);
        }

        throw std::invalid_argument {
            fmt::format("Unknown {} table type \"{}\"",
                        Opm::ParserKeywords::ROCKOPTS::keywordName, tableType)
        };
    }

    std::size_t defaultTargetSizeRock()
    {
        // No TABDIMS => NTSFUN == NTPVT == NTROCC
        return Opm::ParserKeywords::TABDIMS::NTPVT::defaultValue;
    }

    std::size_t targetSizeRock(const Opm::Deck& deck)
    {
        return deck.hasKeyword<Opm::ParserKeywords::TABDIMS>()
            ? targetSizeRockFromTabdims(deck)
            : defaultTargetSizeRock();
    }
}

namespace Opm {

namespace {
namespace str {
const std::string emptystr = "";

struct find_comment {
    /*
     * A note on performance: using a function to plug functionality into
     * find_terminator rather than plain functions because it almost ensures
     * inlining, where the plain function can reduce to a function pointer.
     */
    template< typename Itr >
    Itr operator()( Itr begin, Itr end ) const {
        auto itr = std::find( begin, end, '-' );
        for( ; itr != end; itr = std::find( itr + 1, end, '-' ) )
            if( (itr + 1) != end &&  *( itr + 1 ) == '-' ) return itr;

        return end;
    }
};

template< typename Itr, typename Term >
inline Itr find_terminator( Itr begin, Itr end, Term terminator ) {

    auto pos = terminator( begin, end );

    if( pos == begin || pos == end) return pos;

    auto qbegin = std::find_if( begin, end, RawConsts::is_quote() );

    if( qbegin == end || qbegin > pos )
        return pos;

    auto qend = std::find( qbegin + 1, end, *qbegin );

    // Quotes are not balanced - probably an error?!
    if( qend == end ) return end;

    return find_terminator( qend + 1, end, terminator );
}

/**
    This function will return a copy of the input string where all
    characters following '--' are removed. The copy is a view and relies on
    the source string to remain alive. The function handles quoting with
    single quotes and double quotes:

    ABC --Comment                =>  ABC
    ABC '--Comment1' --Comment2  =>  ABC '--Comment1'
    ABC "-- Not balanced quote?  =>  ABC "-- Not balanced quote?
*/
static inline std::string_view strip_comments( std::string_view str ) {
    auto terminator = find_terminator( str.begin(), str.end(), find_comment() );
    std::size_t size = std::distance(str.begin(), terminator);
    return str.substr(0, size);
}

template< typename Itr >
inline Itr trim_left( Itr begin, Itr end ) {
    return std::find_if_not( begin, end, RawConsts::is_separator() );
}

template< typename Itr >
inline Itr trim_right( Itr begin, Itr end ) {

    std::reverse_iterator< Itr > rbegin( end );
    std::reverse_iterator< Itr > rend( begin );

    return std::find_if_not( rbegin, rend, RawConsts::is_separator() ).base();
}

inline std::string_view trim( std::string_view str ) {
    auto fst = trim_left( str.begin(), str.end() );
    auto lst = trim_right( fst, str.end() );
    const std::size_t start = std::distance(str.begin(), fst);
    const std::size_t size = std::distance(fst, lst);
    return str.substr(start, size);
}

inline std::string_view del_after_first_slash( std::string_view view ) {
    using itr = std::string_view::const_iterator;
    const auto term = []( itr begin, itr end ) {
        return std::find( begin, end, '/' );
    };

    auto begin = view.begin();
    auto end = view.end();
    auto slash = find_terminator( begin, end, term );

    /* we want to preserve terminating slashes */
    if (slash != end) {
        ++slash;
    }
    const std::size_t size = std::distance(begin, slash);
    return view.substr(0, size);
}

inline std::string_view del_after_last_slash( std::string_view view ) {
    if (view.empty()) {
        return view;
    }
    auto slash = view.find_last_of('/');
    if (slash == std::string_view::npos) {
        slash = view.size();
    }
    else {
        ++slash;
    }

    return view.substr(0, slash);
}

inline std::string_view del_after_slash(std::string_view view, bool raw_strings) {
    if (raw_strings)
        return del_after_last_slash(view);
    else
        return del_after_first_slash(view);
}


inline bool getline( std::string_view& input, std::string_view& line ) {
    if (input.empty()) {
        return false;
    }

    /* we know that we always append a newline onto the input string, so we can
     * safely assume that pos+1 will either be end-of-input (i.e. empty range)
     * or the start of the next line
     */

    const auto pos = input.find_first_of('\n');
    line = input.substr(0, pos);
    input = input.substr(pos+1);

    return true;
}

/*
 * Read the input file and remove everything that isn't interesting data,
 * including stripping comments, removing leading/trailing whitespaces and
 * everything after (terminating) slashes. Manually copying into the string for
 * performance.
 */
inline std::string fast_clean( const std::string& str ) {
    std::string dst;
    dst.resize( str.size() );

    std::string_view input( str ), line;
    auto dsti = dst.begin();
    while( true ) {

        if ( getline( input, line ) ) {
            line = trim( strip_comments(line));

            std::copy( line.begin(), line.end(), dsti );
            dsti += std::distance( line.begin(), line.end() );
            *dsti++ = '\n';
        } else
            break;
    }

    dst.resize( std::distance( dst.begin(), dsti ) );
    return dst;
}

inline bool starts_with(const std::string_view& view, const std::string& str) {
    auto str_size = str.size();
    if (str_size > view.size())
        return false;

    auto str_data = str.data();
    auto pos = view.begin();

    std::size_t si = 0;
    while (true) {
        if (*pos != str_data[si])
            return false;

        ++pos;
        ++si;

        if (si == str_size)
            return true;
    }
}

inline std::string clean( const std::vector<std::pair<std::string, std::string>>& code_keywords, const std::string& str ) {
    auto count = std::count_if(code_keywords.begin(), code_keywords.end(), [&str](const std::pair<std::string, std::string>& code_pair)
                                                                  {
                                                                     return str.find(code_pair.first) != std::string::npos;
                                                                   });

    if (count == 0)
        return fast_clean(str);
    else {
        std::string dst;
        dst.resize( str.size() );

        std::string_view input( str ), line;
        auto dsti = dst.begin();
        while( true ) {
            for (const auto& code_pair : code_keywords) {
                const auto& keyword = code_pair.first;

                if (starts_with(input, keyword)) {
                    std::string end_string = code_pair.second;
                    auto end_pos = input.find(end_string);
                    if (end_pos == std::string::npos) {
                        std::copy(input.begin(), input.end(), dsti);
                        dsti += std::distance( input.begin(), input.end() );
                        input = {};
                        break;
                    } else {
                        end_pos += end_string.size();
                        std::copy(input.begin(), input.begin() + end_pos, dsti);
                        dsti += end_pos;
                        *dsti++ = '\n';
                        input.remove_prefix(end_pos + 1);
                        break;
                    }
                }
            }

            if ( getline( input, line ) ) {
                line = trim( strip_comments(line));

                std::copy( line.begin(), line.end(), dsti );
                dsti += std::distance( line.begin(), line.end() );
                *dsti++ = '\n';
            } else
                break;
        }

        dst.resize( std::distance( dst.begin(), dsti ) );
        return dst;
    }
}




inline std::string make_deck_name(const std::string_view& str) {
    auto first_sep = std::find_if( str.begin(), str.end(), RawConsts::is_separator() );
    return uppercase( std::string( str.substr( 0, first_sep - str.begin()) ));
}


inline std::string_view update_record_buffer(const std::string_view& record_buffer,
                                             const std::string_view& line) {
    if (record_buffer.empty())
        return line;
    else {
        const std::size_t size = line.data() + line.size() - record_buffer.data();
        // intentionally not using substr since that will clamp the size
        return {record_buffer.data(), size};
    }
}


inline bool isTerminator(const std::string_view& line) {
    return (line.size() == 1 && line.back() == RawConsts::slash);
}

inline bool isTerminatedRecordString(const std::string_view& line) {
    return (line.back() == RawConsts::slash);
}

}

struct file {
    file( std::filesystem::path p, const std::string& in ) :
        input( in ), path( p )
    {}

    std::string_view input;
    size_t lineNR = 0;
    std::filesystem::path path;
};


class InputStack : public std::stack< file, std::vector< file > > {
    public:
        void push( std::string&& input, std::filesystem::path p = "<memory string>" );

    private:
        std::list< std::string > string_storage;
        using base = std::stack< file, std::vector< file > >;
};

void InputStack::push( std::string&& input, std::filesystem::path p ) {
    this->string_storage.push_back( std::move( input ) );
    this->emplace( p, this->string_storage.back() );
}

class ParserState {
    public:
        ParserState( const std::vector<std::pair<std::string,std::string>>&,
                     const ParseContext&, ErrorGuard&,
                     const std::set<Opm::Ecl::SectionType>& ignore = {});

        ParserState( const std::vector<std::pair<std::string,std::string>>&,
                     const ParseContext&, ErrorGuard&,
                     std::filesystem::path, const std::set<Opm::Ecl::SectionType>& ignore = {});

        void loadString( const std::string& );
        void loadFile( const std::filesystem::path& );
        void openRootFile( const std::filesystem::path& );

        void handleRandomText(const std::string_view& ) const;
        std::optional<std::filesystem::path> getIncludeFilePath( std::string ) const;
        void addPathAlias( const std::string& alias, const std::string& path );

        const std::filesystem::path& current_path() const;
        size_t line() const;

        bool done() const;
        std::string_view getline();
        void ungetline(const std::string_view& ln);
        void closeFile();

        const std::set<Opm::Ecl::SectionType>& get_ignore() {return ignore_sections; };
        bool check_section_keywords(bool& has_edit, bool& has_regions, bool& has_summary);

    private:
        const std::vector<std::pair<std::string, std::string>> code_keywords;
        InputStack input_stack;

        std::set<Opm::Ecl::SectionType> ignore_sections;
        std::map< std::string, std::string > pathMap;

    public:
        ParserKeywordSizeEnum lastSizeType = SLASH_TERMINATED;
        std::string lastKeyWord;

        Deck deck;
        std::filesystem::path rootPath;
        std::unique_ptr<Python> python;
        const ParseContext& parseContext;
        ErrorGuard& errors;
        bool unknown_keyword = false;
};

const std::filesystem::path& ParserState::current_path() const {
    return this->input_stack.top().path;
}

size_t ParserState::line() const {
    return this->input_stack.top().lineNR;
}

bool ParserState::done() const {

    while( !this->input_stack.empty() &&
            this->input_stack.top().input.empty() )
        const_cast< ParserState* >( this )->input_stack.pop();

    return this->input_stack.empty();
}

std::string_view ParserState::getline() {
    std::string_view ln;

    str::getline( this->input_stack.top().input, ln );
    this->input_stack.top().lineNR++;

    return ln;
}



void ParserState::ungetline(const std::string_view& line) {
    auto& file_view = this->input_stack.top().input;
    if (line.data() + line.size() + 1 != file_view.data())
        throw std::invalid_argument("line view does not immediately proceed file_view");

    const auto* end = file_view.data() + file_view.size();
    const std::size_t size = std::distance(line.data(), end);
    // intentionally not using substr since that will clamp the size
    file_view = {line.data(), size};
    this->input_stack.top().lineNR--;
}




void ParserState::closeFile() {
    this->input_stack.pop();
}

ParserState::ParserState(const std::vector<std::pair<std::string, std::string>>& code_keywords_arg,
                         const ParseContext& __parseContext,
                         ErrorGuard& errors_arg,
                         const std::set<Opm::Ecl::SectionType>& ignore) :
    code_keywords(code_keywords_arg),
    ignore_sections(ignore),
    python( std::make_unique<Python>() ),
    parseContext( __parseContext ),
    errors( errors_arg )
{}

ParserState::ParserState( const std::vector<std::pair<std::string, std::string>>& code_keywords_arg,
                          const ParseContext& context,
                          ErrorGuard& errors_arg,
                          std::filesystem::path p,
                          const std::set<Opm::Ecl::SectionType>& ignore ) :
    code_keywords(code_keywords_arg),
    ignore_sections(ignore),
    rootPath( std::filesystem::canonical( p ).parent_path() ),
    python( std::make_unique<Python>() ),
    parseContext( context ),
    errors( errors_arg )
{
    openRootFile( p );
}

bool ParserState::check_section_keywords(bool& has_edit, bool& has_regions, bool& has_summary) {

    std::string_view root_file_str = this->input_stack.top().input;

    has_edit = false;
    has_regions = false;
    has_summary = false;

    int n = 0;
    auto p0 = root_file_str.find_first_not_of(" \t\n");

    while (p0 != std::string::npos){

        auto p1 = root_file_str.find_first_of(" \t\n", p0 + 1);

        if (root_file_str.substr(p0, p1-p0) == "RUNSPEC")
            n++;
        else if (root_file_str.substr(p0, p1-p0) == "GRID")
            n++;
        else if (root_file_str.substr(p0, p1-p0) == "EDIT")
            has_edit = true;
        else if (root_file_str.substr(p0, p1-p0) == "PROPS")
            n++;
        else if (root_file_str.substr(p0, p1-p0) == "REGIONS")
            has_regions = true;
        else if (root_file_str.substr(p0, p1-p0) == "SOLUTION")
            n++;
        else if (root_file_str.substr(p0, p1-p0) == "SUMMARY")
            has_summary = true;
        else if (root_file_str.substr(p0, p1-p0) == "SCHEDULE")
            n++;

        p0 = root_file_str.find_first_not_of(" \t\n", p1);
    }

    if (n < 5)
        return false;
    else
        return true;
}

void ParserState::loadString(const std::string& input) {
    this->input_stack.push( str::clean( this->code_keywords, input + "\n" ) );
}

void ParserState::loadFile(const std::filesystem::path& inputFile) {

    const auto closer = []( std::FILE* f ) { std::fclose( f ); };
    std::unique_ptr<std::FILE, decltype(closer)> ufp{
        std::fopen( inputFile.generic_string().c_str(), "rb" ),
        closer
    };

    // make sure the file we'd like to parse is readable
    if( !ufp ) {
        std::string msg = "Could not read from file: " + inputFile.string();
        parseContext.handleError( ParseContext::PARSE_MISSING_INCLUDE , msg, {}, errors);
        return;
    }

    /*
     * read the input file C-style. This is done for performance
     * reasons, as streams are slow
     */

    auto* fp = ufp.get();
    std::string buffer;
    std::fseek( fp, 0, SEEK_END );
    buffer.resize( std::ftell( fp ) + 1 );
    std::rewind( fp );
    const auto readc = std::fread( &buffer[ 0 ], 1, buffer.size() - 1, fp );
    buffer.back() = '\n';

    if( std::ferror( fp ) || readc != buffer.size() - 1 )
        throw std::runtime_error( "Error when reading input file '"
                                  + inputFile.string() + "'" );

    this->input_stack.push( str::clean( this->code_keywords, buffer ), inputFile );
}

/*
 * We have encountered 'random' characters in the input file which
 * are not correctly formatted as a keyword heading, and not part
 * of the data section of any keyword.
 */

void ParserState::handleRandomText(const std::string_view& keywordString) const
{
    const std::string trimmedCopy { keywordString };
    const KeywordLocation location {
        lastKeyWord, this->current_path().generic_string(), this->line()
    };

    std::string errorKey{};
    std::string msg{};
    if (trimmedCopy[0] == '/') {
        errorKey = ParseContext::PARSE_RANDOM_SLASH;
        msg = "Extra line starting with '/' detected in {file} line {line}";
    }
    else if (lastSizeType == OTHER_KEYWORD_IN_DECK) {
        errorKey = ParseContext::PARSE_EXTRA_RECORDS;
        msg = "Too many records in keyword {keyword}\n"
            "In {file} line {line}";
    }
    else {
        errorKey = ParseContext::PARSE_RANDOM_TEXT;
        msg = fmt::format("String {} not formatted as valid keyword\n"
                          "In {{file}} line {{line}}.", keywordString);
    }

    parseContext.handleError(errorKey, msg, location, errors);
}


void ParserState::openRootFile( const std::filesystem::path& inputFile) {

    this->loadFile( inputFile );
    this->deck.setDataFile( inputFile.string() );
    const std::filesystem::path& inputFileCanonical = std::filesystem::canonical(inputFile);
    this->rootPath = inputFileCanonical.parent_path();
}

std::optional<std::filesystem::path> ParserState::getIncludeFilePath( std::string path ) const {
    static const std::string pathKeywordPrefix("$");
    static const std::string validPathNameCharacters("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_");

    size_t positionOfPathName = path.find(pathKeywordPrefix);

    if ( positionOfPathName != std::string::npos) {
        std::string stringStartingAtPathName = path.substr(positionOfPathName+1);
        size_t cutOffPosition = stringStartingAtPathName.find_first_not_of(validPathNameCharacters);
        std::string stringToFind = stringStartingAtPathName.substr(0, cutOffPosition);
        std::string stringToReplace = this->pathMap.at( stringToFind );
        replaceAll(path, pathKeywordPrefix + stringToFind, stringToReplace);
    }

    // Check if there are any backslashes in the path...
    if (path.find('\\') != std::string::npos) {
        // ... if so, replace with slashes and create a warning.
        std::replace(path.begin(), path.end(), '\\', '/');
        OpmLog::warning("Replaced one or more backslash with a slash in an INCLUDE path.");
    }

    // trim leading and trailing whitespace just like the other simulator
    std::regex trim_regex("^\\s+|\\s+$");
    const auto& trimmed_path = std::regex_replace(path, trim_regex, "");
    std::filesystem::path includeFilePath(trimmed_path);

    if (includeFilePath.is_relative())
        includeFilePath = this->rootPath / includeFilePath;

    try {
        includeFilePath = std::filesystem::canonical(includeFilePath);
    } catch (const std::filesystem::filesystem_error& fs_error) {
        parseContext.handleError( ParseContext::PARSE_MISSING_INCLUDE ,
                                  fmt::format("File '{}' included via INCLUDE"
                                              " directive does not exist.",
                                              trimmed_path),
                                  {}, errors);
        return {};
    }

    return includeFilePath;
}

void ParserState::addPathAlias( const std::string& alias, const std::string& path ) {
    this->pathMap.emplace( alias, path );
}

RawKeyword*
newRawKeyword(const ParserKeyword& parserKeyword,
              const std::string&   keywordString,
              ParserState&         parserState,
              const Parser&        parser)
{
    for (const auto& keyword : parserKeyword.prohibitedKeywords()) {
        if (parserState.deck.hasKeyword(keyword)) {
            parserState
                .parseContext
                .handleError(
                    ParseContext::PARSE_INVALID_KEYWORD_COMBINATION,
                    fmt::format("Incompatible keyword combination: {} declared "
                                "when {} is already present.", keywordString, keyword),
                    KeywordLocation{
                        keywordString,
                        parserState.current_path().generic_string(),
                        parserState.line()
                    },
                    parserState.errors
                );
        }
    }

    for (const auto& keyword : parserKeyword.requiredKeywords()) {
        if (!parserState.deck.hasKeyword(keyword)) {
            parserState
                .parseContext
                .handleError(
                    ParseContext::PARSE_INVALID_KEYWORD_COMBINATION,
                    fmt::format("Incompatible keyword combination: {} declared, "
                                "but {} is missing.", keywordString, keyword),
                    KeywordLocation{
                        keywordString,
                        parserState.current_path().generic_string(),
                        parserState.line()
                    },
                    parserState.errors
                );
        }
    }

    const bool raw_string_keyword = parserKeyword.rawStringKeyword();

    if ((parserKeyword.getSizeType() == SLASH_TERMINATED) ||
        (parserKeyword.getSizeType() == UNKNOWN) ||
        (parserKeyword.getSizeType() == DOUBLE_SLASH_TERMINATED))
    {
        const auto size_type = parserKeyword.getSizeType();
        Raw::KeywordSizeEnum rawSizeType;

        switch (size_type) {
            case SLASH_TERMINATED:        rawSizeType = Raw::SLASH_TERMINATED; break;
            case UNKNOWN:                 rawSizeType = Raw::UNKNOWN; break;
            case DOUBLE_SLASH_TERMINATED: rawSizeType = Raw::DOUBLE_SLASH_TERMINATED; break;
            default:
                throw std::logic_error("Should not be here!");
        }

        return new RawKeyword(keywordString,
                              parserState.current_path().string(),
                              parserState.line(),
                              raw_string_keyword,
                              rawSizeType);
    }

    if (parserKeyword.getSizeType() == SPECIAL_CASE_ROCK) {
        if (parserKeyword.getName() != ParserKeywords::ROCK::keywordName) {
            throw std::logic_error {
                fmt::format("Special case size handling for ROCK "
                            "cannot be applied to keyword {}",
                            parserKeyword.getName())
            };
        }

        return new RawKeyword {
            keywordString,
            parserState.current_path().string(),
            parserState.line(),
            raw_string_keyword,
            Raw::FIXED,
            parserKeyword.min_size(),
            targetSizeRock(parserState.deck)
        };
    }

    if (parserKeyword.hasFixedSize()) {
        auto size_type = Raw::FIXED;
        if (parserKeyword.isCodeKeyword()) {
            size_type = Raw::CODE;
        }

        return new RawKeyword(keywordString,
                              parserState.current_path().string(),
                              parserState.line(),
                              raw_string_keyword,
                              size_type,
                              parserKeyword.min_size(),
                              parserKeyword.getFixedSize());
    }

    const auto& keyword_size = parserKeyword.getKeywordSize();
    const auto& deck = parserState.deck;
    const auto size_type = parserKeyword.isTableCollection()
        ? Raw::TABLE_COLLECTION : Raw::FIXED;

    if (deck.hasKeyword(keyword_size.keyword())) {
        const auto& sizeDefinitionKeyword = deck[keyword_size.keyword()].back();
        const auto& record = sizeDefinitionKeyword.getRecord(0);

        auto targetSize = record.getItem(keyword_size.item()).get<int>(0) + keyword_size.size_shift();
        if (parserKeyword.isAlternatingKeyword()) {
            targetSize *= std::distance(parserKeyword.begin(), parserKeyword.end());
        }

        return new RawKeyword(keywordString,
                              parserState.current_path().string(),
                              parserState.line(),
                              raw_string_keyword,
                              size_type,
                              parserKeyword.min_size(),
                              targetSize);
    }

    const auto msg_fmt = fmt::format("Problem with {{keyword}} - missing {0}\n"
                                     "In {{file}} line {{line}}\n"
                                     "For the keyword {{keyword}} we expect "
                                     "to read the number of records from "
                                     "keyword {0}, {0} was not found",
                                     keyword_size.keyword());

    parserState.parseContext.handleError(ParseContext::PARSE_MISSING_DIMS_KEYWORD,
                                         msg_fmt,
                                         KeywordLocation {
                                             keywordString,
                                             parserState.current_path().string(),
                                             parserState.line()
                                         },
                                         parserState.errors);

    const auto& keyword  = parser.getKeyword(keyword_size.keyword());
    const auto& record   = keyword.getRecord(0);
    const auto& int_item = record.get(keyword_size.item());

    const auto targetSize = int_item.getDefault<int>() + keyword_size.size_shift();
    return new RawKeyword(keywordString,
                          parserState.current_path().string(),
                          parserState.line(),
                          raw_string_keyword,
                          size_type,
                          parserKeyword.min_size(),
                          targetSize);
}

RawKeyword*
newRawKeyword(const std::string&      deck_name,
              ParserState&            parserState,
              const Parser&           parser,
              const std::string_view& line)
{
    if (deck_name.size() > RawConsts::maxKeywordLength) {
        const std::string keyword8 = deck_name.substr(0, RawConsts::maxKeywordLength);
        if (parser.isRecognizedKeyword(keyword8)) {
            const auto msg = std::string {
                "Keyword {keyword} to long - only eight "
                "first characters recognized\n"
                "In {file} line {line}\n"
            };

            parserState.parseContext.handleError(ParseContext::PARSE_LONG_KEYWORD,
                                                 msg,
                                                 KeywordLocation {
                                                     deck_name,
                                                     parserState.current_path().string(),
                                                     parserState.line()
                                                 },
                                                 parserState.errors);

            parserState.unknown_keyword = false;
            const auto& parserKeyword = parser.getParserKeywordFromDeckName(keyword8);

            return newRawKeyword(parserKeyword, keyword8, parserState, parser);
        }
        else if (parser.isBaseRecognizedKeyword(deck_name)) {
            // Typically an OPM extended keyword such as STRESSEQUILNUM.
            parserState.unknown_keyword = false;
            const auto& parserKeyword = parser.getParserKeywordFromDeckName(deck_name);

            return newRawKeyword(parserKeyword, deck_name, parserState, parser);
        }
        else {
            parserState.parseContext.handleUnknownKeyword(deck_name,
                                                          KeywordLocation {
                                                              deck_name,
                                                              parserState.current_path().string(),
                                                              parserState.line()
                                                          }, parserState.errors);
            parserState.unknown_keyword = true;

            return nullptr;
        }
    }

    if (parser.isRecognizedKeyword(deck_name)) {
        parserState.unknown_keyword = false;
        const auto& parserKeyword = parser.getParserKeywordFromDeckName(deck_name);

        return newRawKeyword(parserKeyword, deck_name, parserState, parser);
    }

    if (ParserKeyword::validDeckName(deck_name)) {
        parserState.parseContext.handleUnknownKeyword(deck_name,
                                                      KeywordLocation{
                                                          deck_name,
                                                          parserState.current_path().string(),
                                                          parserState.line()},
                                                      parserState.errors);
        parserState.unknown_keyword = true;

        return nullptr;
    }

    if (! parserState.unknown_keyword) {
        parserState.handleRandomText(line);
    }

    return nullptr;
}

std::unique_ptr<RawKeyword> tryParseKeyword( ParserState& parserState, const Parser& parser) {
    bool is_title = false;
    bool skip = false;
    std::unique_ptr<RawKeyword> rawKeyword;
    std::string_view record_buffer(str::emptystr);
    std::optional<ParserKeyword> parserKeyword;
    const bool silent = parser.silent();
    while( !parserState.done() ) {
        auto line = parserState.getline();

        if( line.empty() && !rawKeyword ) continue;
        if( line.empty() && !is_title ) continue;

        std::string deck_name = str::make_deck_name( line );
        if (parserState.parseContext.isActiveSkipKeyword(deck_name)) {
            skip = true;
            auto msg = fmt::format("{:5} Reading {:<8} in {} line {} \n      ... ignoring everything until 'ENDSKIP' ... ", "", "SKIP", parserState.current_path().string(), parserState.line());
            if (!silent) OpmLog::info(msg);
        } else if (deck_name == "ENDSKIP") {
            skip = false;
            auto msg = fmt::format("{:5} Reading {:<8} in {} line {}", "", "ENDSKIP", parserState.current_path().string(), parserState.line());
            if (!silent) OpmLog::info(msg);
            continue;
        }
        if (skip) continue;

        if( !rawKeyword ) {
            /*
              Extracting a possible keywordname from a line of deck input
              involves several steps.

              1. The make_deck_name() function will strip off everyhing
                 following the first white-space separator and uppercase the
                 string.

              2. The ParserKeyword::validDeckName() verifies that the keyword
                 candidate only contains valid characters.

              3. In the newRawKeyword() function the first 8 characters of the
                 deck_name is used to look for the keyword in the Parser
                 container.
            */
            if (ParserKeyword::validDeckName(deck_name)) {
                const auto ptr = newRawKeyword( deck_name, parserState, parser, line );
                if (ptr) {
                    rawKeyword.reset( ptr );
                    parserKeyword = parser.getParserKeywordFromDeckName(rawKeyword->getKeywordName());
                    parserState.lastSizeType = parserKeyword->getSizeType();
                    parserState.lastKeyWord = deck_name;
                    if (rawKeyword->isFinished())
                        return rawKeyword;

                    if (deck_name == "TITLE")
                        is_title = true;
                }
            } else {
                /* We are looking at some random gibberish?! */
                if (!parserState.unknown_keyword)
                    parserState.handleRandomText( line );
            }
        } else {
            if (rawKeyword->getSizeType() == Raw::CODE) {
                auto end_pos = line.find(parserKeyword->codeEnd());
                if (end_pos != std::string::npos) {
                    std::string_view line_content = line.substr(0, end_pos);
                    record_buffer = str::update_record_buffer(record_buffer, line_content);

                    RawRecord record(record_buffer, rawKeyword->location(), true);
                    rawKeyword->addRecord(record);
                    return rawKeyword;
                } else
                    record_buffer = str::update_record_buffer(record_buffer, line);

                continue;
            }

            if (rawKeyword->can_complete()) {
                /*
                  When we are spinning through a keyword of size type UNKNOWN it
                  is essential to recognize a string as the next keyword. The
                  line starting a new keyword can have arbitrary rubbish
                  following the keyword name - i.e. this  line

                    PORO  Here comes some random gibberish which should be ignored
                       10000*0.15 /

                  To ensure the keyword 'PORO' is recognized in the example
                  above we remove everything following the first space in the
                  line variable before we check if it is the start of a new
                  keyword.
                */
                if( parser.isRecognizedKeyword( deck_name ) ) {
                    rawKeyword->terminateKeyword();
                    parserState.ungetline(line);
                    return rawKeyword;
                }
            }

            line = str::del_after_slash(line, rawKeyword->rawStringKeyword());
            record_buffer = str::update_record_buffer(record_buffer, line);
            if (is_title) {
                if (record_buffer.empty()) {
                    RawRecord record("opm/flow simulation", rawKeyword->location());
                    rawKeyword->addRecord(record);
                } else {
                    RawRecord record(record_buffer, rawKeyword->location());
                    rawKeyword->addRecord(record);
                }
                return rawKeyword;
            }


            if (str::isTerminator(record_buffer)) {
                if (rawKeyword->terminateKeyword())
                    return rawKeyword;
            }


            if (str::isTerminatedRecordString(record_buffer)) {
                const std::size_t size = record_buffer.size() - 1;
                RawRecord record(record_buffer.substr(0, size), rawKeyword->location());
                if (rawKeyword->addRecord(record))
                    return rawKeyword;

                record_buffer = str::emptystr;
            }
        }
    }

    if (rawKeyword) {
        if (rawKeyword->can_complete())
            rawKeyword->terminateKeyword();

        if (!rawKeyword->isFinished()) {
            throw OpmInputError {
                "Keyword is not properly terminated.",
                rawKeyword->location()
            };
        }
    }

    return rawKeyword;
}

std::string_view advance_parser_state( ParserState& parserState, const std::string& to_keyw )
{

    auto line = parserState.getline();

    while (line != to_keyw) {
        line = parserState.getline();
    }

    return line;
}


void addSectionKeyword(ParserState& parserState, const std::string& keyw)
{
    if (!parserState.deck.hasKeyword(keyw)){

        Opm::ParserKeyword section_keyw(keyw);
        Opm::DeckKeyword dk_keyw(section_keyw);

        parserState.deck.addKeyword(dk_keyw);
    }
}

void cleanup_deck_keyword_list(ParserState& parserState, const std::set<Opm::Ecl::SectionType>& ignore)
{
    bool ignore_runspec = ignore.find(Opm::Ecl::RUNSPEC) !=ignore.end()  ? true : false;
    bool ignore_grid = ignore.find(Opm::Ecl::GRID) !=ignore.end()  ? true : false;
    bool ignore_edit = ignore.find(Opm::Ecl::EDIT) !=ignore.end()  ? true : false;
    bool ignore_props = ignore.find(Opm::Ecl::PROPS) !=ignore.end()  ? true : false;
    bool ignore_regions = ignore.find(Opm::Ecl::REGIONS) !=ignore.end()  ? true : false;
    bool ignore_solution = ignore.find(Opm::Ecl::SOLUTION) !=ignore.end()  ? true : false;
    bool ignore_summary = ignore.find(Opm::Ecl::SUMMARY) !=ignore.end()  ? true : false;
    bool ignore_schedule = ignore.find(Opm::Ecl::SCHEDULE) !=ignore.end()  ? true : false;

    std::vector<std::string> keyw_names;
    keyw_names.reserve(parserState.deck.size());

    std::transform(parserState.deck.begin(), parserState.deck.end(),
                   std::back_inserter(keyw_names),
                   [](const auto& dk_keyw) { return dk_keyw.name(); });

    if (ignore_runspec){

        auto iter_from = std::find(keyw_names.begin(), keyw_names.end(), "RUNSPEC");
        auto iter_to = std::find(keyw_names.begin(), keyw_names.end(), "GRID");

        auto n1 = std::distance(keyw_names.begin(), iter_from);
        auto n2 = std::distance(keyw_names.begin(), iter_to);

        parserState.deck.remove_keywords(n1, n2);
        keyw_names.erase(keyw_names.begin() + n1, keyw_names.begin() + n2);
    }

    if (ignore_grid){

        auto iter_from = std::find(keyw_names.begin(), keyw_names.end(), "GRID");
        auto iter_to = std::find(keyw_names.begin(), keyw_names.end(), "EDIT");

        if (iter_to == keyw_names.end())
            iter_to = std::find(keyw_names.begin(), keyw_names.end(), "PROPS");

        auto n1 = std::distance(keyw_names.begin(), iter_from);
        auto n2 = std::distance(keyw_names.begin(), iter_to);

        parserState.deck.remove_keywords(n1, n2);
        keyw_names.erase(keyw_names.begin() + n1, keyw_names.begin() + n2);
    }

    if (ignore_edit){

        auto iter_from = std::find(keyw_names.begin(), keyw_names.end(), "EDIT");

        if (iter_from != keyw_names.end()){
            auto iter_to = std::find(keyw_names.begin(), keyw_names.end(), "PROPS");
            auto n1 = std::distance(keyw_names.begin(), iter_from);
            auto n2 = std::distance(keyw_names.begin(), iter_to);

            parserState.deck.remove_keywords(n1, n2);
            keyw_names.erase(keyw_names.begin() + n1, keyw_names.begin() + n2);

        }
    }

    if (ignore_props){

        auto iter_from = std::find(keyw_names.begin(), keyw_names.end(), "PROPS");
        auto iter_to = std::find(keyw_names.begin(), keyw_names.end(), "REGIONS");

        if (iter_to == keyw_names.end())
            iter_to = std::find(keyw_names.begin(), keyw_names.end(), "SOLUTION");

        auto n1 = std::distance(keyw_names.begin(), iter_from);
        auto n2 = std::distance(keyw_names.begin(), iter_to);

        parserState.deck.remove_keywords(n1, n2);
        keyw_names.erase(keyw_names.begin() + n1, keyw_names.begin() + n2);
    }

    if (ignore_regions){

        auto iter_from = std::find(keyw_names.begin(), keyw_names.end(), "REGIONS");

        if (iter_from != keyw_names.end()){
            auto iter_to = std::find(keyw_names.begin(), keyw_names.end(), "SOLUTION");

            auto n1 = std::distance(keyw_names.begin(), iter_from);
            auto n2 = std::distance(keyw_names.begin(), iter_to);

            parserState.deck.remove_keywords(n1, n2);
            keyw_names.erase(keyw_names.begin() + n1, keyw_names.begin() + n2);
        }
    }

    if (ignore_solution){

        auto iter_from = std::find(keyw_names.begin(), keyw_names.end(), "SOLUTION");
        auto iter_to = std::find(keyw_names.begin(), keyw_names.end(), "SUMMARY");

        if (iter_to == keyw_names.end())
            iter_to = std::find(keyw_names.begin(), keyw_names.end(), "SCHEDULE");

        auto n1 = std::distance(keyw_names.begin(), iter_from);
        auto n2 = std::distance(keyw_names.begin(), iter_to);

        parserState.deck.remove_keywords(n1, n2);
        keyw_names.erase(keyw_names.begin() + n1, keyw_names.begin() + n2);
    }

    if (ignore_summary){

        auto iter_from = std::find(keyw_names.begin(), keyw_names.end(), "SUMMARY");

        if (iter_from != keyw_names.end()){
            auto iter_to = std::find(keyw_names.begin(), keyw_names.end(), "SCHEDULE");
            auto n1 = std::distance(keyw_names.begin(), iter_from);
            auto n2 = std::distance(keyw_names.begin(), iter_to);

            parserState.deck.remove_keywords(n1, n2);
            keyw_names.erase(keyw_names.begin() + n1, keyw_names.begin() + n2);
        }
    }

    if (ignore_schedule){

        auto iter_from = std::find(keyw_names.begin(), keyw_names.end(), "SCHEDULE");

        auto n1 = std::distance(keyw_names.begin(), iter_from);
        auto n2 = std::distance(keyw_names.begin(), keyw_names.end());

        parserState.deck.remove_keywords(n1, n2);
        keyw_names.erase(keyw_names.begin() + n1, keyw_names.begin() + n2);
    }
}


bool parseState( ParserState& parserState, const Parser& parser ) {
    std::string filename = parserState.current_path().string();

    auto ignore = parserState.get_ignore();

    bool has_edit = true;
    bool has_regions = true;
    bool has_summary = true;

    if (ignore.size() > 0)
        if (!parserState.check_section_keywords(has_edit, has_regions, has_summary))
            throw std::runtime_error("Parsing individual sections not possible when section keywords in root input file");

    bool ignore_grid = ignore.find(Opm::Ecl::GRID) !=ignore.end()  ? true : false;
    bool ignore_edit = ignore.find(Opm::Ecl::EDIT) !=ignore.end()  ? true : false;
    bool ignore_props = ignore.find(Opm::Ecl::PROPS) !=ignore.end()  ? true : false;
    bool ignore_regions = ignore.find(Opm::Ecl::REGIONS) !=ignore.end()  ? true : false;
    bool ignore_solution = ignore.find(Opm::Ecl::SOLUTION) !=ignore.end()  ? true : false;
    bool ignore_summary = ignore.find(Opm::Ecl::SUMMARY) !=ignore.end()  ? true : false;
    bool ignore_schedule = ignore.find(Opm::Ecl::SCHEDULE) !=ignore.end()  ? true : false;

    if ((ignore_grid) && (!has_edit) && (!ignore_edit))
        ignore_grid = false;

    if ((ignore_props) && (!has_regions) && (!ignore_regions))
        ignore_props = false;

    if ((ignore_solution) && (!has_summary) && (!ignore_summary))
        ignore_solution = false;

    while( !parserState.done() ) {

        auto rawKeyword = tryParseKeyword( parserState, parser);
        bool do_not_add = false;

        if( !rawKeyword )
            continue;

        std::string_view keyw = rawKeyword->getKeywordName();
        if ((ignore_grid) && (keyw == "GRID")){

            do_not_add = true;
            addSectionKeyword(parserState, "GRID");

            if (has_edit){
                keyw = advance_parser_state( parserState, "EDIT" );
                addSectionKeyword(parserState, "EDIT");
            } else {
                keyw = advance_parser_state( parserState, "PROPS" );
                addSectionKeyword(parserState, "PROPS");
            }
        }

        if ((ignore_edit) && (keyw=="EDIT")){
            do_not_add = true;
            addSectionKeyword(parserState, "EDIT");
            keyw = advance_parser_state( parserState, "PROPS" );
            addSectionKeyword(parserState, "PROPS");
        }

        if ((ignore_props) && (keyw=="PROPS")){
            do_not_add = true;
            addSectionKeyword(parserState, "PROPS");

            if (has_regions){
                keyw = advance_parser_state( parserState, "REGIONS" );
                addSectionKeyword(parserState, "REGIONS");
            } else {
                keyw = advance_parser_state( parserState, "SOLUTION" );
                addSectionKeyword(parserState, "SOLUTION");
            }

        }

        if ((ignore_regions) && (keyw=="REGIONS")){
            do_not_add = true;
            addSectionKeyword(parserState, "REGIONS");
            keyw = advance_parser_state( parserState, "SOLUTION" );
            addSectionKeyword(parserState, "SOLUTION");
        }

        if ((ignore_solution) && (keyw=="SOLUTION")){
            do_not_add = true;
            addSectionKeyword(parserState, "SOLUTION");

            if (has_summary){
                keyw = advance_parser_state( parserState, "SUMMARY" );
                addSectionKeyword(parserState, "SUMMARY");
            } else {
                keyw = advance_parser_state( parserState, "SCHEDULE" );
                addSectionKeyword(parserState, "SCHEDULE");
            }
        }

        if ((ignore_summary) && (keyw=="SUMMARY")){
            do_not_add = true;
            addSectionKeyword(parserState, "SUMMARY");
            keyw = advance_parser_state( parserState, "SCHEDULE" );
            addSectionKeyword(parserState, "SCHEDULE");
        }

        if ((ignore_schedule) && (keyw=="SCHEDULE")){
            addSectionKeyword(parserState, "SCHEDULE");
            return true;
        }

        if (rawKeyword->getKeywordName() == Opm::RawConsts::end)
            return true;

        if (rawKeyword->getKeywordName() == Opm::RawConsts::endinclude) {
            parserState.closeFile();
            continue;
        }

        if (rawKeyword->getKeywordName() == Opm::RawConsts::paths) {
            for( const auto& record : *rawKeyword ) {
                std::string pathName = readValueToken<std::string>(record.getItem(0));
                std::string pathValue = readValueToken<std::string>(record.getItem(1));
                parserState.addPathAlias( pathName, pathValue );
            }

            continue;
        }

        if (rawKeyword->getKeywordName() == Opm::RawConsts::include) {
            auto& firstRecord = rawKeyword->getFirstRecord( );
            std::string includeFileAsString = readValueToken<std::string>(firstRecord.getItem(0));
            const auto& includeFile = parserState.getIncludeFilePath( includeFileAsString );

            if (includeFile.has_value()) {
                auto& deck_tree = parserState.deck.tree();
                deck_tree.add_include(std::filesystem::absolute(parserState.current_path()).generic_string(),
                                      includeFile.value().generic_string());
                parserState.loadFile(includeFile.value());
            }
            continue;
        }

        if( parser.isRecognizedKeyword( rawKeyword->getKeywordName() ) ) {
            const auto& kwname = rawKeyword->getKeywordName();
            const auto& parserKeyword = parser.getParserKeywordFromDeckName( kwname );
            {
                const auto& location = rawKeyword->location();
                auto msg = fmt::format("{:5} Reading {:<8} in {} line {}", parserState.deck.size(), rawKeyword->getKeywordName(), location.filename, location.lineno);
                if (!parser.silent()) OpmLog::info(msg);
            }
            try {
                if (rawKeyword->getKeywordName() ==  Opm::RawConsts::pyinput) {
                    if (parserState.python) {
                        std::string python_string = rawKeyword->getFirstRecord().getRecordString();
                        parserState.python->exec(python_string, parser, parserState.deck);
                    }
                    else
                        throw std::logic_error("Cannot yet embed Python while still running Python.");
                }
                else {
                    auto deck_keyword = parserKeyword.parse( parserState.parseContext,
                                                             parserState.errors,
                                                             *rawKeyword,
                                                             parserState.deck.getActiveUnitSystem(),
                                                             parserState.deck.getDefaultUnitSystem());

                    if (deck_keyword.name() == ParserKeywords::IMPORT::keywordName) {
                        bool formatted = deck_keyword.getRecord(0).getItem(1).get<std::string>(0)[0] == 'F';
                        const auto& import_file = parserState.getIncludeFilePath(deck_keyword.getRecord(0).getItem(0).getTrimmedString(0));

                        ImportContainer import(parser, parserState.deck.getActiveUnitSystem(), import_file.value().string(), formatted, parserState.deck.size());
                        for (auto kw : import)
                            parserState.deck.addKeyword(std::move(kw));
                    } else
                        if (!do_not_add)
                            parserState.deck.addKeyword( std::move(deck_keyword) );
                }
            } catch (const OpmInputError& opm_error) {
                throw;
            } catch (const std::exception& e) {
                /*
                  This catch-all of parsing errors is to be able to write a good
                  error message; the parser is quite confused at this state and
                  we should not be tempted to continue the parsing.

                  We log a error message with the name of the problematic
                  keyword and the location in the input deck. We rethrow the
                  same exception without updating the what() message of the
                  exception.
                */
                const OpmInputError opm_error { e, rawKeyword->location() } ;

                OpmLog::error(opm_error.what());

                std::throw_with_nested(opm_error);
            }
        } else {
            const std::string msg = "The keyword " + rawKeyword->getKeywordName() + " is not recognized - ignored";
            KeywordLocation location(rawKeyword->getKeywordName(), parserState.current_path().string(), parserState.line());
            OpmLog::warning(Log::fileMessage(location, msg));
        }
    }

    return true;
}

}


    /* stripComments only exists so that the unit tests can verify it.
     * strip_comment is the actual (internal) implementation
     */
    std::string Parser::stripComments( const std::string& str ) {
        return { str.begin(),
                 str::find_terminator( str.begin(), str.end(), str::find_comment() ) };
    }

    Parser::Parser(bool addDefault) {
        // The addDefaultKeywords() method is implemented in a source file
        // ${PROJECT_BINARY_DIR}/ParserInit.cpp which is generated by the build
        // system.

        if (addDefault)
            this->addDefaultKeywords();
    }


    /*
     About INCLUDE: Observe that the ECLIPSE parser is slightly unlogical
     when it comes to nested includes; the path to an included file is always
     interpreted relative to the filesystem location of the DATA file, and
     not the location of the file issuing the INCLUDE command. That behaviour
     is retained in the current implementation.
     */

    inline void assertFullDeck(const ParseContext& context) {
        if (context.hasKey(ParseContext::PARSE_MISSING_SECTIONS))
            throw new std::logic_error("Cannot construct a state in partial deck context");
    }

    EclipseState Parser::parse(const std::string &filename, const ParseContext& context, ErrorGuard& errors) {
        assertFullDeck(context);
        return EclipseState( Parser{}.parseFile( filename, context, errors ));
    }

    EclipseState Parser::parse(const Deck& deck, const ParseContext& context) {
        assertFullDeck(context);
        return EclipseState(deck);
    }

    EclipseState Parser::parseData(const std::string &data, const ParseContext& context, ErrorGuard& errors) {
        assertFullDeck(context);
        Parser p;
        auto deck = p.parseString(data, context, errors);
        return parse(deck, context);
    }

    EclipseGrid Parser::parseGrid(const std::string &filename, const ParseContext& context , ErrorGuard& errors) {
        if (context.hasKey(ParseContext::PARSE_MISSING_SECTIONS))
            return EclipseGrid{ filename };
        return parse(filename, context, errors).getInputGrid();
    }

    EclipseGrid Parser::parseGrid(const Deck& deck, const ParseContext& context)
    {
        if (context.hasKey(ParseContext::PARSE_MISSING_SECTIONS))
            return EclipseGrid{ deck };
        return parse(deck, context).getInputGrid();
    }

    EclipseGrid Parser::parseGridData(const std::string &data, const ParseContext& context, ErrorGuard& errors) {
        Parser parser;
        auto deck = parser.parseString(data, context, errors);
        if (context.hasKey(ParseContext::PARSE_MISSING_SECTIONS)) {
            return EclipseGrid{ deck };
        }
        return parse(deck, context).getInputGrid();
    }

    Deck Parser::parseFile(const std::string &dataFileName, const ParseContext& parseContext,
                           ErrorGuard& errors, const std::vector<Opm::Ecl::SectionType>& sections) const {


        std::set<Opm::Ecl::SectionType> ignore_sections;

        if (sections.size() > 0) {

            std::set<Opm::Ecl::SectionType> all_sections;
            all_sections = {Opm::Ecl::RUNSPEC, Opm::Ecl::GRID, Opm::Ecl::EDIT, Opm::Ecl::PROPS, Opm::Ecl::REGIONS,
                            Opm::Ecl::SOLUTION, Opm::Ecl::SUMMARY, Opm::Ecl::SCHEDULE};

            std::set<Opm::Ecl::SectionType> read_sections;

            for (auto sec : sections)
                 read_sections.insert(sec);

            std::set_difference(all_sections.begin(), all_sections.end(), read_sections.begin(), read_sections.end(),
                            std::inserter(ignore_sections, ignore_sections.end()));
        }

        /*
          The following rules apply to the .DATA file argument which is
          internalized in the deck:

           1. It is normalized by removing uneccessary '.' characters and
              resolving symlinks.

           2. The relative/abolute status of the path is retained.
        */

        std::string data_file;
        if (dataFileName[0] == '/')
            data_file = std::filesystem::canonical(dataFileName).generic_string();
        else
            data_file = std::filesystem::proximate(std::filesystem::canonical(dataFileName)).generic_string();

        ParserState parserState( this->codeKeywords(), parseContext, errors, data_file, ignore_sections);
        parseState( parserState, *this );

        auto ignore = parserState.get_ignore();

        if (ignore.size() > 0)
            cleanup_deck_keyword_list(parserState, ignore);

        return std::move( parserState.deck );
    }

    Deck Parser::parseFile(const std::string& dataFileName,
                           const ParseContext& parseContext) const {
        ErrorGuard errors;
        return this->parseFile(dataFileName, parseContext, errors, {});
    }

    Deck Parser::parseFile(const std::string& dataFileName,
                           const ParseContext& parseContext,
                           const std::vector<Opm::Ecl::SectionType>& sections) const {
        ErrorGuard errors;
        return this->parseFile(dataFileName, parseContext, errors, sections);
    }

    Deck Parser::parseFile(const std::string& dataFileName) const {
        ErrorGuard errors;
        return this->parseFile(dataFileName, ParseContext(), errors);
    }




    Deck Parser::parseString(const std::string &data, const ParseContext& parseContext, ErrorGuard& errors) const {
        ParserState parserState( this->codeKeywords(), parseContext, errors );
        parserState.loadString( data );
        parseState( parserState, *this );
        return std::move( parserState.deck );
    }

    Deck Parser::parseString(const std::string &data, const ParseContext& parseContext) const {
        ErrorGuard errors;
        return this->parseString(data, parseContext, errors);
    }

    Deck Parser::parseString(const std::string &data) const {
        ErrorGuard errors;
        return this->parseString(data, ParseContext(), errors);
    }

    size_t Parser::size() const {
        return m_deckParserKeywords.size();
    }

    const ParserKeyword* Parser::matchingKeyword(const std::string_view& name) const
    {
        const auto it = std::find_if(m_wildCardKeywords.begin(),
                                     m_wildCardKeywords.end(),
                                     [&name](const auto& wild)
                                     {
                                         return wild.second->matches(name);
                                     });
        return it != m_wildCardKeywords.end() ? it->second : nullptr;
    }

    bool Parser::hasWildCardKeyword(const std::string& internalKeywordName) const {
        return (m_wildCardKeywords.count(internalKeywordName) > 0);
    }

    bool Parser::isRecognizedKeyword(std::string_view name) const
    {
        if (! ParserKeyword::validDeckName(name)) {
            return false;
        }

        return (this->m_deckParserKeywords.find(name) != this->m_deckParserKeywords.end())
            || (this->matchingKeyword(name) != nullptr);
    }

    bool Parser::isBaseRecognizedKeyword(std::string_view name) const
    {
        return ParserKeyword::validDeckName(name)
            && (this->m_deckParserKeywords.find(name) != this->m_deckParserKeywords.end());
    }

void Parser::addParserKeyword( ParserKeyword parserKeyword ) {
    /* Store the keywords in the keyword storage. They aren't free'd until the
     * parser gets destroyed, even if there is no reasonable way to reach them
     * (effectively making them leak). This is not a big problem because:
     *
     * * A keyword can be added that overwrites some *but not all* deckname ->
     *   keyword mappings. Keeping track of this is more hassle than worth for
     *   what is essentially edge case usage.
     * * We can store (and search) via std::string_view's from the keyword added
     *   first because we know that it will be kept around, i.e. we don't have to
     *   deal with subtle lifetime issues.
     * * It means we aren't reliant on some internal name mapping, and can only
     * be concerned with interesting behaviour.
     * * Finally, these releases would in practice never happen anyway until
     *   the parser went out of scope, and now they'll also be cleaned up in the
     *   same sweep.
     */

    this->keyword_storage.push_back( std::move( parserKeyword ) );
    const ParserKeyword * ptr = std::addressof(this->keyword_storage.back());
    for (const auto& deck_name : ptr->deck_names())
    {
        m_deckParserKeywords[deck_name] = ptr;
    }

    if (ptr->hasMatchRegex()) {
        std::string_view name( ptr->getName() );
        m_wildCardKeywords[ name ] = ptr;
    }

    if (ptr->isCodeKeyword())
        this->code_keywords.emplace_back( ptr->getName(), ptr->codeEnd() );
}


void Parser::addParserKeyword(const Json::JsonObject& jsonKeyword) {
    addParserKeyword( ParserKeyword( jsonKeyword ) );
}

bool Parser::hasKeyword( const std::string& name ) const {
    return this->m_deckParserKeywords.find( std::string_view( name ) )
        != this->m_deckParserKeywords.end();
}

const ParserKeyword& Parser::getKeyword( const std::string& name ) const {
    return getParserKeywordFromDeckName( std::string_view( name ) );
}

const ParserKeyword& Parser::getParserKeywordFromDeckName(const std::string_view& name ) const {
    auto candidate = m_deckParserKeywords.find( name );

    if( candidate != m_deckParserKeywords.end() ) return *candidate->second;

    const auto* wildCardKeyword = matchingKeyword( name );

    if ( !wildCardKeyword )
        throw std::invalid_argument( "Do not have parser keyword for parsing: " + std::string(name) );

    return *wildCardKeyword;
}

std::vector<std::string> Parser::getAllDeckNames () const {
    std::vector<std::string> keywords;
    for (auto iterator = m_deckParserKeywords.begin(); iterator != m_deckParserKeywords.end(); iterator++) {
        keywords.push_back(std::string(iterator->first));
    }
    for (auto iterator = m_wildCardKeywords.begin(); iterator != m_wildCardKeywords.end(); iterator++) {
        keywords.push_back(std::string(iterator->first));
    }
    return keywords;
}


    void Parser::loadKeywords(const Json::JsonObject& jsonKeywords) {
        if (jsonKeywords.is_array()) {
            for (size_t index = 0; index < jsonKeywords.size(); index++) {
                Json::JsonObject jsonKeyword = jsonKeywords.get_array_item(index);
                addParserKeyword( ParserKeyword( jsonKeyword ) );
            }
        } else
            throw std::invalid_argument("Input JSON object is not an array");
    }

    bool Parser::loadKeywordFromFile(const std::filesystem::path& configFile) {

        try {
            Json::JsonObject jsonKeyword(configFile);
            addParserKeyword( ParserKeyword( jsonKeyword ) );
            return true;
        }
        catch (...) {
            return false;
        }
    }


    void Parser::loadKeywordsFromDirectory(const std::filesystem::path& directory, bool recursive) {
        if (!std::filesystem::exists(directory))
            throw std::invalid_argument("Directory: " + directory.string() + " does not exist.");
        else {
            std::filesystem::directory_iterator end;
            for (std::filesystem::directory_iterator iter(directory); iter != end; iter++) {
                if (std::filesystem::is_directory(*iter)) {
                    if (recursive)
                        loadKeywordsFromDirectory(*iter, recursive);
                } else {
                    if (ParserKeyword::validInternalName(iter->path().filename().string())) {
                        if (!loadKeywordFromFile(*iter))
                            std::cerr << "** Warning: failed to load keyword from file:" << iter->path() << std::endl;
                    }
                }
            }
        }
    }

    const std::vector<std::pair<std::string,std::string>> Parser::codeKeywords() const {
        return this->code_keywords;
    }


#if 0
    void Parser::applyUnitsToDeck(Deck& deck) const {

        for( auto& deckKeyword : deck ) {

            if( !isRecognizedKeyword( deckKeyword.name() ) ) continue;

            const auto& parserKeyword = getParserKeywordFromDeckName( deckKeyword.name() );
            if( !parserKeyword.hasDimension() ) continue;

            parserKeyword.applyUnitsToDeck(deck , deckKeyword);
        }
    }
#endif


    static bool isSectionDelimiter( const DeckKeyword& keyword )
    {
        const auto& name = keyword.name();
        using namespace std::string_view_literals;
        const auto delimiters = std::array {
            "RUNSPEC"sv,
            "GRID"sv,
            "EDIT"sv,
            "PROPS"sv,
            "REGIONS"sv,
            "SOLUTION"sv,
            "SUMMARY"sv,
            "SCHEDULE"sv
        };

        return std::find(delimiters.begin(), delimiters.end(), name) != delimiters.end();
    }

    bool DeckSection::checkSectionTopology(const Deck& deck,
                                           const Parser& parser,
                                           Opm::ErrorGuard& errorGuard,
                                           bool ensureKeywordSectionAffiliation)
    {
        if( deck.size() == 0 ) {
            std::string msg = "empty decks are invalid\n";
            OpmLog::warning(msg);
            return false;
        }

        bool deckValid = true;
        const std::string errorKey = "SECTION_TOPOLOGY_ERROR";
        // We put errors on the top level to the end of the list to make them
        // more prominent
        std::vector<std::string> topLevelErrors;
        size_t curKwIdx = 0;

        while(curKwIdx < deck.size() && isGlobalKeyword(deck[curKwIdx])) {
            ++curKwIdx;
        }

        bool validKwIdx = curKwIdx < deck.size();
        if( !validKwIdx || deck[curKwIdx].name() != "RUNSPEC" ) {
            constexpr std::string_view msg = "The first keyword of a valid deck must be RUNSPEC (is {})\n";
            auto curKeyword = deck[0];
            topLevelErrors.push_back(Log::fileMessage(curKeyword.location(),
                                                      fmt::format(msg,
                                                                  (validKwIdx)? deck[curKwIdx].name() : "")));
            deckValid = false;
        }

        std::string curSectionName = validKwIdx? deck[curKwIdx].name() : "";

        for (++curKwIdx; curKwIdx < deck.size(); ++curKwIdx) {
                auto checker = [&errorGuard, &deckValid, &parser, curSectionName,
                                ensureKeywordSectionAffiliation, errorKey]
                    (const std::string& curKeywordName, const KeywordLocation& location)
                {
                    const auto& parserKeyword =
                        parser.getParserKeywordFromDeckName( curKeywordName );
                    if (ensureKeywordSectionAffiliation && !parserKeyword.isValidSection(curSectionName)) {
                        constexpr std::string_view msg =
                            "The keyword '{}' is located in the '{}' section where it is invalid";
                        errorGuard.addError(errorKey,
                                            Log::fileMessage(location,
                                                             fmt::format(msg,
                                                                         curKeywordName,
                                                                         curSectionName)
                                                         ) );
                        deckValid = false;
                    }
                };

            const auto& curKeyword = deck[curKwIdx];
            const std::string& curKeywordName = curKeyword.name();

            if (!isSectionDelimiter( curKeyword )) {
                if( !parser.isRecognizedKeyword( curKeywordName ) )
                    // ignore unknown keywords for now (i.e. they can appear in any section)
                    continue;

                const bool isOperateKeyword =
                    Fieldprops::keywords::is_oper_keyword(curKeywordName);

                if (isOperateKeyword) {
                    for (const auto& record : curKeyword) {
                        const auto& operName = record.getItem(0).getTrimmedString(0);
                        if (!parser.isRecognizedKeyword(operName)) {
                            // ignore unknown keywords
                            continue;
                        }
                        checker(operName, curKeyword.location());
                    }
                }
                else {
                    checker(curKeyword.name(), curKeyword.location());
                }

                continue;
            }

            if (curSectionName == "RUNSPEC") {
                if (curKeywordName != "GRID") {
                    std::string msg =
                        "The RUNSPEC section must be followed by GRID instead of "+curKeywordName;
                    topLevelErrors.push_back(Log::fileMessage(curKeyword.location(), msg) );
                    deckValid = false;
                }

                curSectionName = curKeywordName;
            }
            else if (curSectionName == "GRID") {
                if (curKeywordName != "EDIT" && curKeywordName != "PROPS") {
                    std::string msg =
                        "The GRID section must be followed by EDIT or PROPS instead of "+curKeywordName;
                    topLevelErrors.push_back(Log::fileMessage(curKeyword.location(), msg) );
                    deckValid = false;
                }

                curSectionName = curKeywordName;
            }
            else if (curSectionName == "EDIT") {
                if (curKeywordName != "PROPS") {
                    std::string msg =
                        "The EDIT section must be followed by PROPS instead of "+curKeywordName;
                    topLevelErrors.push_back(Log::fileMessage(curKeyword.location(), msg) );
                    deckValid = false;
                }

                curSectionName = curKeywordName;
            }
            else if (curSectionName == "PROPS") {
                if (curKeywordName != "REGIONS" && curKeywordName != "SOLUTION") {
                    std::string msg =
                        "The PROPS section must be followed by REGIONS or SOLUTION instead of "+curKeywordName;
                    topLevelErrors.push_back(Log::fileMessage(curKeyword.location(), msg) );
                    deckValid = false;
                }

                curSectionName = curKeywordName;
            }
            else if (curSectionName == "REGIONS") {
                if (curKeywordName != "SOLUTION") {
                    std::string msg =
                        "The REGIONS section must be followed by SOLUTION instead of "+curKeywordName;
                    topLevelErrors.push_back(Log::fileMessage(curKeyword.location(), msg) );
                    deckValid = false;
                }

                curSectionName = curKeywordName;
            }
            else if (curSectionName == "SOLUTION") {
                if (curKeywordName != "SUMMARY" && curKeywordName != "SCHEDULE") {
                    std::string msg =
                        "The SOLUTION section must be followed by SUMMARY or SCHEDULE instead of "+curKeywordName;
                    topLevelErrors.push_back(Log::fileMessage(curKeyword.location(), msg) );
                    deckValid = false;
                }

                curSectionName = curKeywordName;
            }
            else if (curSectionName == "SUMMARY") {
                if (curKeywordName != "SCHEDULE") {
                    std::string msg =
                        "The SUMMARY section must be followed by SCHEDULE instead of "+curKeywordName;
                    topLevelErrors.push_back(Log::fileMessage(curKeyword.location(), msg) );
                    deckValid = false;
                }

                curSectionName = curKeywordName;
            }
            else if (curSectionName == "SCHEDULE") {
                // schedule is the last section, so every section delimiter after it is wrong...
                std::string msg =
                    "The SCHEDULE section must be the last one ("
                    +curKeywordName+" specified after SCHEDULE)";
                topLevelErrors.push_back(Log::fileMessage(curKeyword.location(), msg) );
                deckValid = false;
            }
        }

        // SCHEDULE is the last section and it is mandatory, so make sure it is there
        if (curSectionName != "SCHEDULE") {
            const auto& curKeyword = deck[deck.size() - 1];
            std::string msg =
                "The last section of a valid deck must be SCHEDULE (is "+curSectionName+")";
            topLevelErrors.push_back(Log::fileMessage(curKeyword.location(), msg) );
            deckValid = false;
        }

        for(const auto& err: topLevelErrors) {
            errorGuard.addError(errorKey, err);
        }

        return deckValid;
    }

} // namespace Opm
