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

#ifndef FIELDPROPS_HPP
#define FIELDPROPS_HPP

#include <opm/common/utility/OpmInputError.hpp>

#include <opm/input/eclipse/EclipseState/Grid/Box.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FieldData.hpp>
#include <opm/input/eclipse/EclipseState/Grid/Keywords.hpp>
#include <opm/input/eclipse/EclipseState/Grid/SatfuncPropertyInitializers.hpp>
#include <opm/input/eclipse/EclipseState/Grid/TranCalculator.hpp>
#include <opm/input/eclipse/EclipseState/Runspec.hpp>
#include <opm/input/eclipse/EclipseState/Util/OrderedMap.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>

#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <opm/input/eclipse/Deck/DeckSection.hpp>
#include <opm/input/eclipse/Deck/value_status.hpp>

#include <cstddef>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace Opm {

class Deck;
class EclipseGrid;
class NumericalAquifers;

namespace Fieldprops
{

namespace keywords {

/*
  Regarding global keywords
  =========================

  It turns out that when the option 'ALL' is used for the PINCH keyword we
  require the MULTZ keyword specified for all cells, also the inactive cells.
  The premise for the FieldProps implementation has all the way been that only
  the active cells should be stored.

  In order to support the ALL option of the PINCH keyword we have bolted on a
  limited support for global storage. By setting .global = true in the
  keyword_info describing the keyword you get:

  1. Normal deck assignment like

       MULTZ
           ..... /

  2. Scalar operations like EQUALS and MULTIPLY.

  These operations also support the full details of the BOX behavior.

  The following operations do not work
  ------------------------------------

  1. Operations involving multiple keywords like

       COPY
         MULTX  MULTZ /
       /

    this also includes the OPERATE which involves multiple keywords for some
    of its operations.

  2. All region operatins like EQUALREG and MULTREG.

  The operations which are not properly implemented will be intercepted and a
  std::logic_error() exception will be thrown.
*/



inline bool isFipxxx(const std::string& keyword) {
    // FIPxxxx can be any keyword, e.g. FIPREG or FIPXYZ that has the pattern "FIP.+"
    // However, it can not be FIPOWG as that is an actual keyword.
    if (keyword.size() < 4 || keyword == "FIPOWG") {
        return false;
    }
    return keyword[0] == 'F' && keyword[1] == 'I' && keyword[2] == 'P';
}


/*
  The aliased_keywords map defines aliases for other keywords. The FieldProps
  objects will translate those keywords before further processing. The aliases
  will also be exposed by the FieldPropsManager object.

  However, the following methods of FieldProps do not fully support aliases:
  - FieldProps::keys() does not return the aliases.
  - FieldProps::erase() and FieldProps::extract() do not support aliases. Using
    them with an aliased keyword will also remove the alias.

  Note that the aliases are also added to GRID::double_keywords.

  The PERMR and PERMTHT keywords are aliases for PERMX and PERMY, respectively.
*/
namespace ALIAS {
    static const std::unordered_map<std::string, std::string> aliased_keywords = {{"PERMR", "PERMX"},
                                                                                  {"PERMTHT", "PERMY"}};
}


namespace GRID {
static const std::unordered_map<std::string, keyword_info<double>> double_keywords = {{"DISPERC",keyword_info<double>{}.unit_string("Length")},
                                                                                      {"MINPVV",  keyword_info<double>{}.init(0.0).unit_string("ReservoirVolume").global_kw(true)},
                                                                                      {"MULTPV",  keyword_info<double>{}.init(1.0).mult(true)},
                                                                                      {"NTG",     keyword_info<double>{}.init(1.0)},
                                                                                      {"PORO",    keyword_info<double>{}.distribute_top(true)},
                                                                                      {"PERMX",   keyword_info<double>{}.unit_string("Permeability").distribute_top(true).global_kw_until_edit()},
                                                                                      {"PERMY",   keyword_info<double>{}.unit_string("Permeability").distribute_top(true).global_kw_until_edit()},
                                                                                      {"PERMZ",   keyword_info<double>{}.unit_string("Permeability").distribute_top(true).global_kw_until_edit()},
                                                                                      {"PERMR",   keyword_info<double>{}.unit_string("Permeability").distribute_top(true).global_kw_until_edit()},
                                                                                      {"PERMTHT",   keyword_info<double>{}.unit_string("Permeability").distribute_top(true).global_kw_until_edit()},
                                                                                      {"TEMPI",   keyword_info<double>{}.unit_string("Temperature")},
                                                                                      {"THCONR",  keyword_info<double>{}.unit_string("Energy/AbsoluteTemperature*Length*Time")},
                                                                                      {"THCONSF", keyword_info<double>{}},
                                                                                      {"HEATCR",  keyword_info<double>{}.unit_string("Energy/ReservoirVolume*AbsoluteTemperature")},
                                                                                      {"HEATCRT", keyword_info<double>{}.unit_string("Energy/ReservoirVolume*AbsoluteTemperature*AbsoluteTemperature")},
                                                                                      {"THCROCK", keyword_info<double>{}.unit_string("Energy/AbsoluteTemperature*Length*Time")},
                                                                                      {"THCOIL",  keyword_info<double>{}.unit_string("Energy/AbsoluteTemperature*Length*Time")},
                                                                                      {"THCGAS",  keyword_info<double>{}.unit_string("Energy/AbsoluteTemperature*Length*Time")},
                                                                                      {"THCWATER",keyword_info<double>{}.unit_string("Energy/AbsoluteTemperature*Length*Time")},
                                                                                      {"YMODULE", keyword_info<double>{}.unit_string("Giga*Pascal")},
                                                                                      {"PRATIO", keyword_info<double>{}.unit_string("1")},
                                                                                      {"BIOTCOEF", keyword_info<double>{}.unit_string("1")},
                                                                                      {"POELCOEF", keyword_info<double>{}.unit_string("1")},
                                                                                      {"THERMEXR", keyword_info<double>{}.unit_string("1/AbsoluteTemperature")},
                                                                                      {"THELCOEF", keyword_info<double>{}.unit_string("Pressure/AbsoluteTemperature")},
                                                                                      {"MULTX",   keyword_info<double>{}.init(1.0).mult(true)},
                                                                                      {"MULTX-",  keyword_info<double>{}.init(1.0).mult(true)},
                                                                                      {"MULTY",   keyword_info<double>{}.init(1.0).mult(true)},
                                                                                      {"MULTY-",  keyword_info<double>{}.init(1.0).mult(true)},
                                                                                      {"MULTZ",   keyword_info<double>{}.init(1.0).mult(true).global_kw(true)},
                                                                                      {"MULTZ-",  keyword_info<double>{}.init(1.0).mult(true).global_kw(true)}};

static const std::unordered_map<std::string, keyword_info<int>> int_keywords = {{"ACTNUM",  keyword_info<int>{}.init(1)},
                                                                                {"FLUXNUM", keyword_info<int>{}},
                                                                                {"ISOLNUM", keyword_info<int>{}.init(1)},
                                                                                {"MULTNUM", keyword_info<int>{}.init(1)},
                                                                                {"OPERNUM", keyword_info<int>{}},
                                                                                {"ROCKNUM", keyword_info<int>{}}};

}

namespace EDIT {

/*
  The TRANX, TRANY and TRANZ properties are handled very differently from the
  other properties. It is important that these fields are not entered into the
  double_keywords list of the EDIT section, that way we risk silent failures
  due to the special treatment of the TRAN fields.
*/

static const std::unordered_map<std::string, keyword_info<double>> double_keywords = {{"MULTPV",  keyword_info<double>{}.init(1.0).mult(true)},
                                                                                      {"PORV",    keyword_info<double>{}.unit_string("ReservoirVolume")},
                                                                                      {"MULTX",   keyword_info<double>{}.init(1.0).mult(true)},
                                                                                      {"MULTX-",  keyword_info<double>{}.init(1.0).mult(true)},
                                                                                      {"MULTY",   keyword_info<double>{}.init(1.0).mult(true)},
                                                                                      {"MULTY-",  keyword_info<double>{}.init(1.0).mult(true)},
                                                                                      {"MULTZ",   keyword_info<double>{}.init(1.0).mult(true).global_kw(true)},
                                                                                      {"MULTZ-",  keyword_info<double>{}.init(1.0).mult(true).global_kw(true)}};

static const std::unordered_map<std::string, keyword_info<int>> int_keywords = {};
}

namespace PROPS {
static const std::unordered_map<std::string, keyword_info<double>> double_keywords = {{"SWATINIT", keyword_info<double>{}},
                                                                                      {"PCG",      keyword_info<double>{}.unit_string("Pressure")},
                                                                                      {"IPCG",     keyword_info<double>{}.unit_string("Pressure")},
                                                                                      {"PCW",      keyword_info<double>{}.unit_string("Pressure")},
                                                                                      {"IPCW",     keyword_info<double>{}.unit_string("Pressure")}};
static const std::unordered_map<std::string, keyword_info<int>> int_keywords = {};

#define dirfunc(base) base, base "X", base "X-", base "Y", base "Y-", base "Z", base "Z-"

static const std::set<std::string> satfunc = {"SWLPC", "ISWLPC", "SGLPC", "ISGLPC",
                                              dirfunc("SGL"),
                                              dirfunc("ISGL"),
                                              dirfunc("SGU"),
                                              dirfunc("ISGU"),
                                              dirfunc("SWL"),
                                              dirfunc("ISWL"),
                                              dirfunc("SWU"),
                                              dirfunc("ISWU"),
                                              dirfunc("SGCR"),
                                              dirfunc("ISGCR"),
                                              dirfunc("SOWCR"),
                                              dirfunc("ISOWCR"),
                                              dirfunc("SOGCR"),
                                              dirfunc("ISOGCR"),
                                              dirfunc("SWCR"),
                                              dirfunc("ISWCR"),
                                              dirfunc("KRW"),
                                              dirfunc("IKRW"),
                                              dirfunc("KRWR"),
                                              dirfunc("IKRWR"),
                                              dirfunc("KRO"),
                                              dirfunc("IKRO"),
                                              dirfunc("KRORW"),
                                              dirfunc("IKRORW"),
                                              dirfunc("KRORG"),
                                              dirfunc("IKRORG"),
                                              dirfunc("KRG"),
                                              dirfunc("IKRG"),
                                              dirfunc("KRGR"),
                                              dirfunc("IKRGR")};

#undef dirfunc
}

namespace REGIONS {

static const std::unordered_map<std::string, keyword_info<int>> int_keywords = {{"ENDNUM",   keyword_info<int>{}.init(1)},
                                                                                {"EOSNUM",   keyword_info<int>{}.init(1)},
                                                                                {"EQLNUM",   keyword_info<int>{}.init(1)},
                                                                                {"FIPNUM",   keyword_info<int>{}.init(1)},
                                                                                {"IMBNUM",   keyword_info<int>{}.init(1)},
                                                                                {"OPERNUM",  keyword_info<int>{}},
                                                                                {"STRESSEQUILNUM", keyword_info<int>{}.init(1)},
                                                                                {"MISCNUM",  keyword_info<int>{}},
                                                                                {"MISCNUM",  keyword_info<int>{}},
                                                                                {"PVTNUM",   keyword_info<int>{}.init(1)},
                                                                                {"SATNUM",   keyword_info<int>{}.init(1)},
                                                                                {"LWSLTNUM", keyword_info<int>{}},
                                                                                {"ROCKNUM",  keyword_info<int>{}},
                                                                                {"KRNUMX",   keyword_info<int>{}},
                                                                                {"KRNUMY",   keyword_info<int>{}},
                                                                                {"KRNUMZ",   keyword_info<int>{}},
                                                                                {"IMBNUMX",   keyword_info<int>{}},
                                                                                {"IMBNUMY",   keyword_info<int>{}},
                                                                                {"IMBNUMZ",   keyword_info<int>{}},
                                                                                };
}

namespace SOLUTION {

static const std::unordered_map<std::string, keyword_info<double>> double_keywords = {{"PRESSURE", keyword_info<double>{}.unit_string("Pressure")},
                                                                                      {"SPOLY",    keyword_info<double>{}.unit_string("Density")},
                                                                                      {"SPOLYMW",  keyword_info<double>{}},
                                                                                      {"SSOL",     keyword_info<double>{}},
                                                                                      {"SWAT",     keyword_info<double>{}},
                                                                                      {"SGAS",     keyword_info<double>{}},
                                                                                      {"SMICR",    keyword_info<double>{}.unit_string("Density")},
                                                                                      {"SOXYG",    keyword_info<double>{}.unit_string("Density")},
                                                                                      {"SUREA",    keyword_info<double>{}.unit_string("Density")},
                                                                                      {"SBIOF",    keyword_info<double>{}},
                                                                                      {"SCALC",    keyword_info<double>{}},
                                                                                      {"SALTP",    keyword_info<double>{}},
                                                                                      {"SALT",     keyword_info<double>{}.unit_string("Salinity")},
                                                                                      {"TEMPI",    keyword_info<double>{}.unit_string("Temperature")},
                                                                                      {"RS",       keyword_info<double>{}.unit_string("GasDissolutionFactor")},
                                                                                      {"RV",       keyword_info<double>{}.unit_string("OilDissolutionFactor")},
                                                                                      {"RVW",      keyword_info<double>{}.unit_string("OilDissolutionFactor")}
                                                                                      };

static const std::unordered_map<std::string, keyword_info<double>> composition_keywords = {{"XMF", keyword_info<double>{}},
                                                                                           {"YMF", keyword_info<double>{}},
                                                                                           {"ZMF", keyword_info<double>{}},
                                                                                          };
}

namespace SCHEDULE {

static const std::unordered_map<std::string, keyword_info<double>> double_keywords = {{"MULTX",   keyword_info<double>{}.init(1.0).mult(true)},
                                                                                      {"MULTX-",  keyword_info<double>{}.init(1.0).mult(true)},
                                                                                      {"MULTY",   keyword_info<double>{}.init(1.0).mult(true)},
                                                                                      {"MULTY-",  keyword_info<double>{}.init(1.0).mult(true)},
                                                                                      {"MULTZ",   keyword_info<double>{}.init(1.0).mult(true).global_kw(true)},
                                                                                      {"MULTZ-",  keyword_info<double>{}.init(1.0).mult(true).global_kw(true)}};

static const std::unordered_map<std::string, keyword_info<int>> int_keywords = {{"ROCKNUM",   keyword_info<int>{}}};

}

template <typename T>
keyword_info<T> global_kw_info(const std::string& name, bool allow_unsupported = false);

bool is_oper_keyword(const std::string& name);
} // end namespace keywords

} // end namespace Fieldprops

class FieldProps {
public:

    using ScalarOperation = Fieldprops::ScalarOperation;

    struct MultregpRecord {
        int region_value;
        double multiplier;
        std::string region_name;


        MultregpRecord(int rv, double m, const std::string& rn) :
            region_value(rv),
            multiplier(m),
            region_name(rn)
        {}


        bool operator==(const MultregpRecord& other) const {
            return this->region_value == other.region_value &&
                   this->multiplier == other.multiplier &&
                   this->region_name == other.region_name;
        }
    };

    /// Property array existence status
    enum class GetStatus {
        /// Property exists and its property data is fully defined.
        OK = 1,

        /// Property array has not been fully initialised.
        ///
        /// At least one data element is uninitialised or has an empty
        /// default.
        ///
        /// Will throw an exception of type \code std::runtime_error
        /// \endcode in FieldDataManager::verify_status().
        INVALID_DATA = 2,

        /// Property has not yet been defined in the input file.
        ///
        /// Will throw an exception of type \code std::out_of_range \endcode
        /// in FieldDataManager::verify_status().
        MISSING_KEYWORD = 3,

        /// Named property is not known to the internal handling mechanism.
        ///
        /// Might for instance be because the client code requested a
        /// property of type \c int, but the name is only known as property
        /// of type \c double, or because there was a misprint in the
        /// property name.
        ///
        /// Will throw an exception of type \code std::logic_error \endcode
        /// in FieldDataManager::verify_status().
        NOT_SUPPPORTED_KEYWORD = 4,
    };

    /// Wrapper type for field properties
    ///
    /// \tparam T Property element type.  Typically \c double or \c int.
    template<typename T>
    struct FieldDataManager
    {
        /// Property name.
        const std::string& keyword;

        /// Request status.
        GetStatus status;

        /// Property data.
        const Fieldprops::FieldData<T>* data_ptr;

        /// Constructor
        ///
        /// \param[in] k Property name
        /// \param[in] s Request status
        /// \param[in] d Property data.  Pass \c nullptr for missing property data.
        FieldDataManager(const std::string& k, GetStatus s, const Fieldprops::FieldData<T>* d)
            : keyword(k)
            , status(s)
            , data_ptr(d)
        {}

        /// Validate result of \code try_get<>() \endcode request
        ///
        /// Throws an exception of type \code OpmInputError \endcode if
        /// \code this->status \endcode is not \code GetStatus::OK \endcode.
        /// Does nothing otherwise.
        ///
        /// \param[in] Input keyword which prompted request.
        ///
        /// \param[in] descr Textual description of context in which request
        ///   occurred.
        ///
        /// \param[in] operation Name of operation which prompted request.
        void verify_status(const KeywordLocation& loc,
                           const std::string&     descr,
                           const std::string&     operation) const
        {
            switch (this->status) {
            case FieldProps::GetStatus::OK:
                return;

            case FieldProps::GetStatus::INVALID_DATA:
                throw OpmInputError {
                    descr + " " + this->keyword +
                    " is not fully initialised for " + operation,
                    loc
                };

            case FieldProps::GetStatus::MISSING_KEYWORD:
                throw OpmInputError {
                    descr + " " + this->keyword +
                    " does not exist in input deck for " + operation,
                    loc
                };

            case FieldProps::GetStatus::NOT_SUPPPORTED_KEYWORD:
                throw OpmInputError {
                    descr + " " + this->keyword +
                    " is not supported for " + operation,
                    loc
                };
            }
        }

        /// Validate result of \code try_get<>() \endcode request
        ///
        /// Throws an exception as outlined in \c GetStatus if \code
        /// this->status \endcode is not \code GetStatus::OK \endcode.  Does
        /// nothing otherwise.
        void verify_status() const
        {
            switch (status) {
            case FieldProps::GetStatus::OK:
                return;

            case FieldProps::GetStatus::INVALID_DATA:
                throw std::runtime_error("The keyword: " + keyword + " has not been fully initialized");

            case FieldProps::GetStatus::MISSING_KEYWORD:
                throw std::out_of_range("No such keyword in deck: " + keyword);

            case FieldProps::GetStatus::NOT_SUPPPORTED_KEYWORD:
                throw std::logic_error("The keyword  " + keyword + " is not supported");
            }
        }

        /// Access underlying property data elements
        ///
        /// Returns \c nullptr if property data is not fully defined.
        const std::vector<T>* ptr() const
        {
            return (this->data_ptr != nullptr)
                ? &this->data_ptr->data
                : nullptr;
        }

        /// Access underlying property data elements
        ///
        /// Throws an exception as outlined in \c GetStatus if property data
        /// is not fully defined.
        const std::vector<T>& data() const
        {
            this->verify_status();
            return this->data_ptr->data;
        }

        /// Read-only access to contained FieldData object
        ///
        /// Throws an exception as outlined in \c GetStatus if property data
        /// is not fully defined.
        const Fieldprops::FieldData<T>& field_data() const
        {
            this->verify_status();
            return *this->data_ptr;
        }

        /// Property validity predicate.
        ///
        /// Returns true if property exists and has fully defined data
        /// elements.  False otherwise.
        bool valid() const
        {
            return this->status == GetStatus::OK;
        }
    };

    /// Options to restrict or relax a try_get() request.
    ///
    /// For instance, the source array of a COPY operation must exist and be
    /// fully defined, whereas the target array need not be either.  On the
    /// other hand, certain use cases might want to look up property names
    /// of an unsupported keyword type--e.g., transmissibility keywords
    /// among the integer properties--and those requests should specify the
    /// AllowUnsupported flag.
    ///
    /// Note: We would ideally use strong enums (i.e., "enum class") for
    /// this, but those don't mesh very well with bitwise operations.
    enum TryGetFlags : unsigned int {
        /// Whether or not to permit looking up property names of unmatching
        /// types.
        AllowUnsupported = (1u << 0),

        /// Whether or not the property must already exist.
        MustExist = (1u << 1),
    };

    /// Normal constructor for FieldProps.
    FieldProps(const Deck& deck,
               const Phases& phases,
               EclipseGrid& grid,
               const TableManager& table_arg,
               const std::size_t ncomps);

    /// Special case constructor used to process ACTNUM only.
    FieldProps(const Deck& deck, const EclipseGrid& grid);

    void reset_actnum(const std::vector<int>& actnum);

    void prune_global_for_schedule_run();

    void apply_numerical_aquifers(const NumericalAquifers& numerical_aquifers);

    const std::string& default_region() const;

    std::vector<int> actnum();
    const std::vector<int>& actnumRaw() const;

    template <typename T>
    static bool supported(const std::string& keyword);

    template <typename T>
    bool has(const std::string& keyword) const;

    template <typename T>
    std::vector<std::string> keys() const;

    /// Request read-only property array from internal cache
    ///
    /// Will create property array if permitted and possible.
    ///
    /// \tparam T Property element type.  Typically \c double or \c int.
    ///
    /// \param[in] keyword Property name
    ///
    /// \param[in] flags Options to limit or relax request processing.
    ///   Treated as a bitwise mask of \c TryGetFlags.
    ///
    /// \return Access structure for requested property array.
    template <typename T>
    FieldDataManager<T>
    try_get(const std::string& keyword, const unsigned int flags = 0u)
    {
        const auto allow_unsupported =
            (flags & TryGetFlags::AllowUnsupported) != 0u;

        if (!allow_unsupported && !FieldProps::template supported<T>(keyword)) {
            return { keyword, GetStatus::NOT_SUPPPORTED_KEYWORD, nullptr };
        }

        const auto has0 = this->template has<T>(keyword);
        if (!has0 && ((flags & TryGetFlags::MustExist) != 0)) {
            // Client requested a property which must exist, e.g., as a
            // source array for a COPY operation, but the property has not
            // (yet) been defined in the run's input.
            return { keyword, GetStatus::MISSING_KEYWORD, nullptr };
        }

        const auto& field_data = this->template
            init_get<T>(keyword, std::is_same_v<T, double> && allow_unsupported);

        if (field_data.valid() || allow_unsupported) {
            // Note: FieldDataManager depends on init_get<>() producing a
            // long-lived FieldData instance.
            return { keyword, GetStatus::OK, &field_data };
        }

        if (! has0) {
            // Client requested a property which did not exist and which
            // could not be created from a default description.
            this->template erase<T>(keyword);

            return { keyword, GetStatus::MISSING_KEYWORD, nullptr };
        }

        // If we get here then the property exists but has not been fully
        // defined yet.
        return { keyword, GetStatus::INVALID_DATA, nullptr };
    }

    template <typename T>
    const std::vector<T>& get(const std::string& keyword)
    {
        return this->template try_get<T>(keyword).data();
    }

    template <typename T>
    std::vector<T> get_global(const std::string& keyword)
    {
        const auto managed_field_data = this->template try_get<T>(keyword);
        const auto& field_data = managed_field_data.field_data();

        const auto& kw_info = Fieldprops::keywords::
            template global_kw_info<T>(keyword);

        return kw_info.global
            ? *field_data.global_data
            : this->global_copy(field_data.data, kw_info.scalar_init);
    }

    template <typename T>
    std::vector<T> get_copy(const std::string& keyword, bool global)
    {
        const auto has0 = this->template has<T>(keyword);

        // Recall: FieldDataManager::field_data() will throw various
        // exception types if the 'status' is anything other than 'OK'.
        //
        // Get_copy() depends on this behaviour to not proceed to extracting
        // values in such cases.  In other words, get_copy() uses exceptions
        // for control flow, and we cannot move this try_get() call into the
        // 'has0' branch even though the actual 'field_data' object returned
        // from try_get() is only needed/used there.
        const auto& field_data = this->template try_get<T>(keyword).field_data();

        if (has0) {
            return this->get_copy(field_data.data, field_data.kw_info.scalar_init, global);
        }

        const auto initial_value = Fieldprops::keywords::
            template global_kw_info<T>(keyword).scalar_init;

        return this->get_copy(this->template extract<T>(keyword), initial_value, global);
    }

    template <typename T>
    std::vector<bool> defaulted(const std::string& keyword)
    {
        const auto& field = this->template init_get<T>(keyword);
        std::vector<bool> def(field.numCells());

        for (std::size_t i = 0; i < def.size(); ++i) {
            def[i] = value::defaulted(field.value_status[i]);
        }

        return def;
    }

    template <typename T>
    std::vector<T> global_copy(const std::vector<T>&   data,
                               const std::optional<T>& default_value) const
    {
        const T fill_value = default_value.has_value() ? *default_value : 0;

        std::vector<T> global_data(this->global_size, fill_value);

        std::size_t i = 0;
        for (std::size_t g = 0; g < this->global_size; g++) {
            if (this->m_actnum[g]) {
                global_data[g] = data[i];
                ++i;
            }
        }

        return global_data;
    }

    std::size_t active_size;
    std::size_t global_size;

    std::size_t num_int() const
    {
        return this->int_data.size();
    }

    std::size_t num_double() const
    {
        return this->double_data.size();
    }

    void handle_schedule_keywords(const std::vector<DeckKeyword>& keywords);
    bool tran_active(const std::string& keyword) const;
    void apply_tran(const std::string& keyword, std::vector<double>& data);
    void apply_tranz_global(const std::vector<size_t>& indices, std::vector<double>& data) const;
    bool operator==(const FieldProps& other) const;
    static bool rst_cmp(const FieldProps& full_arg, const FieldProps& rst_arg);

    const std::unordered_map<std::string,Fieldprops::TranCalculator>& getTran() const
    {
        return tran;
    }

    std::vector<std::string> fip_regions() const;

    void deleteMINPVV();

    void set_active_indices(const std::vector<int>& indices);

private:
    void processMULTREGP(const Deck& deck);
    void scanGRIDSection(const GRIDSection& grid_section);
    void scanGRIDSectionOnlyACTNUM(const GRIDSection& grid_section);
    void scanEDITSection(const EDITSection& edit_section);
    void scanPROPSSection(const PROPSSection& props_section);
    void scanREGIONSSection(const REGIONSSection& regions_section);
    void scanSOLUTIONSection(const SOLUTIONSection& solution_section, const std::size_t ncomps);
    double getSIValue(const std::string& keyword, double raw_value) const;
    double getSIValue(ScalarOperation op, const std::string& keyword, double raw_value) const;

    template <typename T>
    void erase(const std::string& keyword);

    template <typename T>
    std::vector<T> extract(const std::string& keyword);

    template <typename T>
    std::vector<T> get_copy(const std::vector<T>&   x,
                            const std::optional<T>& initial_value,
                            const bool              global) const
    {
        return (! global) ? x : this->global_copy(x, initial_value);
    }

    template <typename T>
    std::vector<T> get_copy(std::vector<T>&&        x,
                            const std::optional<T>& initial_value,
                            const bool              global) const
    {
        return (! global) ? std::move(x) : this->global_copy(x, initial_value);
    }

    template <typename T>
    void operate(const DeckRecord& record,
                 Fieldprops::FieldData<T>& target_data,
                 const Fieldprops::FieldData<T>& src_data,
                 const std::vector<Box::cell_index>& index_list,
                 const bool global = false);

    template <typename T>
    Fieldprops::FieldData<T>&
    init_get(const std::string& keyword, bool allow_unsupported = false);

    template <typename T>
    Fieldprops::FieldData<T>&
    init_get(const std::string& keyword,
             const Fieldprops::keywords::keyword_info<T>& kw_info,
             const bool multiplier_in_edit = false);

    std::string region_name(const DeckItem& region_item) const;

    std::pair<std::vector<Box::cell_index>,bool>
    region_index(const std::string& region_name, int region_value);

    void handle_OPERATE(const DeckKeyword& keyword, Box box);
    void handle_operation(Section section, const DeckKeyword& keyword, Box box,
                          const bool onlyACTNUM = false);
    void handle_operateR(const DeckKeyword& keyword);
    void handle_region_operation(const DeckKeyword& keyword);
    void handle_COPY(const DeckKeyword& keyword, Box box, bool region);
    void distribute_toplayer(Fieldprops::FieldData<double>& field_data,
                             const std::vector<double>& deck_data,
                             const Box& box);

    double get_beta(const std::string& func_name, const std::string& target_array, double raw_beta);
    double get_alpha(const std::string& func_name, const std::string& target_array, double raw_alpha);

    void handle_keyword(Section section, const DeckKeyword& keyword, Box& box,
                        const bool onlyACTNUM = false);
    void handle_double_keyword(Section section,
                               const Fieldprops::keywords::keyword_info<double>& kw_info,
                               const DeckKeyword& keyword,
                               const std::string& keyword_name,
                               const Box& box);

    void handle_double_keyword(Section section,
                               const Fieldprops::keywords::keyword_info<double>& kw_info,
                               const DeckKeyword& keyword,
                               const Box& box);

    void handle_int_keyword(const Fieldprops::keywords::keyword_info<int>& kw_info,
                            const DeckKeyword& keyword,
                            const Box& box);

    void init_satfunc(const std::string& keyword, Fieldprops::FieldData<double>& satfunc);
    void init_porv(Fieldprops::FieldData<double>& porv);
    void init_tempi(Fieldprops::FieldData<double>& tempi);

    std::string canonical_fipreg_name(const std::string& fipreg);
    const std::string& canonical_fipreg_name(const std::string& fipreg) const;

    /// \brief Apply multipliers of the EDIT section
    ///
    /// Multipliers are stored intermediately in FieldPropsManager::getMultiplierPrefix()+keyword_name FieldData arrays
    /// to prevent EQUALS MULT* in the EDIT section from overwriting values from the GRID
    /// section. This method will now apply them to keyword_name FieldData arrays and throw away the intermediate storage
    void apply_multipliers();

    static constexpr std::string_view getMultiplierPrefix()
    {
        using namespace std::literals;
        return "__MULT__"sv;
    }

    const UnitSystem unit_system;
    std::size_t nx,ny,nz;
    Phases m_phases;
    SatFuncControls m_satfuncctrl;
    std::vector<int> m_actnum;
    std::unordered_map<int,int> m_active_index;
    std::vector<double> cell_volume;
    std::vector<double> cell_depth;
    const std::string m_default_region;
    const EclipseGrid * grid_ptr;      // A bit undecided whether to properly use the grid or not ...
    TableManager tables;
    std::optional<satfunc::RawTableEndPoints> m_rtep;
    std::vector<MultregpRecord> multregp;
    std::unordered_map<std::string, Fieldprops::FieldData<int>> int_data;
    std::unordered_map<std::string, Fieldprops::FieldData<double>> double_data;
    std::unordered_map<std::string, std::string> fipreg_shortname_translation{};

    std::unordered_map<std::string,Fieldprops::TranCalculator> tran;

    /// \brief A map of multiplier keywords found in the EDIT/SCHEDULE section
    ///
    /// key ist the original keyword name. and value is the keyword information.
    /// This list is used in apply_multipliers, where the multipliers will actually
    /// be applied.
    std::unordered_map<std::string,Fieldprops::keywords::keyword_info<double>> multiplier_kw_infos_;
};

}
#endif
