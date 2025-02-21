/*
  Copyright 2019 Equinor.
  Copyright 2017 Statoil ASA.
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

#include <opm/output/eclipse/Tables.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Phase.hpp>
#include <opm/input/eclipse/EclipseState/Tables/FlatTable.hpp> // PVTW, PVCDO, SWOFLET, SGOFLET
#include <opm/input/eclipse/EclipseState/Tables/GsfTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PvdgTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PvdoTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PvtgTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PvtoTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SgfnTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SgofTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SgwfnTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SlgofTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/Sof2Table.hpp>
#include <opm/input/eclipse/EclipseState/Tables/Sof3Table.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SwfnTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SwofTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/Tabdims.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableContainer.hpp>
#include <opm/input/eclipse/EclipseState/Tables/WsfTable.hpp>

#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <opm/output/eclipse/VectorItems/tabdims.hpp>
#include <opm/output/eclipse/LinearisedOutputTable.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iterator>
#include <numeric>
#include <tuple>
#include <utility>
#include <vector>

using Ix = ::Opm::RestartIO::Helpers::VectorItems::TabDims::index;

namespace {
    /// Convenience type alias for a callable entity that extracts the
    /// independent and primary dependent variates of a single property
    /// function table into a linearised output table.  Calling the function
    /// will assign the independent variate of the sub-table primID within
    /// the table identified as tableID to column zero of 'table' and all
    /// dependent variates to columns one &c of 'table'.  The function call
    /// must return the number of active (used) rows within the sub-table
    /// through its return value.
    using BuildDep = std::function<
        std::size_t(const std::size_t             tableID,
                    const std::size_t             primID,
                    ::Opm::LinearisedOutputTable& table)
        >;

    /// Create linearised, padded TAB vector entries for a collection of
    /// tabulated saturation functions corresponding to a single input
    /// keyword.
    ///
    /// \param[in] numTab Number of tables in this table collection.
    ///
    /// \param[in] numPrim Number of primary look-up keys for each table.
    ///    Mostly relevant (i.e., greater than one) for miscible oil
    ///    ("PVTO") or miscible gas ("PVTG") tables and typically
    ///    corresponds to the number of Rs/Rv entries from the TABDIMS
    ///    keyword.
    ///
    /// \param[in] numRows Number of rows to allocate for each padded
    ///    output table.
    ///
    /// \param[in] numDep Number of dependent variates (columns) in each
    ///    table of this collection of input tables.  Total number of
    ///    columns in the result vector will be 1 + 2*numDep to account
    ///    for the independent variate, the dependent variates and the
    ///    derivatives of the dependent variates with respect to the
    ///    independent variate.
    ///
    /// \param[in] buildDeps Function object that implements the
    ///    protocol outlined for \code BuildDep::operator()()
    ///    \endcode.  Typically initialised from a lambda expression.
    ///
    /// \return Linearised, padded TAB vector entries for a collection of
    ///    tabulated property functions (e.g., saturation functions or PVT
    ///    functions) corresponding to a single input keyword.  Derivatives
    ///    included as additional columns.
    std::vector<double>
    createPropfuncTable(const std::size_t numTab,
                        const std::size_t numPrim,
                        const std::size_t numRows,
                        const std::size_t numDep,
                        const double      fillVal,
                        const BuildDep&   buildDeps)
    {
        const auto numCols = 1 + 2*numDep;

        auto descr = ::Opm::DifferentiateOutputTable::Descriptor{};

        auto linTable = ::Opm::LinearisedOutputTable {
            numTab, numPrim, numRows, numCols, fillVal
        };

        // Note unusual loop counters here.  Specifically, we mutate the
        // members of 'descr' rather than using separate counters, which is
        // the custom construction.  The reason for this is that these
        // members are also used in calcSlopes() to indicate which table
        // subsection to analyse and compute derivatives for.  With separate
        // loop counters, we would need to form a new table descriptor for
        // each call to calcSlopes().
        for (descr.tableID = 0*numTab;
             descr.tableID < 1*numTab; ++descr.tableID)
        {
            for (descr.primID = 0*numPrim;
                 descr.primID < 1*numPrim; ++descr.primID)
            {
                descr.numActRows =
                    buildDeps(descr.tableID, descr.primID, linTable);

                // Derivatives.  Use values already stored in linTable to
                // take advantage of any unit conversion already applied.
                // We don't have to do anything special for the units here.
                //
                // Note: argument 'descr' implies argument-dependent lookup
                //    whence we unambiguously invoke function calcSlopes()
                //    from namespace ::Opm::DifferentiateOutputTable.
                calcSlopes(numDep, descr, linTable);
            }
        }

        return linTable.getDataDestructively();
    }
} // Anonymous

/// Functions to facilitate generating TAB vector entries for tabulated
/// saturation functions.
namespace { namespace SatFunc {
    namespace detail {
        /// Create linearised, padded TAB vector entries for a collection of
        /// tabulated saturation functions corresponding to a single input
        /// keyword.
        ///
        /// \tparam BuildDependent Callable entity that extracts the
        ///    independent and primary dependent variates of a single
        ///    saturation function table into a linearised output table.
        ///    Must implement a function call operator of the form
        ///    \code
        ///       std::size_t
        ///       operator()(const std::size_t      tableID,
        ///                  const std::size_t      primID,
        ///                  LinearisedOutputTable& table);
        ///    \endcode
        ///    that will assign the independent variate of the sub-table
        ///    primID within the table identified as tableID to column zero
        ///    of 'table' and all dependent variates to columns one &c of
        ///    'table'.  The function call operator must return the number
        ///    of active (used) rows within the sub-table through its return
        ///    value.
        ///
        /// \param[in] numTab Number of tables in this table collection.
        ///
        /// \param[in] numRows Number of rows to allocate for each padded
        ///    output table.
        ///
        /// \param[in] numDep Number of dependent variates (columns) in each
        ///    table of this collection of input tables.  Total number of
        ///    columns in the result vector will be 1 + 2*numDep to account
        ///    for the independent variate, the dependent variates and the
        ///    derivatives of the dependent variates with respect to the
        ///    independent variate.
        ///
        /// \param[in] buildDeps Function object that implements the
        ///    protocol outlined for \code BuildDependent::operator()()
        ///    \endcode.  Typically a lambda expression.
        ///
        /// \return Linearised, padded TAB vector entries for a collection
        ///    of tabulated saturation functions corresponding to a single
        ///    input keyword.  Derivatives included as additional columns.
        template <class BuildDependent>
        std::vector<double>
        createSatfuncTable(const std::size_t numTab,
                           const std::size_t numRows,
                           const std::size_t numDep,
                           BuildDependent&&  buildDeps)
        {
            // Saturation functions do not have sub-tables so the number of
            // primary look-up keys is one.
            const auto numPrim = std::size_t{1};
            const auto fillVal = 1.0e20;

            return createPropfuncTable(numTab, numPrim, numRows, numDep, fillVal,
                                       std::forward<BuildDependent>(buildDeps));
        }

        /// Normalise and output relative permeability values to destination range
        ///
        /// \tparam InIt Input iterator type.
        ///
        /// \tparam OutIt Output/destination iterator type.
        ///
        /// \param begin Start of input range of relative permeability values
        ///
        /// \param end End of input range of relative permeability values.
        ///    Must be reachable from \p begin.
        ///
        /// \param tolcrit Minimum relative permeability threshold value for
        ///    phase to be considered mobile.  Values less than this threshold
        ///    are output as zero.
        ///
        /// \param dest Beginning of output range for transformed relative
        ///    permeability values.  Must be able to accommodate at least
        ///    \code distance(begin, end) \endcode output values.
        template <typename InIt, typename OutIt>
        void outputRelperm(InIt         begin,
                           InIt         end,
                           const double tolcrit,
                           OutIt        dest)
        {
            std::transform(begin, end, dest,
                [tolcrit](const double kr) -> double
            {
                return (kr > tolcrit) ? kr : 0.0;
            });
        }

        /// Normalise and output relative permeability values to destination range
        ///
        /// \tparam RelPermColumn Container type for relative permeability
        ///    values.  Typically an \code Opm::TableColumn \endcode.
        ///
        /// \tparam OutIt Output/destination iterator type.
        ///
        /// \param krCol Range of relative permeability values.
        ///
        /// \param tolcrit Minimum relative permeability threshold value for
        ///    phase to be considered mobile.  Values less than this threshold
        ///    are output as zero.
        ///
        /// \param dest Beginning of output range for transformed relative
        ///    permeability values.  Must be able to accommodate at least
        ///    \code distance(begin, end) \endcode output values.
        template <typename RelPermColumn, typename OutIt>
        void outputRelperm(const RelPermColumn& krCol,
                           const double         tolcrit,
                           OutIt&&              dest)
        {
            outputRelperm(std::begin(krCol), std::end(krCol),
                          tolcrit, std::forward<OutIt>(dest));
        }
    } // namespace detail

    /// Quick query structure for the run's active phases
    struct ActivePhases
    {
        /// Constructor
        ///
        /// \param[in] es Run's master configuration object
        explicit ActivePhases(const Opm::EclipseState& es);

        /// Three-phase predicate
        ///
        /// \return Whether or not current run is three-phase
        bool isThreePhase() const
        {
            return this->oil && this->gas && this->wat;
        }

        /// Whether or not OIL is active in the current run
        bool oil{false};

        /// Whether or not GAS is active in the current run
        bool gas{false};

        /// Whether or not WATER is active in the current run
        bool wat{false};
    };

    ActivePhases::ActivePhases(const Opm::EclipseState& es)
    {
        const auto& phases = es.runspec().phases();

        this->oil = phases.active(Opm::Phase::OIL);
        this->gas = phases.active(Opm::Phase::GAS);
        this->wat = phases.active(Opm::Phase::WATER);
    }

    /// Query structure for the run's saturation function family
    class Family
    {
    public:
        /// Constructor
        ///
        /// \param[in] es Run's master configuration object.
        ///
        /// \param[in] active Query object for run's active phases.
        explicit Family(const Opm::EclipseState& es, const ActivePhases& active);

        /// Whether or not run uses saturation functions of family I--i.e.,
        /// SGOF/SWOF/SLGOF/SGOFLET/SWOFLET.
        [[nodiscard]] bool one() const noexcept
        {
            return this->is(Type::One);
        }

        /// Whether or not run uses saturation functions of family II--i.e.,
        /// SGFN/SOF{2,3}/SWFN/SWGFN.
        [[nodiscard]] bool two() const noexcept
        {
            return this->is(Type::Two);
        }

        /// Whether or not run uses saturation functions of family
        /// III--i.e., GSF and WSF.  Two-phase gas/water runs and,
        /// especially, CO2STORE.
        [[nodiscard]] bool three() const noexcept
        {
            return this->is(Type::Three);
        }

        /// Whether or not run uses an ambiguous set of saturation function
        /// tables--i.e., tables from at least two families or no family at
        /// all.
        [[nodiscard]] bool ambiguous() const noexcept
        {
            return ! (this->one() || this->two() || this->three());
        }

    private:
        /// Categories for run's saturation function family.
        enum class Type {
            /// Run uses saturation function family I (SGOF/SWOF &c)
            One,

            /// Run uses saturation function family II (SGFN/SOF*/SWFN &c)
            Two,

            /// Run uses saturation function family III (GSF/WSF)
            Three,

            /// No saturation function family specified (input error)
            None,

            /// Run uses saturation functions from at least two families
            /// (input error).
            Ambiguous,
        };

        /// Run's active saturation function family.
        Type type_{ Type::None };

        /// Predicate for run's active saturation function family
        ///
        /// \param[in] t Category against which we'll check the run's active
        ///    saturation function family.
        ///
        /// \return Whether or not the run's active saturation function
        ///    family is \p t.
        [[nodiscard]] bool is(const Type t) const noexcept
        {
            return this->type_ == t;
        }
    };

    namespace detail {
        /// Helper predicate for whether or not current run uses saturation
        /// functions from family I.
        ///
        /// \param[in] es Run's master configuration object.
        ///
        /// \param[in] active Query object for run's active phases.
        bool isFamOne(const ActivePhases&      active,
                      const Opm::TableManager& tabMgr)
        {
            return
                (active.gas && (tabMgr.hasTables("SGOF") ||
                                !tabMgr.getSgofletTable().empty()))

                ||

                (active.wat && (tabMgr.hasTables("SWOF") ||
                                !tabMgr.getSwofletTable().empty()))
                ;
        }

        /// Helper predicate for whether or not current run uses saturation
        /// functions from family II.
        ///
        /// \param[in] es Run's master configuration object.
        ///
        /// \param[in] active Query object for run's active phases.
        bool isFamTwo(const ActivePhases&      active,
                      const Opm::TableManager& tabMgr)
        {
            return
                (active.gas &&
                 (tabMgr.hasTables("SGFN") || tabMgr.hasTables("SGWFN")))

                ||

                (active.oil &&
                 ((active.isThreePhase() && tabMgr.hasTables("SOF3"))
                  || tabMgr.hasTables("SOF2")))

                ||

                (active.wat && tabMgr.hasTables("SWFN"))
                ;
        }

        /// Helper predicate for whether or not current run uses saturation
        /// functions from family III.
        ///
        /// \param[in] es Run's master configuration object.
        ///
        /// \param[in] active Query object for run's active phases.
        bool isFamThree(const ActivePhases&      active,
                        const Opm::TableManager& tabMgr)
        {
            return
                (active.gas && tabMgr.hasTables("GSF")) &&
                (active.wat && tabMgr.hasTables("WSF"))
                ;
        }
    } // namespace detail

    Family::Family(const Opm::EclipseState& es, const ActivePhases& active)
    {
        const auto& tabMgr = es.getTableManager();

        // SGOF and/or SWOF.
        const auto I = detail::isFamOne(active, tabMgr);

        // SGFN, SOF{2,3}, SWFN, SGWFN.
        const auto II = detail::isFamTwo(active, tabMgr);

        // GSF, WSF.
        const auto III = detail::isFamThree(active, tabMgr);

        if (I + II + III > 1) {
            this->type_ = Type::Ambiguous;
        }
        else if (! (I || II || III)) {
            this->type_ = Type::None;
        }
        else if (I) {
            this->type_ = Type::One;
        }
        else if (II) {
            this->type_ = Type::Two;
        }
        else {
            this->type_ = Type::Three;
        }
    }

    namespace detail {
        /// Create linearised and padded 'TAB' vector entries of normalised
        /// SGFN/SWFN tables for all saturation function regions from Family
        /// One table data.
        ///
        /// \param[in] numRows Number of rows to allocate in the output
        ///    vector for each table.  Expected to be equal to the number of
        ///    declared saturation nodes in the simulation run's TABDIMS
        ///    keyword (Item 3).
        ///
        /// \param tolcrit Minimum relative permeability threshold value for
        ///    phase to be considered mobile.  Values less than this threshold
        ///    are output as zero.
        ///
        /// \param[in] units Active unit system.  Needed to convert SI
        ///    convention capillary pressure values (Pascal) to declared
        ///    conventions of the run specification.
        ///
        /// \param[in] paramLET LET parameters for all saturation regions.
        ///
        /// \return Linearised and padded 'TAB' vector values sampled from
        ///    LET.  Unit-converted capillary pressure.  Derivatives for all
        ///    curves.
        template <typename ParamLET>
        std::vector<double>
        sampleLET(const std::size_t          numRows,
                  const double               tolcrit,
                  const Opm::UnitSystem&     units,
                  const ParamLET&            paramLET,
                  const std::vector<double>& scon)
        {
            std::vector<double> letTab(5 * numRows * paramLET.size(), +2.0e20);

            const auto offset = numRows * paramLET.size();
            const auto uPress = ::Opm::UnitSystem::measure::pressure;

            for (std::size_t k = 0; k < paramLET.size(); ++k) {
                const auto& let = paramLET[k];

                const double sMin = let.s1_residual;
                const double sMax = 1.0 - let.s2_residual - scon[k];
                const double ds = (sMax - sMin) / (numRows - 1);

                for (std::size_t i = 0; i < numRows; ++i) {
                    const double S = sMin + i*ds;

                    double Kr = 0.0;
                    if (S > let.s1_critical) {
                        const double Ss =
                            std::min(1.0, (S - let.s1_critical) / (1.0 - scon[k] - let.s1_critical - let.s2_critical));

                        const double powS   = std::pow(Ss,       let.l1_relperm);
                        const double pow1mS = std::pow(1.0 - Ss, let.t1_relperm);

                        Kr = let.krt1_relperm * powS / (powS + pow1mS*let.e1_relperm);
                        Kr = (Kr > tolcrit) ? Kr : 0.0;
                    }

                    const double Ss =
                        std::min(1.0, (S - let.s1_residual) / (1.0 - scon[k] - let.s1_residual - let.s2_residual));

                    const double powS   = std::pow(Ss,       let.t_pc);
                    const double pow1mS = std::pow(1.0 - Ss, let.l_pc);
                    const double Pc = let.pct_pc +
                        (let.pcir_pc - let.pct_pc) * pow1mS / (pow1mS + powS*let.e_pc);

                    letTab[0*offset + k*numRows + i] = S;
                    letTab[1*offset + k*numRows + i] = Kr;
                    letTab[2*offset + k*numRows + i] = units.from_si(uPress, Pc); // Pc;

                    if (i > 0) {
                        letTab[3*offset + k*numRows + i] =
                            (letTab[1*offset + k*numRows + i] - letTab[1*offset + k*numRows + i - 1]) / ds;

                        letTab[4*offset + k*numRows + i] =
                            (letTab[2*offset + k*numRows + i] - letTab[2*offset + k*numRows + i - 1]) / ds;
                    }
                }
            }

            return letTab;
        }

        /// Create linearised and padded 'TAB' vector entries of normalised
        /// two/three-phase SOFN tables for all saturation function regions
        /// from Family One LET data.
        ///
        /// \param[in] numRows Number of rows to allocate in the output
        ///    vector for each table.  If both SWOFLET and SGOFLET are
        ///    active, it is expected to be 2x the number of declared
        ///    saturation nodes in the simulation run's TABDIMS keyword
        ///    (Item 3), otherwise 1x the number.
        ///
        /// \param tolcrit Minimum relative permeability threshold value
        ///    for phase to be considered mobile.  Values less than this
        ///    threshold are output as zero.
        ///
        /// \param[in] swofLET Collection of SWOFLET parameters for all
        ///    saturation regions (compute krow).
        ///
        /// \param[in] sgofLET Collection of SGOFLET parameters for all
        ///    saturation regions (compute krog).
        ///
        /// \return Linearised and padded 'TAB' vector values for two and
        ///    three-phase SOFN tables.  If both SWOFLET and SGOFLET are
        ///    active, columns represent
        ///
        ///       | s | krow(s) | krog(s) | krow'(s) | krog'(s) |
        ///
        ///    If only one is active, columns are removed accordingly.
        std::vector<double>
        sampleLET(const std::size_t        numRows,
                  const double             tolcrit,
                  const Opm::SwofletTable& swofLET,
                  const Opm::SgofletTable& sgofLET)
        {
            std::size_t nCol = 1;
            std::size_t numTab = 0;
            int dOff = -1;

            if (! swofLET.empty()) {
                numTab = swofLET.size();
                nCol += 2;
                dOff += 1;
            }

            if (! sgofLET.empty()) {
                numTab = sgofLET.size();
                nCol += 2;
                dOff += 1;
            }

            std::vector<double> letTab(nCol * numRows * numTab, +2.0e20);

            const auto offset = numRows*numTab;

            for (std::size_t k = 0; k < numTab; ++k) {
                double sMin = 1.0;
                double sMax = 0.0;
                double Swco = 0.0;

                if (!swofLET.empty()) {
                    Swco = swofLET[k].s1_residual;
                    sMin = std::min(sMin,       swofLET[k].s2_residual);
                    sMax = std::max(sMax, 1.0 - swofLET[k].s1_residual);
                }

                if (!sgofLET.empty()) {
                    sMin = std::min(sMin,       sgofLET[k].s2_residual + Swco);
                    sMax = std::max(sMax, 1.0 - sgofLET[k].s1_residual - Swco);
                }

                const double ds = (sMax - sMin) / (numRows - 1);

                for (std::size_t i = 0; i < numRows; ++i) {
                    const double S = sMin + i*ds;
                    letTab[k*numRows + i] = S;

                    if (! swofLET.empty()) {
                        double Krow = 0.0;

                        const auto& letw = swofLET[k];
                        if (S > letw.s2_critical) {
                            const double Sow =
                                std::min(1.0, (S - letw.s2_critical) / (1.0 - letw.s1_critical - letw.s2_critical));

                            const double powS = std::pow(Sow, letw.l2_relperm);
                            const double pow1mS = std::pow(1.0 - Sow, letw.t2_relperm);

                            Krow = letw.krt2_relperm * powS / (powS + pow1mS*letw.e2_relperm);
                            Krow = (Krow > tolcrit) ? Krow : 0.0;
                        }

                        letTab[offset + k*numRows + i] = Krow;

                        if (i > 0) {
                            letTab[(2 + dOff)*offset + k*numRows + i] =
                                (letTab[offset + k*numRows + i] - letTab[offset + k*numRows + i-1]) / ds;
                        }
                    }

                    if (! sgofLET.empty()) {
                        double Krog = 0.0;
                        const auto& letg = sgofLET[k];
                        if (S > letg.s2_critical) {
                            const double Sog =
                                std::min(1.0, (S - letg.s2_critical) / (1.0 - Swco - letg.s1_critical - letg.s2_critical));

                            const double powS   = std::pow(Sog,       letg.l2_relperm);
                            const double pow1mS = std::pow(1.0 - Sog, letg.t2_relperm);

                            Krog = letg.krt2_relperm * powS / (powS + pow1mS*letg.e2_relperm);
                            Krog = (Krog > tolcrit) ? Krog : 0.0;
                        }

                        letTab[(1 + dOff)*offset + k*numRows + i] = Krog;

                        if (i > 0) {
                            letTab[(3 + dOff)*offset + k*numRows + i] =
                                (letTab[(1 + dOff)*offset + k*numRows + i] -
                                 letTab[(1 + dOff)*offset + k*numRows + i - 1]) / ds;
                        }
                    }
                }
            }

            return letTab;
        }

        /// Infer number of allocated rows in normalised saturation function
        /// tables.
        ///
        /// \param[in] nssfunc Run's declared number of saturation nodes (==
        ///    TABDIMS(3)).
        ///
        /// \param[in] tables Tabulated saturation functions corresponding
        ///    to a single input keyword.
        ///
        /// \return Maximum of declared number of saturation nodes and the
        ///    maximum number of active table rows across all saturation
        ///    regions in \p tables.
        std::size_t allocatedRows(const std::size_t          nssfun,
                                  const Opm::TableContainer& tables)
        {
            const std::size_t maxActiveRows =
                std::accumulate(tables.tables().begin(), tables.tables().end(), 0UL,
                                [](const std::size_t acc, const auto& table)
                                {
                                    return std::max(acc, table.second->numRows());
                                });

            return std::max(nssfun, maxActiveRows);
        }
    } // namespace detail

    // -----------------------------------------------------------------------
    // Return type from normalised*() and specific keyword functions:
    //
    //   - Element at index 0 (get<0>(t)) is the number of allocated rows in
    //     each individual table.  This is usually NSSFUN (== TABDIMS(3)),
    //     but we allow input tables larger than this value.  In that case,
    //     the number of allocated rows will correspond to the maximum
    //     number of active table rows in the keyword's input data across
    //     all saturation regions.
    //
    //   - Element at index 1 (get<1>(t)) is the number of individual
    //     tables--i.e., the number of saturation regions.  This is *mostly*
    //     redundant information but provided as a simplification for client
    //     code.
    //
    //   - Element at index 2 (get<2>(t)) is normalised and linearised table
    //     data for the corresponding saturation functions.
    // -----------------------------------------------------------------------

    /// Functions to create linearised, padded, and normalised SGFN output
    /// tables from various input saturation function keywords.
    namespace Gas {
        /// Create linearised and padded 'TAB' vector entries of normalised
        /// SGFN tables for all saturation function regions from Family Two
        /// table data (SGFN keyword).
        ///
        /// \param[in] numRows Number of rows to allocate in the output
        ///    vector for each table.  Expected to be equal to the number of
        ///    declared saturation nodes in the simulation run's TABDIMS
        ///    keyword (Item 3).
        ///
        /// \param tolcrit Minimum relative permeability threshold value for
        ///    phase to be considered mobile.  Values less than this threshold
        ///    are output as zero.
        ///
        /// \param[in] units Active unit system.  Needed to convert SI
        ///    convention capillary pressure values (Pascal) to declared
        ///    conventions of the run specification.
        ///
        /// \param[in] sgfn Collection of SGFN tables for all saturation
        ///    regions.
        ///
        /// \return Descriptor and linearised/normalised/padded 'TAB' vector
        ///    values for output SGFN tables.  A unit-converted copy of the
        ///    input table \p sgfn with added derivatives.
        std::tuple<int, int, std::vector<double>>
        fromSGFN(const std::size_t          numRows,
                 const double               tolcrit,
                 const Opm::UnitSystem&     units,
                 const Opm::TableContainer& sgfn)
        {
            using SGFN = ::Opm::SgfnTable;

            const auto tableSize = detail::allocatedRows(numRows, sgfn);

            const auto numTab = sgfn.size();
            const auto numDep = std::size_t{2}; // Krg, Pcgo

            return {
                static_cast<int>(tableSize),
                static_cast<int>(numTab),
                detail::createSatfuncTable(numTab, tableSize, numDep,
                    [tolcrit, &units, &sgfn](const std::size_t           tableID,
                                             const std::size_t           primID,
                                             Opm::LinearisedOutputTable& linTable)
                {
                    const auto& t = sgfn.getTable<SGFN>(tableID);

                    auto numActRows = std::size_t{0};

                    // Sg
                    {
                        const auto& Sg = t.getSgColumn();

                        numActRows = Sg.size();
                        std::copy(std::begin(Sg), std::end(Sg),
                                  linTable.column(tableID, primID, 0));
                    }

                    // Krg(Sg)
                    detail::outputRelperm(t.getKrgColumn(), tolcrit,
                                          linTable.column(tableID, primID, 1));

                    // Pcgo(Sg)
                    {
                        constexpr auto uPress = ::Opm::UnitSystem::measure::pressure;

                        const auto& pc = t.getPcogColumn();
                        std::transform(std::begin(pc), std::end(pc),
                                       linTable.column(tableID, primID, 2),
                                       [&units](const double Pc) -> double
                                       {
                                           return units.from_si(uPress, Pc);
                                       });
                    }

                    // Inform createSatfuncTable() of number of active rows
                    // in this table.  Needed to compute slopes of piecewise
                    // linear interpolants.
                    return numActRows;
                })
            };
        }

        /// Create linearised and padded 'TAB' vector entries of normalised
        /// SGFN tables for all saturation function regions from Family Two
        /// table data (SGWFN keyword).
        ///
        /// \param[in] numRows Number of rows to allocate in the output
        ///    vector for each table.  Expected to be equal to the number of
        ///    declared saturation nodes in the simulation run's TABDIMS
        ///    keyword (Item 3).
        ///
        /// \param tolcrit Minimum relative permeability threshold value for
        ///    phase to be considered mobile.  Values less than this threshold
        ///    are output as zero.
        ///
        /// \param[in] units Active unit system.  Needed to convert SI
        ///    convention capillary pressure values (Pascal) to declared
        ///    conventions of the run specification.
        ///
        /// \param[in] sgwfn Collection of SGWFN tables for all saturation
        ///    regions.
        ///
        /// \return Descriptor and linearised/normalised/padded 'TAB' vector
        ///    values for output SGFN tables.  Corresponds to unit-converted
        ///    copies of columns 1, 2, and 4--with added derivatives--of the
        ///    input SGWFN tables.
        std::tuple<int, int, std::vector<double>>
        fromSGWFN(const std::size_t          numRows,
                  const double               tolcrit,
                  const Opm::UnitSystem&     units,
                  const Opm::TableContainer& sgwfn)
        {
            using Sgwfn = ::Opm::SgwfnTable;

            const auto tableSize = detail::allocatedRows(numRows, sgwfn);

            const auto numTab = sgwfn.size();
            const auto numDep = std::size_t{2}; // Krg, Pcgw

            return {
                static_cast<int>(tableSize),
                static_cast<int>(numTab),
                detail::createSatfuncTable(numTab, tableSize, numDep,
                    [tolcrit, &units, &sgwfn](const std::size_t           tableID,
                                              const std::size_t           primID,
                                              Opm::LinearisedOutputTable& linTable)
                {
                    const auto& t = sgwfn.getTable<Sgwfn>(tableID);

                    auto numActRows = std::size_t{0};

                    // Sg
                    {
                        const auto& Sg = t.getSgColumn();

                        numActRows = Sg.size();
                        std::copy(std::begin(Sg), std::end(Sg),
                                  linTable.column(tableID, primID, 0));
                    }

                    // Krg(Sg)
                    detail::outputRelperm(t.getKrgColumn(), tolcrit,
                                          linTable.column(tableID, primID, 1));

                    // Pcgw(Sg)
                    {
                        constexpr auto uPress = ::Opm::UnitSystem::measure::pressure;

                        const auto& pc = t.getPcgwColumn();
                        std::transform(std::begin(pc), std::end(pc),
                                       linTable.column(tableID, primID, 2),
                                       [&units](const double Pc) -> double
                                       {
                                           return units.from_si(uPress, Pc);
                                       });
                    }

                    // Inform createSatfuncTable() of number of active rows
                    // in this table.  Needed to compute slopes of piecewise
                    // linear interpolants.
                    return numActRows;
                })
            };
        }

        /// Create linearised and padded 'TAB' vector entries of normalised
        /// SGFN tables for all saturation function regions from Family Two
        /// table data.
        ///
        /// Distinguishes between keywords SGFN and SGWFN.
        ///
        /// \param[in] numRows Number of rows to allocate in the output
        ///    vector for each table.  Expected to be equal to the number of
        ///    declared saturation nodes in the simulation run's TABDIMS
        ///    keyword (Item 3).
        ///
        /// \param tolcrit Minimum relative permeability threshold value for
        ///    phase to be considered mobile.  Values less than this threshold
        ///    are output as zero.
        ///
        /// \param[in] units Active unit system.  Needed to convert SI
        ///    convention capillary pressure values (Pascal) to declared
        ///    conventions of the run specification.
        ///
        /// \param[in] tabMgr Collection of run's tabulated functions.
        ///
        /// \return Descriptor and linearised/normalised/padded 'TAB' vector
        ///    values for output SGFN tables.
        std::tuple<int, int, std::vector<double>>
        normalisedTable_Two(const std::size_t        numRows,
                            const double             tolcrit,
                            const Opm::UnitSystem&   units,
                            const Opm::TableManager& tabMgr)
        {
            if (const auto& sgfn = tabMgr.getSgfnTables(); ! sgfn.empty()) {
                // Typical case: Saturation functions for gas in SGFN table
                return fromSGFN(numRows, tolcrit, units, sgfn);
            }

            // Saturation functions for gas given in SGWFN tables
            return fromSGWFN(numRows, tolcrit, units, tabMgr.getSgwfnTables());
        }

        /// Create linearised and padded 'TAB' vector entries of normalised
        /// SGFN tables for all saturation function regions from Family One
        /// table data (SGOF keyword).
        ///
        /// \param[in] numRows Number of rows to allocate in the output
        ///    vector for each table.  Expected to be equal to the number of
        ///    declared saturation nodes in the simulation run's TABDIMS
        ///    keyword (Item 3).
        ///
        /// \param tolcrit Minimum relative permeability threshold value for
        ///    phase to be considered mobile.  Values less than this threshold
        ///    are output as zero.
        ///
        /// \param[in] units Active unit system.  Needed to convert SI
        ///    convention capillary pressure values (Pascal) to declared
        ///    conventions of the run specification.
        ///
        /// \param[in] sgof Collection of SGOF tables for all saturation
        ///    regions.
        ///
        /// \return Descriptor and linearised/normalised/padded 'TAB' vector
        ///    values for output SGFN tables.  Corresponds to unit-converted
        ///    copies of columns 1, 2, and 4--with added derivatives--of the
        ///    input SGOF tables.
        std::tuple<int, int, std::vector<double>>
        fromSGOF(const std::size_t          numRows,
                 const double               tolcrit,
                 const Opm::UnitSystem&     units,
                 const Opm::TableContainer& sgof)
        {
            using SGOF = ::Opm::SgofTable;

            const auto tableSize = detail::allocatedRows(numRows, sgof);

            const auto numTab = sgof.size();
            const auto numDep = std::size_t{2}; // Krg, Pcgo

            return {
                static_cast<int>(tableSize),
                static_cast<int>(numTab),
                detail::createSatfuncTable(numTab, tableSize, numDep,
                    [tolcrit, &units, &sgof](const std::size_t           tableID,
                                             const std::size_t           primID,
                                             Opm::LinearisedOutputTable& linTable)
                {
                    const auto& t = sgof.getTable<SGOF>(tableID);

                    auto numActRows = std::size_t{0};

                    // Sg
                    {
                        const auto& Sg = t.getSgColumn();

                        numActRows = Sg.size();
                        std::copy(std::begin(Sg), std::end(Sg),
                                  linTable.column(tableID, primID, 0));
                    }

                    // Krg(Sg)
                    detail::outputRelperm(t.getKrgColumn(), tolcrit,
                                          linTable.column(tableID, primID, 1));

                    // Pcgo(Sg)
                    {
                        constexpr auto uPress = ::Opm::UnitSystem::measure::pressure;

                        const auto& pc = t.getPcogColumn();
                        std::transform(std::begin(pc), std::end(pc),
                                       linTable.column(tableID, primID, 2),
                                       [&units](const double Pc) -> double
                                       {
                                           return units.from_si(uPress, Pc);
                                       });
                    }

                    // Inform createSatfuncTable() of number of active rows
                    // in this table.  Needed to compute slopes of piecewise
                    // linear interpolants.
                    return numActRows;
                })
            };
        }

        /// Create linearised and padded 'TAB' vector entries of normalised
        /// SGFN tables for all saturation function regions from Family One
        /// table data (SLGOF keyword).
        ///
        /// \param[in] numRows Number of rows to allocate in the output
        ///    vector for each table.  Expected to be equal to the number of
        ///    declared saturation nodes in the simulation run's TABDIMS
        ///    keyword (Item 3).
        ///
        /// \param tolcrit Minimum relative permeability threshold value for
        ///    phase to be considered mobile.  Values less than this threshold
        ///    are output as zero.
        ///
        /// \param[in] units Active unit system.  Needed to convert SI
        ///    convention capillary pressure values (Pascal) to declared
        ///    conventions of the run specification.
        ///
        /// \param[in] slgof Collection of SLGOF tables for all saturation
        ///    regions.
        ///
        /// \return Descriptor and linearised/normalised/padded 'TAB' vector
        ///    values for output SGFN tables.  Corresponds to unit-converted
        ///    copies of columns 1, 2, and 4--with added derivatives--of the
        ///    input SLGOF tables.
        std::tuple<int, int, std::vector<double>>
        fromSLGOF(const std::size_t          numRows,
                 const double               tolcrit,
                 const Opm::UnitSystem&     units,
                 const Opm::TableContainer& slgof)
        {
            using SLGOF = ::Opm::SlgofTable;

            const auto tableSize = detail::allocatedRows(numRows, slgof);

            const auto numTab = slgof.size();
            const auto numDep = std::size_t{2}; // Krg, Pcgo

            return {
                static_cast<int>(tableSize),
                static_cast<int>(numTab),
                detail::createSatfuncTable(numTab, tableSize, numDep,
                    [tolcrit, &units, &slgof](const std::size_t           tableID,
                                                    const std::size_t           primID,
                                                    Opm::LinearisedOutputTable& linTable)
                {
                    const auto& t = slgof.getTable<SLGOF>(tableID);

                    auto numActRows = std::size_t{0};

                    /// Reverse order in all columns, since Sl is increasing (Sg descreasing)
                    // Sg = 1.0 - Sl
                    {
                        const auto& Sl = t.getSlColumn();

                        numActRows = Sl.size();
                        std::vector<double> Sg(std::rbegin(Sl), std::rend(Sl));
                        std::transform(Sg.begin(), Sg.end(), Sg.begin(),
                                       [](double x){ return (1.0 - x); });
                        std::copy(std::begin(Sg), std::end(Sg),
                                  linTable.column(tableID, primID, 0));
                    }

                    // Krg(Sg)
                    const auto& Krg_rev = t.getKrgColumn();
                    std::vector<double> Krg(Krg_rev.rbegin(), Krg_rev.rend());
                    detail::outputRelperm(Krg, tolcrit,
                                          linTable.column(tableID, primID, 1));

                    // Pcgo(Sg)
                    {
                        constexpr auto uPress = ::Opm::UnitSystem::measure::pressure;

                        const auto& pc = t.getPcogColumn();
                        std::transform(std::rbegin(pc), std::rend(pc),
                                       linTable.column(tableID, primID, 2),
                                       [&units](const double Pc) -> double
                                       {
                                           return units.from_si(uPress, Pc);
                                       });
                    }

                    // Inform createSatfuncTable() of number of active rows
                    // in this table.  Needed to compute slopes of piecewise
                    // linear interpolants.
                    return numActRows;
                })
            };
        }

        /// Create linearised and padded 'TAB' vector entries of normalised
        /// SGFN tables for all saturation function regions from Family One
        /// table data.
        ///
        /// Distinguishes between keywords SGOF, SLGOF and SGOFLET.
        ///
        /// \param[in] numRows Number of rows to allocate in the output
        ///    vector for each table.  Expected to be equal to the number of
        ///    declared saturation nodes in the simulation run's TABDIMS
        ///    keyword (Item 3).
        ///
        /// \param tolcrit Minimum relative permeability threshold value for
        ///    phase to be considered mobile.  Values less than this threshold
        ///    are output as zero.
        ///
        /// \param[in] units Active unit system.  Needed to convert SI
        ///    convention capillary pressure values (Pascal) to declared
        ///    conventions of the run specification.
        ///
        /// \param[in] tabMgr Collection of run's tabulated functions.
        ///
        /// \return Descriptor and linearised/normalised/padded 'TAB' vector
        ///    values for output SGFN tables.
        std::tuple<int, int, std::vector<double>>
        normalisedTable_One(const std::size_t        numRows,
                            const double             tolcrit,
                            const Opm::UnitSystem&   units,
                            const Opm::TableManager& tabMgr)
        {
            const auto& letTables = tabMgr.getSgofletTable();

            if (letTables.empty()) {
                // Typical case: Not LET-based saturation functions.
                const auto& sgof = tabMgr.getSgofTables();
                if (!sgof.empty()) {
                    return fromSGOF(numRows, tolcrit, units, sgof);
                }
                const auto& slgof = tabMgr.getSlgofTables();
                if (!slgof.empty()) {
                    return fromSLGOF(numRows, tolcrit, units, slgof);
                }
                // Warn here?
            }

            // LET-based saturation functions for gas.
            const auto swco = [ntab = letTables.size(), &swofLET = tabMgr.getSwofletTable()]()
            {
                auto connsat = std::vector<double>(ntab, 0.0);

                std::transform(swofLET.begin(),
                               swofLET.end(),
                               connsat.begin(),
                               [](const auto& letRecord)
                               { return letRecord.s1_residual; });

                return connsat;
            }();

            return {
                static_cast<int>(numRows),
                static_cast<int>(letTables.size()),
                detail::sampleLET(numRows, tolcrit, units, letTables, swco)
            };
        }

        /// Create linearised and padded 'TAB' vector entries of normalised
        /// SGFN tables for all saturation function regions from Family
        /// Three table data (GSF keyword).
        ///
        /// \param[in] numRows Number of rows to allocate in the output
        ///    vector for each table.  Expected to be equal to the number of
        ///    declared saturation nodes in the simulation run's TABDIMS
        ///    keyword (Item 3).
        ///
        /// \param tolcrit Minimum relative permeability threshold value for
        ///    phase to be considered mobile.  Values less than this
        ///    threshold are output as zero.
        ///
        /// \param[in] units Active unit system.  Needed to convert SI
        ///    convention capillary pressure values (Pascal) to declared
        ///    conventions of the run specification.
        ///
        /// \param[in] gsf Collection of GSF tables for all saturation
        ///    regions.
        ///
        /// \return Descriptor and linearised/normalised/padded 'TAB' vector
        ///    values for output SGFN tables.  A unit-converted copy of the
        ///    input table \p gsf with added derivatives.
        std::tuple<int, int, std::vector<double>>
        fromGSF(const std::size_t          numRows,
                const double               tolcrit,
                const Opm::UnitSystem&     units,
                const Opm::TableContainer& gsf)
        {
            using Gsf = ::Opm::GsfTable;

            const auto tableSize = detail::allocatedRows(numRows, gsf);

            const auto numTab = gsf.size();
            const auto numDep = std::size_t{2}; // Krg, Pcgw

            return {
                static_cast<int>(tableSize),
                static_cast<int>(numTab),
                detail::createSatfuncTable(numTab, tableSize, numDep,
                    [tolcrit, &units, &gsf](const std::size_t           tableID,
                                            const std::size_t           primID,
                                            Opm::LinearisedOutputTable& linTable)
                {
                    const auto& t = gsf.getTable<Gsf>(tableID);

                    auto numActRows = std::size_t{0};

                    // Sg
                    {
                        const auto& Sg = t.getSgColumn();

                        numActRows = Sg.size();
                        std::copy(std::begin(Sg), std::end(Sg),
                                  linTable.column(tableID, primID, 0));
                    }

                    // Krg(Sg)
                    detail::outputRelperm(t.getKrgColumn(), tolcrit,
                                          linTable.column(tableID, primID, 1));

                    // Pcgw(Sg)
                    {
                        constexpr auto uPress = ::Opm::UnitSystem::measure::pressure;

                        const auto& pc = t.getPcgwColumn();
                        std::transform(std::begin(pc), std::end(pc),
                                       linTable.column(tableID, primID, 2),
                                       [&units](const double Pc) -> double
                                       {
                                           return units.from_si(uPress, Pc);
                                       });
                    }

                    // Inform createSatfuncTable() of number of active rows
                    // in this table.  Needed to compute slopes of piecewise
                    // linear interpolants.
                    return numActRows;
                })
            };
        }

        /// Create linearised and padded 'TAB' vector entries of normalised
        /// SGFN tables for all saturation function regions from table data
        /// for all saturation function families.
        ///
        /// \param[in] family Predicate for current run's saturation
        ///    function family.
        ///
        /// \param[in] numRows Number of rows to allocate in the output
        ///    vector for each table.  Expected to be equal to the number of
        ///    declared saturation nodes in the simulation run's TABDIMS
        ///    keyword (Item 3).
        ///
        /// \param tolcrit Minimum relative permeability threshold value for
        ///    phase to be considered mobile.  Values less than this threshold
        ///    are output as zero.
        ///
        /// \param[in] units Active unit system.  Needed to convert SI
        ///    convention capillary pressure values (Pascal) to declared
        ///    conventions of the run specification.
        ///
        /// \param[in] tabMgr Collection of run's tabulated functions.
        ///
        /// \return Descriptor and linearised/normalised/padded 'TAB' vector
        ///    values for output SGFN tables.
        std::tuple<int, int, std::vector<double>>
        normalisedTable(const Family&            family,
                        const std::size_t        numRows,
                        const double             tolcrit,
                        const Opm::UnitSystem&   units,
                        const Opm::TableManager& tabMgr)
        {
            if (family.one()) {
                return normalisedTable_One(numRows, tolcrit, units, tabMgr);
            }

            if (family.two()) {
                return normalisedTable_Two(numRows, tolcrit, units, tabMgr);
            }

            return fromGSF(numRows, tolcrit, units, tabMgr.getGsfTables());
        }
    } // Gas

    /// Functions to create linearised, padded, and normalised SOFN output
    /// tables from various input saturation function keywords, depending on
    /// number of active phases.
    namespace Oil {
        /// Form normalised SOFN output tables for two-phase runs.
        namespace TwoPhase {
            /// Create linearised and padded 'TAB' vector entries of
            /// normalised two-phase SOFN tables for all saturation function
            /// regions from Family Two table data (SOF2 keyword).
            ///
            /// \param[in] numRows Number of rows to allocate in the output
            ///    vector for each table.  Expected to be equal to the
            ///    number of declared saturation nodes in the simulation
            ///    run's TABDIMS keyword (Item 3).
            ///
            /// \param tolcrit Minimum relative permeability threshold value
            ///    for phase to be considered mobile.  Values less than this
            ///    threshold are output as zero.
            ///
            /// \param[in] sof2 Collection of SOF2 tables for all saturation
            ///   regions.
            ///
            /// \return Descriptor and linearised/normalised/padded 'TAB'
            ///    vector values for output SOFN tables.  Essentially just a
            ///    padded copy of the input SOF2 table--with added
            ///    derivatives.
            std::tuple<int, int, std::vector<double>>
            fromSOF2(const std::size_t          numRows,
                     const double               tolcrit,
                     const Opm::TableContainer& sof2)
            {
                using SOF2 = ::Opm::Sof2Table;

                const auto tableSize = detail::allocatedRows(numRows, sof2);

                const auto numTab = sof2.size();
                const auto numDep = std::size_t{1}; // Kro

                return {
                    static_cast<int>(tableSize),
                    static_cast<int>(numTab),
                    detail::createSatfuncTable(numTab, tableSize, numDep,
                        [tolcrit, &sof2](const std::size_t           tableID,
                                         const std::size_t           primID,
                                         Opm::LinearisedOutputTable& linTable)
                    {
                        const auto& t = sof2.getTable<SOF2>(tableID);

                        auto numActRows = std::size_t{0};

                        // So
                        {
                            const auto& So = t.getSoColumn();

                            numActRows = So.size();
                            std::copy(std::begin(So), std::end(So),
                                      linTable.column(tableID, primID, 0));
                        }

                        // Kro(So)
                        detail::outputRelperm(t.getKroColumn(), tolcrit,
                                              linTable.column(tableID, primID, 1));

                        // Inform createSatfuncTable() of number of active
                        // rows in this table.  Needed to compute slopes of
                        // piecewise linear interpolants.
                        return numActRows;
                    })
                };
            }

            /// Create linearised and padded 'TAB' vector entries of
            /// normalised two-phase SOFN tables for all saturation function
            /// regions from Family One table data (SGOF keyword--G/O System).
            ///
            /// \param[in] numRows Number of rows to allocate in the output
            ///    vector for each table.  Expected to be equal to the
            ///    number of declared saturation nodes in the simulation
            ///    run's TABDIMS keyword (Item 3).
            ///
            /// \param tolcrit Minimum relative permeability threshold value
            ///    for phase to be considered mobile.  Values less than this
            ///    threshold are output as zero.
            ///
            /// \param[in] sgof Collection of SGOF tables for all saturation
            ///    regions.
            ///
            /// \return Descriptor and linearised/normalised/padded 'TAB'
            ///    vector values for output SOFN tables.  Corresponds to
            ///    translated (1-Sg), reverse saturation column (column 1)
            ///    and reverse column of relative permeability for oil
            ///    (column 3) from the input SGOF table--with added
            ///    derivatives.
            std::tuple<int, int, std::vector<double>>
            fromSGOF(const std::size_t          numRows,
                     const double               tolcrit,
                     const Opm::TableContainer& sgof)
            {
                using SGOF = ::Opm::SgofTable;

                const auto tableSize = detail::allocatedRows(numRows, sgof);

                const auto numTab = sgof.size();
                const auto numDep = std::size_t{1}; // Kro

                return {
                    static_cast<int>(tableSize),
                    static_cast<int>(numTab),
                    detail::createSatfuncTable(numTab, tableSize, numDep,
                        [tolcrit, &sgof](const std::size_t           tableID,
                                         const std::size_t           primID,
                                         Opm::LinearisedOutputTable& linTable)
                    {
                        const auto& t = sgof.getTable<SGOF>(tableID);

                        auto numActRows = std::size_t{0};

                        // So
                        {
                            const auto& Sg = t.getSgColumn();
                            numActRows = Sg.size();

                            auto So = std::vector<double>{};
                            So.reserve(numActRows);

                            // Two-phase system => So = 1-Sg
                            std::transform(std::begin(Sg), std::end(Sg),
                                           std::back_inserter(So),
                                           [](const double sg) -> double
                                           {
                                               return 1.0 - sg;
                                           });

                            std::copy(So.rbegin(), So.rend(),
                                      linTable.column(tableID, primID, 0));
                        }

                        // Kro(So)
                        {
                            const auto& kr = t.getKrogColumn();

                            const auto krog = std::vector<double> {
                                std::begin(kr), std::end(kr)
                            };

                            detail::outputRelperm(krog.rbegin(), krog.rend(), tolcrit,
                                                  linTable.column(tableID, primID, 1));
                        }

                        // Inform createSatfuncTable() of number of active
                        // rows in this table.  Needed to compute slopes of
                        // piecewise linear interpolants.
                        return numActRows;
                    })
                };
            }


            /// Create linearised and padded 'TAB' vector entries of
            /// normalised two-phase SOFN tables for all saturation function
            /// regions from Family One table data (SLGOF keyword--G/O System).
            ///
            /// \param[in] numRows Number of rows to allocate in the output
            ///    vector for each table.  Expected to be equal to the
            ///    number of declared saturation nodes in the simulation
            ///    run's TABDIMS keyword (Item 3).
            ///
            /// \param tolcrit Minimum relative permeability threshold value
            ///    for phase to be considered mobile.  Values less than this
            ///    threshold are output as zero.
            ///
            /// \param[in] slgof Collection of SLGOF tables for all saturation
            ///    regions.
            ///
            /// \return Descriptor and linearised/normalised/padded 'TAB'
            ///    vector values for output SOFN tables.  Corresponds to
            ///    translated (1-Sg), reverse saturation column (column 1)
            ///    and reverse column of relative permeability for oil
            ///    (column 3) from the input SLGOF table--with added
            ///    derivatives.
            std::tuple<int, int, std::vector<double>>
            fromSLGOF(const std::size_t          numRows,
                     const double               tolcrit,
                     const Opm::TableContainer& slgof)
            {
                using SLGOF = ::Opm::SlgofTable;

                const auto tableSize = detail::allocatedRows(numRows, slgof);

                const auto numTab = slgof.size();
                const auto numDep = std::size_t{1}; // Kro

                return {
                    static_cast<int>(tableSize),
                    static_cast<int>(numTab),
                    detail::createSatfuncTable(numTab, tableSize, numDep,
                        [tolcrit, &slgof](const std::size_t           tableID,
                                         const std::size_t           primID,
                                         Opm::LinearisedOutputTable& linTable)
                    {
                        const auto& t = slgof.getTable<SLGOF>(tableID);

                        auto numActRows = std::size_t{0};

                        // So = Sl (two-phase)
                        {
                            const auto& Sl = t.getSlColumn();
                            numActRows = Sl.size();
                            std::copy(Sl.begin(), Sl.end(), linTable.column(tableID, primID, 0));
                        }

                        // Kro(So)
                        {
                            const auto& kr = t.getKrogColumn();

                            const auto krog = std::vector<double> {
                                std::begin(kr), std::end(kr)
                            };

                            detail::outputRelperm(krog.begin(), krog.end(), tolcrit,
                                                  linTable.column(tableID, primID, 1));
                        }

                        // Inform createSatfuncTable() of number of active
                        // rows in this table.  Needed to compute slopes of
                        // piecewise linear interpolants.
                        return numActRows;
                    })
                };
            }

            /// Create linearised and padded 'TAB' vector entries of
            /// normalised SOFN tables for all saturation function regions
            /// from Family One table data in the case of an oil/gas
            /// two-phase system.
            ///
            /// \param[in] numRows Number of rows to allocate in the output
            ///    vector for each table.  Expected to be equal to the
            ///    number of declared saturation nodes in the simulation
            ///    run's TABDIMS keyword (Item 3).
            ///
            /// \param tolcrit Minimum relative permeability threshold value for
            ///    phase to be considered mobile.  Values less than this threshold
            ///    are output as zero.
            ///
            /// \param[in] tabMgr Collection of run's tabulated functions.
            ///
            /// \return Descriptor and linearised/normalised/padded 'TAB' vector
            ///    values for output SOFN tables.
            std::tuple<int, int, std::vector<double>>
            normalisedTable_Gas(const std::size_t        numRows,
                                const double             tolcrit,
                                const Opm::TableManager& tabMgr)
            {
                const auto& gasLET = tabMgr.getSgofletTable();

                if (gasLET.empty()) {
                    const auto& sgof = tabMgr.getSgofTables();
                    if (!sgof.empty()) {
                        return fromSGOF(numRows, tolcrit, sgof);
                    }
                    const auto& slgof = tabMgr.getSlgofTables();
                    if (!slgof.empty()) {
                        return fromSLGOF(numRows, tolcrit, slgof);
                    }
                    // Warn here?
                }

                return {
                    static_cast<int>(numRows),
                    static_cast<int>(gasLET.size()),
                    detail::sampleLET(numRows, tolcrit,
                                      tabMgr.getSwofletTable(),
                                      gasLET)
                };
            }

            /// Create linearised and padded 'TAB' vector entries for
            /// normalised two-phase SOFN tables for all saturation function
            /// regions from Family One table data (SWOF keyword--O/W System).
            ///
            /// \param[in] numRows Number of rows to allocate in the output
            ///    vector for each table.  Expected to be equal to the
            ///    number of declared saturation nodes in the simulation
            ///    run's TABDIMS keyword (Item 3).
            ///
            /// \param tolcrit Minimum relative permeability threshold value
            ///    for phase to be considered mobile.  Values less than this
            ///    threshold are output as zero.
            ///
            /// \param[in] swof Collection of SWOF tables for all saturation
            ///    regions.
            ///
            /// \return Linearised and padded 'TAB' vector values for
            ///    two-phase SOFN tables.  Corresponds to translated (1-Sw),
            ///    reverse saturation column (column 1) and reverse column
            ///    of relative permeability for oil (column 3) from the
            ///    input SWOF table--with added derivatives.
            std::tuple<int, int, std::vector<double>>
            fromSWOF(const std::size_t          numRows,
                     const double               tolcrit,
                     const Opm::TableContainer& swof)
            {
                using SWOF = ::Opm::SwofTable;

                const auto tableSize = detail::allocatedRows(numRows, swof);

                const auto numTab = swof.size();
                const auto numDep = std::size_t{1}; // Kro

                return {
                    static_cast<int>(tableSize),
                    static_cast<int>(numTab),
                    detail::createSatfuncTable(numTab, tableSize, numDep,
                        [tolcrit, &swof](const std::size_t           tableID,
                                         const std::size_t           primID,
                                         Opm::LinearisedOutputTable& linTable)
                    {
                        const auto& t = swof.getTable<SWOF>(tableID);

                        auto numActRows = std::size_t{0};

                        // So
                        {
                            const auto& Sw = t.getSwColumn();
                            numActRows = Sw.size();

                            auto So = std::vector<double>{};
                            So.reserve(numActRows);

                            // Two-phase system => So = 1-Sw
                            std::transform(std::begin(Sw), std::end(Sw),
                                           std::back_inserter(So),
                                           [](const double sw)
                                           {
                                               return 1.0 - sw;
                                           });

                            std::copy(So.rbegin(), So.rend(),
                                      linTable.column(tableID, primID, 0));
                        }

                        // Kro(So)
                        {
                            const auto& kr = t.getKrowColumn();

                            const auto krow = std::vector<double> {
                                std::begin(kr), std::end(kr)
                            };

                            detail::outputRelperm(krow.rbegin(), krow.rend(), tolcrit,
                                                  linTable.column(tableID, primID, 1));
                        }

                        // Inform createSatfuncTable() of number of active
                        // rows in this table.  Needed to compute slopes of
                        // piecewise linear interpolants.
                        return numActRows;
                    })
                };
            }

            /// Create linearised and padded 'TAB' vector entries of
            /// normalised SOFN tables for all saturation function regions
            /// from Family One table data in the case of an oil/water
            /// two-phase system.
            ///
            /// \param[in] numRows Number of rows to allocate in the output
            ///    vector for each table.  Expected to be equal to the
            ///    number of declared saturation nodes in the simulation
            ///    run's TABDIMS keyword (Item 3).
            ///
            /// \param tolcrit Minimum relative permeability threshold value for
            ///    phase to be considered mobile.  Values less than this threshold
            ///    are output as zero.
            ///
            /// \param[in] tabMgr Collection of run's tabulated functions.
            ///
            /// \return Descriptor and linearised/normalised/padded 'TAB' vector
            ///    values for output SOFN tables.
            std::tuple<int, int, std::vector<double>>
            normalisedTable_Water(const std::size_t        numRows,
                                  const double             tolcrit,
                                  const Opm::TableManager& tabMgr)
            {
                const auto& watLET = tabMgr.getSwofletTable();

                if (watLET.empty()) {
                    return fromSWOF(numRows, tolcrit, tabMgr.getSwofTables());
                }

                return {
                    static_cast<int>(numRows),
                    static_cast<int>(watLET.size()),
                    detail::sampleLET(numRows, tolcrit, watLET,
                                      tabMgr.getSgofletTable())
                };
            }

            /// Create linearised and padded 'TAB' vector entries of
            /// normalised SOFN tables for all saturation function regions
            /// from table data for all saturation function families in the
            /// case of a two-phase oil/gas or a two-phase oil/water system.
            ///
            /// \param[in] family Predicate for current run's saturation
            ///    function family.
            ///
            /// \param[in] active Predicate for current run's active phases.
            ///    Needed to distinguish which of the other two phases
            ///    (gas/water) are active.
            ///
            /// \param[in] numRows Number of rows to allocate in the output
            ///    vector for each table.  Expected to be equal to the
            ///    number of declared saturation nodes in the simulation
            ///    run's TABDIMS keyword (Item 3).
            ///
            /// \param tolcrit Minimum relative permeability threshold value for
            ///    phase to be considered mobile.  Values less than this threshold
            ///    are output as zero.
            ///
            /// \param[in] tabMgr Collection of run's tabulated functions.
            ///
            /// \return Descriptor and linearised/normalised/padded 'TAB' vector
            ///    values for output SOFN tables.
            std::tuple<int, int, std::vector<double>>
            normalisedTable(const Family&            family,
                            const ActivePhases&      active,
                            const std::size_t        numRows,
                            const double             tolcrit,
                            const Opm::TableManager& tabMgr)
            {
                if (family.two()) {
                    return fromSOF2(numRows, tolcrit, tabMgr.getSof2Tables());
                }

                // Family I (S[GW]OF or S[GW]OFLET).
                if (active.gas) {
                    // SGOF{,LET}
                    return normalisedTable_Gas(numRows, tolcrit, tabMgr);
                }

                // SWOF{,LET}
                return normalisedTable_Water(numRows, tolcrit, tabMgr);
            }
        } // TwoPhase

        /// Form normalised SOFN output tables for three-phase runs.
        namespace ThreePhase {
            /// Facility to provide oil saturation and relative permeability
            /// look-up based on data in a Family One table.
            class DerivedKroFunction
            {
            public:
                /// Constructor
                ///
                /// \param[in] s Phase saturation values.  Increasing Sg in
                ///   the case of SGOF or increasing Sw in the case of SWOF.
                ///
                /// \param[in] kro Relative permeability for oil.  Should be
                ///   the decreasing KrOG column in the case of SGOF or the
                ///   decreasing KrOW column in the case of SWOF.
                ///
                /// \param[in] So_off Oil saturation offset through which to
                ///   convert input phase saturation values to saturation
                ///   values for oil.  Should be (1 - Sw_conn) in the case
                ///   of SGOF and 1.0 for the case of SWOF.
                DerivedKroFunction(std::vector<double> s,
                                   std::vector<double> kro,
                                   const double        So_off)
                    : s_     (std::move(s))
                    , kro_   (std::move(kro))
                    , So_off_(So_off)
                {}

                /// Get oil saturation at particular node.
                ///
                /// \param[in] i Saturation node identifier.
                ///
                /// \return Oil saturation at node \p i.
                double So(const std::size_t i) const
                {
                    return this->So_off_ - this->s_[i];
                }

                /// Get relative permeability for oil at particular node.
                ///
                /// \param[in] i Saturation node identifier.
                ///
                /// \return Relative permeability for oil at node \p i.
                double Kro(const std::size_t i) const
                {
                    return this->kro_[i];
                }

                /// Get relative permeability for oil at particular oil
                /// saturation.
                ///
                /// Uses piecewise linear interpolation in the input KrO
                /// table.
                ///
                /// \param[in] So Oil saturation.
                ///
                /// \return Relative permeability for oil at So.
                double Kro(const double So) const
                {
                    const auto s = this->So_off_ - So;

                    auto b = std::begin(this->s_);
                    auto e = std::end  (this->s_);
                    auto p = std::lower_bound(b, e, s);

                    if (p == b) { return this->kro_.front(); }
                    if (p == e) { return this->kro_.back (); }

                    const auto i = p - b; // *Right-hand* end-point.

                    const auto sl = this->s_  [i - 1];
                    const auto yl = this->kro_[i - 1];

                    const auto sr = this->s_  [i - 0];
                    const auto yr = this->kro_[i - 0];

                    const auto t = (s - sl) / (sr - sl);

                    return t*yr + (1.0 - t)*yl;
                }

                /// Retrieve number of active saturation nodes in this
                /// table.
                std::vector<double>::size_type size() const
                {
                    return this->s_.size();
                }

            private:
                /// Input phase saturation.  Sg or Sw.
                std::vector<double> s_;

                /// Input relative permeabilty for oil.  KrOG or KrOW.
                std::vector<double> kro_;

                /// Oil saturation offset through which to convert between
                /// input phase saturation and oil saturation.
                double So_off_;
            };

            /// Pair of saturation node index and saturation function table.
            struct TableElement
            {
                /// Which numeric table to use for look-up.
                std::size_t function;

                /// Saturation node ID within 'function'.
                std::size_t index;
            };

            // S{G,W}OF tables have KrOX data in terms of increasing phase
            // saturation for Gas and Water, respectively, so we need to
            // traverse those tables in the opposite direction in order to
            // generate the KrOX values in terms of increasing phase
            // saturation for Oil.
            std::vector<TableElement>
            makeReverseRange(const std::size_t function,
                             const std::size_t n)
            {
                auto ret = std::vector<TableElement>{};
                ret.reserve(n);

                for (auto i = n; i > 0; --i) {
                    ret.push_back( TableElement {
                        function, i - 1
                    });
                }

                return ret;
            }

            // Join derived KrO functions on common saturation values for
            // oil.  Heavy lifting by std::set_union() to avoid outputting
            // common oil saturation values more than once.  Relies on input
            // tables having sorted phase saturation values (required by ECL
            // format).
            std::vector<TableElement>
            mergeTables(const std::vector<DerivedKroFunction>& t, const double tolerance = 1.0e-7)
            {
                auto ret = std::vector<TableElement>{};

                const auto t0 = makeReverseRange(0, t[0].size());
                const auto t1 = makeReverseRange(1, t[1].size());

                ret.reserve(t0.size() + t1.size());

                std::set_union(std::begin(t0), std::end(t0),
                               std::begin(t1), std::end(t1),
                               std::back_inserter(ret),
                    [&t,tolerance](const TableElement& e1, const TableElement& e2)
                {
                    const double val1 = t[e1.function].So(e1.index);
                    const double val2 = t[e2.function].So(e2.index);
                    return ( (val1 + tolerance) < val2);
                });

                return ret;
            }

            // Create collection of individual columns of single SOF3 table
            // through joining input SGOF and SWOF tables on increasing oil
            // saturation and appropriate KrOX columns.
            std::array<std::vector<double>, 3>
            makeSOF3Table(const Opm::SgofTable& sgof,
                          const Opm::SwofTable& swof)
            {
                auto ret = std::array<std::vector<double>, 3>{};

                auto tbl = std::vector<DerivedKroFunction>{};
                tbl.reserve(2);

                // Note: Order between Krow(So) and Krog(So) matters
                // here.  This order must match the expected column
                // order in SOF3--i.e. [ So, Krow, Krog ].

                // 1) Krow(So)
                {
                    // So = 1.0 - Sw
                    const auto& Sw     = swof.getSwColumn();
                    const auto& Krow   = swof.getKrowColumn();
                    const auto& So_off = 1.0;

                    tbl.emplace_back(Sw  .vectorCopy(),
                                     Krow.vectorCopy(), So_off);
                }

                // 2) Krog(So)
                {
                    // So = (1.0 - Sw_conn) - Sg
                    const auto& Sg     = sgof.getSgColumn();
                    const auto& Krog   = sgof.getKrogColumn();
                    const auto  So_off =
                        1.0 - swof.getSwColumn()[0];

                    tbl.emplace_back(Sg  .vectorCopy(),
                                     Krog.vectorCopy(), So_off);
                }

                const auto mrg = mergeTables(tbl);

                for (auto& col : ret) { col.reserve(mrg.size()); }

                for (const auto& row : mrg) {
                    const auto self  =     row.function;
                    const auto other = 1 - row.function;

                    // 1) Assign So
                    ret[0].push_back(tbl[self].So(row.index));

                    // 2) Assign Kro for "self" column (the one that got
                    //    picked for this row).
                    ret[1 + self].push_back(tbl[self].Kro(row.index));

                    // 3) Assign Kro for "other" column (the one that
                    //    did not get picked for this row).
                    ret[1 + other].push_back(tbl[other].Kro(ret[0].back()));
                }

                return ret;
            }

            // Create collection of individual columns of single SOF3 table
            // through joining input SLGOF and SWOF tables on increasing oil
            // saturation and appropriate KrOX columns.
            std::array<std::vector<double>, 3>
            makeSOF3Table(const Opm::SlgofTable& slgof,
                          const Opm::SwofTable& swof)
            {
                auto ret = std::array<std::vector<double>, 3>{};

                auto tbl = std::vector<DerivedKroFunction>{};
                tbl.reserve(2);

                // Note: Order between Krow(So) and Krog(So) matters
                // here.  This order must match the expected column
                // order in SOF3--i.e. [ So, Krow, Krog ].

                // 1) Krow(So)
                {
                    // So = 1.0 - Sw
                    const auto& Sw     = swof.getSwColumn();
                    const auto& Krow   = swof.getKrowColumn();
                    const auto& So_off = 1.0;

                    tbl.emplace_back(Sw  .vectorCopy(),
                                     Krow.vectorCopy(), So_off);
                }

                // 2) Krog(So)
                {
                    // Map Krog(Sl) [decreasing Sg] to Krog(Sg) [increasing Sg]
                    const auto swco = swof.getSwColumn()[0];
                    const auto& Sl     = slgof.getSlColumn();
                    const auto& Krog   = slgof.getKrogColumn();
                    const auto  So_off =
                        1.0 - swco;

                    // Sg = 1.0 - Sl;
                    std::vector<double> sg( Sl.rbegin(), Sl.rend() );
                    std::transform(sg.begin(), sg.end(), sg.begin(),
                                   [](double x){ return (1.0 - x); });

                    std::vector<double> krog ( Krog.rbegin(), Krog.rend() );
                    tbl.emplace_back(sg, krog, So_off);
                }

                const auto mrg = mergeTables(tbl);

                for (auto& col : ret) { col.reserve(mrg.size()); }

                for (const auto& row : mrg) {
                    const auto self  =     row.function;
                    const auto other = 1 - row.function;

                    // 1) Assign So
                    ret[0].push_back(tbl[self].So(row.index));

                    // 2) Assign Kro for "self" column (the one that got
                    //    picked for this row).
                    ret[1 + self].push_back(tbl[self].Kro(row.index));

                    // 3) Assign Kro for "other" column (the one that
                    //    did not get picked for this row).
                    ret[1 + other].push_back(tbl[other].Kro(ret[0].back()));
                }

                return ret;
            }

            /// Create linearised and padded 'TAB' vector entries of
            /// normalised three-phase SOFN tables for all saturation
            /// function regions from Family One table data.
            ///
            /// \param[in] numRows Number of rows to allocate in the output
            ///    vector for each table.  Expected to be twice the number
            ///    of declared saturation nodes in the simulation run's
            ///    TABDIMS keyword (Item 3).
            ///
            /// \param tolcrit Minimum relative permeability threshold value
            ///    for phase to be considered mobile.  Values less than this
            ///    threshold are output as zero.
            ///
            /// \param[in] sgof Collection of SGOF tables for all saturation
            ///    regions.
            ///
            /// \param[in] swof Collection of SWOF tables for all saturation
            ///    regions.
            ///
            /// \return Linearised and padded 'TAB' vector values for
            ///    three-phase SOFN tables.  Corresponds to column 1 from
            ///    both of the input SGOF and SWOF tables, as well as column
            ///    3 from the input SWOF table and column 3 from the input
            ///    SGOF table--expanded so as to have values for all oil
            ///    saturation nodes.  Derivatives added in columns 4 and 5.
            std::tuple<int, int, std::vector<double>>
            fromSGOFandSWOF(const std::size_t          numRows,
                            const double               tolcrit,
                            const Opm::TableContainer& sgof,
                            const Opm::TableContainer& swof)
            {
                using SGOF = ::Opm::SgofTable;
                using SWOF = ::Opm::SwofTable;

                const auto tableSizeGas = detail::allocatedRows(numRows, sgof);
                const auto tableSizeWat = detail::allocatedRows(numRows, swof);

                // Merged table has 2*NSSFUN rows to accommodate both the
                // gas and the water table.
                const auto tableSize =
                    2 * std::max(tableSizeGas, tableSizeWat);

                const auto numTab = sgof.size();
                const auto numDep = std::size_t{2}; // Krow, Krog

                return {
                    static_cast<int>(tableSize),
                    static_cast<int>(numTab),
                    detail::createSatfuncTable(numTab, tableSize, numDep,
                        [tolcrit, &sgof, &swof]
                        (const std::size_t           tableID,
                         const std::size_t           primID,
                         Opm::LinearisedOutputTable& linTable)
                    {
                        const auto sof3 =
                            makeSOF3Table(sgof.getTable<SGOF>(tableID),
                                          swof.getTable<SWOF>(tableID));

                        auto numActRows = std::size_t{0};

                        // So
                        {
                            const auto& So = sof3[0];

                            numActRows = So.size();
                            std::copy(std::begin(So), std::end(So),
                                      linTable.column(tableID, primID, 0));
                        }

                        // Krow(So)
                        detail::outputRelperm(sof3[1], tolcrit,
                                              linTable.column(tableID, primID, 1));

                        // Krog(So)
                        detail::outputRelperm(sof3[2], tolcrit,
                                              linTable.column(tableID, primID, 2));

                        // Inform createSatfuncTable() of number of active
                        // rows in this table.  Needed to compute slopes of
                        // piecewise linear interpolants.
                        return numActRows;
                    })
                };
            }


            /// Create linearised and padded 'TAB' vector entries of
            /// normalised three-phase SOFN tables for all saturation
            /// function regions from Family One table data.
            ///
            /// \param[in] numRows Number of rows to allocate in the output
            ///    vector for each table.  Expected to be twice the number
            ///    of declared saturation nodes in the simulation run's
            ///    TABDIMS keyword (Item 3).
            ///
            /// \param tolcrit Minimum relative permeability threshold value
            ///    for phase to be considered mobile.  Values less than this
            ///    threshold are output as zero.
            ///
            /// \param[in] slgof Collection of SLGOF tables for all saturation
            ///    regions.
            ///
            /// \param[in] swof Collection of SWOF tables for all saturation
            ///    regions.
            ///
            /// \return Linearised and padded 'TAB' vector values for
            ///    three-phase SOFN tables.  Corresponds to column 1 from
            ///    both of the input SLGOF and SWOF tables, as well as column
            ///    3 from the input SWOF table and column 3 from the input
            ///    SGOF table--expanded so as to have values for all oil
            ///    saturation nodes.  Derivatives added in columns 4 and 5.
            std::tuple<int, int, std::vector<double>>
            fromSLGOFandSWOF(const std::size_t          numRows,
                             const double               tolcrit,
                             const Opm::TableContainer& slgof,
                             const Opm::TableContainer& swof)
            {
                using SLGOF = ::Opm::SlgofTable;
                using SWOF = ::Opm::SwofTable;

                const auto tableSizeGas = detail::allocatedRows(numRows, slgof);
                const auto tableSizeWat = detail::allocatedRows(numRows, swof);

                // Merged table has 2*NSSFUN rows to accommodate both the
                // gas and the water table.
                const auto tableSize =
                    2 * std::max(tableSizeGas, tableSizeWat);

                const auto numTab = slgof.size();
                const auto numDep = std::size_t{2}; // Krow, Krog

                return {
                    static_cast<int>(tableSize),
                    static_cast<int>(numTab),
                    detail::createSatfuncTable(numTab, tableSize, numDep,
                        [tolcrit, &slgof, &swof]
                        (const std::size_t           tableID,
                         const std::size_t           primID,
                         Opm::LinearisedOutputTable& linTable)
                    {
                        const auto sof3 =
                            makeSOF3Table(slgof.getTable<SLGOF>(tableID),
                                          swof.getTable<SWOF>(tableID));

                        auto numActRows = std::size_t{0};

                        // So
                        {
                            const auto& So = sof3[0];

                            numActRows = So.size();
                            std::copy(std::begin(So), std::end(So),
                                      linTable.column(tableID, primID, 0));
                        }

                        // Krow(So)
                        detail::outputRelperm(sof3[1], tolcrit,
                                              linTable.column(tableID, primID, 1));

                        // Krog(So)
                        detail::outputRelperm(sof3[2], tolcrit,
                                              linTable.column(tableID, primID, 2));

                        // Inform createSatfuncTable() of number of active
                        // rows in this table.  Needed to compute slopes of
                        // piecewise linear interpolants.
                        return numActRows;
                    })
                };
            }


            /// Create linearised and padded 'TAB' vector entries of
            /// normalised SOFN tables for all saturation function regions
            /// from Family One table data for a three-phase oil/gas/water
            /// system.
            ///
            /// Distinguishes between the SGOF/SWOF, SLGOF/SWOF and SGOFLET/SWOFLET
            /// family of keywords.
            ///
            /// \param[in] numRows Number of rows to allocate in the output
            ///    vector for each table.  Expected to be equal to the
            ///    number of declared saturation nodes in the simulation
            ///    run's TABDIMS keyword (Item 3).
            ///
            /// \param tolcrit Minimum relative permeability threshold value for
            ///    phase to be considered mobile.  Values less than this threshold
            ///    are output as zero.
            ///
            /// \param[in] tabMgr Collection of run's tabulated functions.
            ///
            /// \return Descriptor and linearised/normalised/padded 'TAB' vector
            ///    values for output SOFN tables.
            std::tuple<int, int, std::vector<double>>
            normalisedTable_One(const std::size_t        numRows,
                                const double             tolcrit,
                                const Opm::TableManager& tabMgr)
            {
                const auto& gasLET = tabMgr.getSgofletTable();
                const auto& watLET = tabMgr.getSwofletTable();

                if (gasLET.empty() || watLET.empty()) {
                    // Common case: Saturation functions for oil entered in
                    // S[L][GW]OF tables.
                    const auto& sgof = tabMgr.getSgofTables();
                    if (!sgof.empty()) {
                        return fromSGOFandSWOF(numRows, tolcrit,
                                            sgof,
                                            tabMgr.getSwofTables());
                    }
                    const auto& slgof = tabMgr.getSlgofTables();
                    if (!slgof.empty()) {
                        return fromSLGOFandSWOF(numRows, tolcrit,
                                            slgof,
                                            tabMgr.getSwofTables());
                    }
                    // Warn here?
                }

                // Saturation functions for oil entered in S[GW]OFLET tables.
                return {
                    static_cast<int>(2 * numRows),
                    static_cast<int>(gasLET.size()),
                    detail::sampleLET(2 * numRows, tolcrit, watLET, gasLET)
                };
            }

            /// Create linearised and padded 'TAB' vector entries of
            /// normalised three-phase SOFN tables for all saturation
            /// function regions from Family Two table data (SOF3 keyword).
            ///
            /// \param[in] numRows Number of rows to allocate in the output
            ///    vector for each table.  Expected to be equal to the
            ///    number of declared saturation nodes in the simulation
            ///    run's TABDIMS keyword (Item 3).
            ///
            /// \param tolcrit Minimum relative permeability threshold value
            ///    for phase to be considered mobile.  Values less than this
            ///    threshold are output as zero.
            ///
            /// \param[in] sof3 Collection of SOF3 tables for all saturation
            ///    regions.
            ///
            /// \return Descriptor and linearised/normalised/padded 'TAB'
            ///    vector values for output SOFN tables.  Essentially a
            ///    padded copy of the input SOF3 tables, \p sof3, with added
            ///    derivatives.
            std::tuple<int, int, std::vector<double>>
            fromSOF3(const std::size_t          numRows,
                     const double               tolcrit,
                     const Opm::TableContainer& sof3)
            {
                using SOF3 = ::Opm::Sof3Table;

                const auto tableSize = detail::allocatedRows(numRows, sof3);

                const auto numTab = sof3.size();
                const auto numDep = std::size_t{2}; // Krow, Krog

                return {
                    static_cast<int>(tableSize),
                    static_cast<int>(numTab),
                    detail::createSatfuncTable(numTab, tableSize, numDep,
                        [tolcrit, &sof3](const std::size_t           tableID,
                                         const std::size_t           primID,
                                         Opm::LinearisedOutputTable& linTable)
                    {
                        const auto& t = sof3.getTable<SOF3>(tableID);

                        auto numActRows = std::size_t{0};

                        // So
                        {
                            const auto& So = t.getSoColumn();

                            numActRows = So.size();
                            std::copy(std::begin(So), std::end(So),
                                      linTable.column(tableID, primID, 0));
                        }

                        // Krow(So)
                        detail::outputRelperm(t.getKrowColumn(), tolcrit,
                                              linTable.column(tableID, primID, 1));

                        // Krog(So)
                        detail::outputRelperm(t.getKrogColumn(), tolcrit,
                                              linTable.column(tableID, primID, 2));

                        // Inform createSatfuncTable() of number of active
                        // rows in this table.  Needed to compute slopes of
                        // piecewise linear interpolants.
                        return numActRows;
                    })
                };
            }

            /// Create linearised and padded 'TAB' vector entries of
            /// normalised SOFN tables for all saturation function regions
            /// from table data for all saturation function families in a
            /// three-phase oil/gas/water system.
            ///
            /// \param[in] family Predicate for current run's saturation
            ///    function family.
            ///
            /// \param[in] numRows Number of rows to allocate in the output
            ///    vector for each table.  Expected to be equal to the
            ///    number of declared saturation nodes in the simulation
            ///    run's TABDIMS keyword (Item 3).
            ///
            /// \param tolcrit Minimum relative permeability threshold value for
            ///    phase to be considered mobile.  Values less than this threshold
            ///    are output as zero.
            ///
            /// \param[in] tabMgr Collection of run's tabulated functions.
            ///
            /// \return Descriptor and linearised/normalised/padded 'TAB' vector
            ///    values for output SOFN tables.
            std::tuple<int, int, std::vector<double>>
            normalisedTable(const Family&            family,
                            const std::size_t        numRows,
                            const double             tolcrit,
                            const Opm::TableManager& tabMgr)
            {
                if (family.two()) {
                    return fromSOF3(numRows, tolcrit, tabMgr.getSof3Tables());
                }

                // Family I: SWOF + SGOF or S{G,W}OFLET
                return normalisedTable_One(numRows, tolcrit, tabMgr);
            }
        } // ThreePhase

        /// Create linearised and padded 'TAB' vector entries of normalised
        /// SOFN tables for all saturation function regions from table data
        /// for all saturation function families.
        ///
        /// \param[in] family Predicate for current run's saturation
        ///    function family.
        ///
        /// \param[in] active Predicate for current run's active phases.
        ///    Needed to query for whether or not this is a three-phase run
        ///    or, if not, which of the other two phases (gas/water) are
        ///    active.
        ///
        /// \param[in] numRows Number of rows to allocate in the output
        ///    vector for each table.  Expected to be equal to the number of
        ///    declared saturation nodes in the simulation run's TABDIMS
        ///    keyword (Item 3).
        ///
        /// \param tolcrit Minimum relative permeability threshold value for
        ///    phase to be considered mobile.  Values less than this threshold
        ///    are output as zero.
        ///
        /// \param[in] tabMgr Collection of run's tabulated functions.
        ///
        /// \return Descriptor and linearised/normalised/padded 'TAB' vector
        ///    values for output SOFN tables.
        std::tuple<int, int, std::vector<double>>
        normalisedTable(const Family&            family,
                        const ActivePhases&      active,
                        const std::size_t        numRows,
                        const double             tolcrit,
                        const Opm::TableManager& tabMgr)
        {
            if (active.isThreePhase()) {
                return ThreePhase::
                    normalisedTable(family, numRows, tolcrit, tabMgr);
            }

            return TwoPhase::normalisedTable(family, active, numRows,
                                             tolcrit, tabMgr);
        }
    } // Oil

    /// Functions to create linearised, padded, and normalised SWFN output
    /// tables from various input saturation function keywords.
    namespace Water {
        /// Create linearised and padded 'TAB' vector entries of normalised
        /// SWFN tables for all saturation function regions from Family Two
        /// table data (SWFN keyword).
        ///
        /// \param[in] numRows Number of rows to allocate for each table in
        ///    the output vector.  Expected to be equal to the number of
        ///    declared saturation nodes in the simulation run's TABDIMS
        ///    keyword (Item 3).
        ///
        /// \param tolcrit Minimum relative permeability threshold value for
        ///    phase to be considered mobile.  Values less than this threshold
        ///    are output as zero.
        ///
        /// \param[in] units Active unit system.  Needed to convert SI
        ///    convention capillary pressure values (Pascal) to declared
        ///    conventions of the run specification.
        ///
        /// \param[in] swfn Collection of SWFN tables for all saturation
        ///    regions.
        ///
        /// \return Descriptor and linearised/normalised/padded 'TAB' vector
        ///    values for output SWFN tables.  A unit-converted copy of the
        ///    input table \p swfn with added derivatives.
        std::tuple<int, int, std::vector<double>>
        fromSWFN(const std::size_t          numRows,
                 const double               tolcrit,
                 const Opm::UnitSystem&     units,
                 const Opm::TableContainer& swfn)
        {
            using SWFN = ::Opm::SwfnTable;

            const auto tableSize = detail::allocatedRows(numRows, swfn);

            const auto numTab = swfn.size();
            const auto numDep = std::size_t{2}; // Krw, Pcow

            return {
                static_cast<int>(tableSize),
                static_cast<int>(numTab),
                detail::createSatfuncTable(numTab, tableSize, numDep,
                    [tolcrit, &swfn, &units]
                        (const std::size_t           tableID,
                         const std::size_t           primID,
                         Opm::LinearisedOutputTable& linTable)
                {
                    const auto& t = swfn.getTable<SWFN>(tableID);

                    auto numActRows = std::size_t{0};

                    // Sw
                    {
                        const auto& Sw = t.getSwColumn();

                        numActRows = Sw.size();
                        std::copy(std::begin(Sw), std::end(Sw),
                                  linTable.column(tableID, primID, 0));
                    }

                    // Krw(Sw)
                    detail::outputRelperm(t.getKrwColumn(), tolcrit,
                                          linTable.column(tableID, primID, 1));

                    // Pcow(Sw)
                    {
                        constexpr auto uPress = ::Opm::UnitSystem::measure::pressure;

                        const auto& pc = t.getPcowColumn();
                        std::transform(std::begin(pc), std::end(pc),
                                       linTable.column(tableID, primID, 2),
                                       [&units](const double Pc) -> double
                                       {
                                           return units.from_si(uPress, Pc);
                                       });
                    }

                    // Inform createSatfuncTable() of number of active rows
                    // in this table.  Needed to compute slopes of piecewise
                    // linear interpolants.
                    return numActRows;
                })
            };
        }

        /// Create linearised and padded 'TAB' vector entries of normalised
        /// SWFN tables for all saturation function regions from Family Two
        /// table data (SGWFN keyword).
        ///
        /// \param[in] numRows Number of rows to allocate for each table in
        ///    the output vector.  Expected to be equal to the number of
        ///    declared saturation nodes in the simulation run's TABDIMS
        ///    keyword (Item 3).
        ///
        /// \param tolcrit Minimum relative permeability threshold value for
        ///    phase to be considered mobile.  Values less than this threshold
        ///    are output as zero.
        ///
        /// \param[in] swof Collection of SWOF tables for all saturation
        ///    regions.
        ///
        /// \return Descriptor and linearised/normalised/padded 'TAB' vector
        ///    values for output SWFN tables.  Corresponds to unit-converted
        ///    copies of columns 1, 3, and 4--with added derivatives--of the
        ///    input SGWFN tables.
        std::tuple<int, int, std::vector<double>>
        fromSGWFN(const std::size_t          numRows,
                  const double               tolcrit,
                  const Opm::TableContainer& sgwfn)
        {
            using SGWFN = ::Opm::SgwfnTable;

            const auto tableSize = detail::allocatedRows(numRows, sgwfn);

            const auto numTab = sgwfn.size();
            const auto numDep = std::size_t{2}; // Krw, Pcow

            return {
                static_cast<int>(tableSize),
                static_cast<int>(numTab),
                detail::createSatfuncTable(numTab, tableSize, numDep,
                    [tolcrit, &sgwfn]
                        (const std::size_t           tableID,
                         const std::size_t           primID,
                         Opm::LinearisedOutputTable& linTable)
                {
                    const auto& t = sgwfn.getTable<SGWFN>(tableID);

                    auto numActRows = std::size_t{0};

                    // Sw
                    {
                        const auto& Sg = t.getSgColumn();
                        numActRows = Sg.size();

                        auto Sw = std::vector<double>{};
                        Sw.reserve(numActRows);

                        std::transform(std::begin(Sg), std::end(Sg),
                                       std::back_inserter(Sw),
                                       [](const auto sg) { return 1.0 - sg; });

                        std::copy(Sw.rbegin(), Sw.rend(),
                                  linTable.column(tableID, primID, 0));
                    }

                    // Krw(Sw)
                    {
                        const auto& kr = t.getKrgwColumn();

                        const auto krw = std::vector<double> {
                            std::begin(kr), std::end(kr)
                        };

                        detail::outputRelperm(krw.rbegin(), krw.rend(), tolcrit,
                                              linTable.column(tableID, primID, 1));
                    }

                    // Pcow(Sw) = zero
                    {
                        const auto& pc = t.getPcgwColumn();

                        std::transform(std::begin(pc), std::end(pc),
                                       linTable.column(tableID, primID, 2),
                                       [](const double) -> double
                                       {
                                           return 0.0;
                                       });
                    }

                    // Inform createSatfuncTable() of number of active rows
                    // in this table.  Needed to compute slopes of piecewise
                    // linear interpolants.
                    return numActRows;
                })
            };
        }

        /// Create linearised and padded 'TAB' vector entries of normalised
        /// SGFN tables for all saturation function regions from Family Two
        /// table data.
        ///
        /// Distinguishes between keywords SWFN and SGWFN.
        ///
        /// \param[in] numRows Number of rows to allocate in the output
        ///    vector for each table.  Expected to be equal to the number of
        ///    declared saturation nodes in the simulation run's TABDIMS
        ///    keyword (Item 3).
        ///
        /// \param tolcrit Minimum relative permeability threshold value for
        ///    phase to be considered mobile.  Values less than this threshold
        ///    are output as zero.
        ///
        /// \param[in] units Active unit system.  Needed to convert SI
        ///    convention capillary pressure values (Pascal) to declared
        ///    conventions of the run specification.
        ///
        /// \param[in] tabMgr Collection of run's tabulated functions.
        ///
        /// \return Descriptor and linearised/normalised/padded 'TAB' vector
        ///    values for output SWFN tables.
        std::tuple<int, int, std::vector<double>>
        normalisedTable_Two(const std::size_t        numRows,
                            const double             tolcrit,
                            const Opm::UnitSystem&   units,
                            const Opm::TableManager& tabMgr)
        {
            if (const auto& swfn = tabMgr.getSwfnTables(); ! swfn.empty()) {
                // Typical case: Saturation functions for water in SWFN table
                return fromSWFN(numRows, tolcrit, units, swfn);
            }

            // Saturation functions for water given in SGWFN tables
            return fromSGWFN(numRows, tolcrit, tabMgr.getSgwfnTables());
        }

        /// Create linearised and padded 'TAB' vector entries of normalised
        /// SWFN tables for all saturation function regions from Family One
        /// table data (SWOF keyword).
        ///
        /// \param[in] numRows Number of rows to allocate for each table in
        ///    the output vector.  Expected to be equal to the number of
        ///    declared saturation nodes in the simulation run's TABDIMS
        ///    keyword (Item 3).
        ///
        /// \param tolcrit Minimum relative permeability threshold value for
        ///    phase to be considered mobile.  Values less than this threshold
        ///    are output as zero.
        ///
        /// \param[in] units Active unit system.  Needed to convert SI
        ///    convention capillary pressure values (Pascal) to declared
        ///    conventions of the run specification.
        ///
        /// \param[in] swof Collection of SWOF tables for all saturation
        ///    regions.
        ///
        /// \return Descriptor and linearised/normalised/padded 'TAB' vector
        ///    values for output SWFN tables.  Corresponds to unit-converted
        ///    copies of columns 1, 2, and 4--with added derivatives--of the
        ///    input SWOF tables.
        std::tuple<int, int, std::vector<double>>
        fromSWOF(const std::size_t          numRows,
                 const double               tolcrit,
                 const Opm::UnitSystem&     units,
                 const Opm::TableContainer& swof)
        {
            using SWOF = ::Opm::SwofTable;

            const auto tableSize = detail::allocatedRows(numRows, swof);

            const auto numTab = swof.size();
            const auto numDep = std::size_t{2}; // Krw, Pcow

            return {
                static_cast<int>(tableSize),
                static_cast<int>(numTab),
                detail::createSatfuncTable(numTab, tableSize, numDep,
                    [tolcrit, &swof, &units]
                        (const std::size_t           tableID,
                         const std::size_t           primID,
                         Opm::LinearisedOutputTable& linTable)
                {
                    const auto& t = swof.getTable<SWOF>(tableID);

                    auto numActRows = std::size_t{0};

                    // Sw
                    {
                        const auto& Sw = t.getSwColumn();

                        numActRows = Sw.size();
                        std::copy(std::begin(Sw), std::end(Sw),
                                  linTable.column(tableID, primID, 0));
                    }

                    // Krw(Sw)
                    detail::outputRelperm(t.getKrwColumn(), tolcrit,
                                          linTable.column(tableID, primID, 1));

                    // Pcow(Sw)
                    {
                        constexpr auto uPress = ::Opm::UnitSystem::measure::pressure;

                        const auto& pc = t.getPcowColumn();
                        std::transform(std::begin(pc), std::end(pc),
                                       linTable.column(tableID, primID, 2),
                                       [&units](const double Pc) -> double
                                       {
                                           return units.from_si(uPress, Pc);
                                       });
                    }

                    // Inform createSatfuncTable() of number of active rows
                    // in this table.  Needed to compute slopes of piecewise
                    // linear interpolants.
                    return numActRows;
                })
            };
        }

        /// Create linearised and padded 'TAB' vector entries of normalised
        /// SGFN tables for all saturation function regions from Family One
        /// table data.
        ///
        /// Distinguishes between keywords SWOF and SWOFLET.
        ///
        /// \param[in] numRows Number of rows to allocate in the output
        ///    vector for each table.  Expected to be equal to the number of
        ///    declared saturation nodes in the simulation run's TABDIMS
        ///    keyword (Item 3).
        ///
        /// \param tolcrit Minimum relative permeability threshold value for
        ///    phase to be considered mobile.  Values less than this threshold
        ///    are output as zero.
        ///
        /// \param[in] units Active unit system.  Needed to convert SI
        ///    convention capillary pressure values (Pascal) to declared
        ///    conventions of the run specification.
        ///
        /// \param[in] tabMgr Collection of run's tabulated functions.
        ///
        /// \return Descriptor and linearised/normalised/padded 'TAB' vector
        ///    values for output SWFN tables.
        std::tuple<int, int, std::vector<double>>
        normalisedTable_One(const std::size_t        numRows,
                            const double             tolcrit,
                            const Opm::UnitSystem&   units,
                            const Opm::TableManager& tabMgr)
        {
            const auto& letTables = tabMgr.getSwofletTable();

            if (letTables.empty()) {
                // Typical case: Saturation functions for water in SWOF
                return fromSWOF(numRows, tolcrit, units, tabMgr.getSwofTables());
            }

            return { // LET-based saturation functions for water
                static_cast<int>(numRows),
                static_cast<int>(letTables.size()),
                detail::sampleLET(numRows, tolcrit, units, letTables,
                                  std::vector<double>(letTables.size(), 0.0))
            };
        }

        /// Create linearised and padded 'TAB' vector entries of normalised
        /// WSF tables for all saturation function regions from Family Three
        /// table data (WSF keyword).
        ///
        /// \param[in] numRows Number of rows to allocate for each table in
        ///    the output vector.  Expected to be equal to the number of
        ///    declared saturation nodes in the simulation run's TABDIMS
        ///    keyword (Item 3).
        ///
        /// \param tolcrit Minimum relative permeability threshold value for
        ///    phase to be considered mobile.  Values less than this threshold
        ///    are output as zero.
        ///
        /// \param[in] units Active unit system.  Needed to convert SI
        ///    convention capillary pressure values (Pascal) to declared
        ///    conventions of the run specification.
        ///
        /// \param[in] wsf Collection of WSF tables for all saturation
        ///    regions.
        ///
        /// \return Descriptor and linearised/normalised/padded 'TAB' vector
        ///    values for output SWFN tables.  A unit-converted copy of the
        ///    input table \p wsf with added derivatives.
        std::tuple<int, int, std::vector<double>>
        fromWSF(const std::size_t          numRows,
                const double               tolcrit,
                const Opm::TableContainer& wsf)
        {
            using WSF = ::Opm::WsfTable;

            const auto tableSize = detail::allocatedRows(numRows, wsf);

            const auto numTab = wsf.size();
            const auto numDep = std::size_t{2}; // Krw, {zero pc}

            return {
                static_cast<int>(tableSize),
                static_cast<int>(numTab),
                detail::createSatfuncTable(numTab, tableSize, numDep,
                    [tolcrit, &wsf]
                        (const std::size_t           tableID,
                         const std::size_t           primID,
                         Opm::LinearisedOutputTable& linTable)
                {
                    const auto& t = wsf.getTable<WSF>(tableID);

                    auto numActRows = std::size_t{0};

                    // Sw
                    {
                        const auto& Sw = t.getSwColumn();

                        numActRows = Sw.size();
                        std::copy(std::begin(Sw), std::end(Sw),
                                  linTable.column(tableID, primID, 0));
                    }

                    // Krw(Sw)
                    detail::outputRelperm(t.getKrwColumn(), tolcrit,
                                          linTable.column(tableID, primID, 1));

                    // Pcow = zero
                    {
                        const auto& krw = t.getKrwColumn();
                        std::transform(std::begin(krw), std::end(krw),
                                       linTable.column(tableID, primID, 2),
                                       [](const double) { return 0.0; });
                    }

                    // Inform createSatfuncTable() of number of active rows
                    // in this table.  Needed to compute slopes of piecewise
                    // linear interpolants.
                    return numActRows;
                })
            };
        }

        /// Create linearised and padded 'TAB' vector entries of normalised
        /// SWFN tables for all saturation function regions from table data
        /// for all saturation function families.
        ///
        /// \param[in] family Predicate for current run's saturation
        ///    function family.
        ///
        /// \param[in] numRows Number of rows to allocate in the output
        ///    vector for each table.  Expected to be equal to the number of
        ///    declared saturation nodes in the simulation run's TABDIMS
        ///    keyword (Item 3).
        ///
        /// \param tolcrit Minimum relative permeability threshold value for
        ///    phase to be considered mobile.  Values less than this threshold
        ///    are output as zero.
        ///
        /// \param[in] units Active unit system.  Needed to convert SI
        ///    convention capillary pressure values (Pascal) to declared
        ///    conventions of the run specification.
        ///
        /// \param[in] tabMgr Collection of run's tabulated functions.
        ///
        /// \return Descriptor and linearised/normalised/padded 'TAB' vector
        ///    values for output SWFN tables.
        std::tuple<int, int, std::vector<double>>
        normalisedTable(const Family&            family,
                        const std::size_t        numRows,
                        const double             tolcrit,
                        const Opm::UnitSystem&   units,
                        const Opm::TableManager& tabMgr)
        {
            if (family.one()) {
                return normalisedTable_One(numRows, tolcrit, units, tabMgr);
            }

            if (family.two()) {
                return normalisedTable_Two(numRows, tolcrit, units, tabMgr);
            }

            return fromWSF(numRows, tolcrit, tabMgr.getWsfTables());
        }
    } // Water

}} // Anonymous::SatFunc

/// Functions to facilitate generating TAB vector entries for tabulated
/// PVT functions.
namespace { namespace PVTFunc {

    /// Functions to create linearised, padded, and normalised gas PVT
    /// output tables from various input gas PVT function keywords.
    namespace Gas {
        /// Create linearised and padded 'TAB' vector entries of normalised
        /// gas tables for all PVT function regions from PVDG (dry gas)
        /// keyword data.
        ///
        /// \param[in] numPressNodes Number of pressure nodes (rows) to
        ///    allocate in the output vector for each table.  Expected to be
        ///    equal to the number of declared pressure nodes in the
        ///    simulation run's TABDIMS keyword (NPPVT, Item 4).
        ///
        /// \param[in] units Active unit system.  Needed to convert SI
        ///    convention pressure values (Pascal), formation volume factors
        ///    (Rm3/Sm3), and viscosity values (Pascal seconds) to declared
        ///    conventions of the run specification.
        ///
        /// \param[in] pvdg Collection of PVDG tables for all PVT regions.
        ///
        /// \return Linearised and padded 'TAB' vector values for output gas
        ///    PVT tables.  A unit-converted copy of the input table \p pvdg
        ///    with added derivatives.
        std::vector<double>
        fromPVDG(const std::size_t          numPressNodes,
                 const Opm::UnitSystem&     units,
                 const Opm::TableContainer& pvdg)
        {
            // Columns [ Pg, 1/Bg, 1/(Bg*mu_g), derivatives ]
            using PVDG = ::Opm::PvdgTable;

            const auto numTab  = pvdg.size();
            const auto numPrim = std::size_t{1}; // No sub-tables
            const auto numRows = numPressNodes;  // One row per pressure node
            const auto numDep  = std::size_t{2}; // 1/Bg, 1/(Bg*mu_g)

            // PVDG fill value = +2.0e20
            const auto fillVal = +2.0e20;

            return createPropfuncTable(numTab, numPrim, numRows, numDep, fillVal,
                [&units, &pvdg](const std::size_t           tableID,
                                const std::size_t           primID,
                                Opm::LinearisedOutputTable& linTable)
                -> std::size_t
            {
                const auto& t = pvdg.getTable<PVDG>(tableID);

                auto numActRows = std::size_t{0};

                // Column 0: Pg
                {
                    constexpr auto uPress = ::Opm::UnitSystem::measure::pressure;

                    const auto& Pg = t.getPressureColumn();

                    numActRows = Pg.size();

                    std::transform(std::begin(Pg), std::end(Pg),
                                   linTable.column(tableID, primID, 0),
                        [&units](const double p) -> double
                    {
                        return units.from_si(uPress, p);
                    });
                }

                // Column 1: 1/Bg
                {
                    constexpr auto uRecipFVF = ::Opm::UnitSystem::measure::
                        gas_inverse_formation_volume_factor;

                    const auto& Bg = t.getFormationFactorColumn();
                    std::transform(std::begin(Bg), std::end(Bg),
                                   linTable.column(tableID, primID, 1),
                        [&units](const double B) -> double
                    {
                        return units.from_si(uRecipFVF, 1.0 / B);
                    });
                }

                // Column 2: 1/(Bg*mu_g)
                {
                    constexpr auto uRecipFVF = ::Opm::UnitSystem::measure::
                        gas_inverse_formation_volume_factor;

                    constexpr auto uVisc = ::Opm::UnitSystem::measure::viscosity;

                    const auto& Bg   = t.getFormationFactorColumn();
                    const auto& mu_g = t.getViscosityColumn();

                    std::transform(std::begin(Bg), std::end(Bg),
                                   std::begin(mu_g),
                                   linTable.column(tableID, primID, 2),
                        [&units](const double B, const double mu) -> double
                    {
                        return units.from_si(uRecipFVF, 1.0 / B)
                            /  units.from_si(uVisc    , mu);
                    });
                }

                // Inform createPropfuncTable() of number of active rows in
                // this table.  Needed to compute slopes of piecewise linear
                // interpolants.
                return numActRows;
            });
        }

        /// Create linearised and padded 'TAB' vector entries of normalised
        /// gas tables for all PVT function regions from PVTG (wet gas with
        /// volatile/vaporised oil) keyword data.
        ///
        /// \param[in] numCompNodes Number of composition nodes (rows per
        ///    sub-table) to allocate in the output vector for each table.
        ///    Expected to be equal to the number of declared composition
        ///    nodes in the simulation run's TABDIMS keyword (NRPVT, Item
        ///    6).
        ///
        /// \param[in] numPressNodes Number of pressure nodes (sub-tables)
        ///    to allocate in the output vector for each gas PVT table.
        ///    Expected to be equal to the number of declared pressure nodes
        ///    in the simulation run's TABDIMS keyword (NPPVT, Item 4).
        ///
        /// \param[in] units Active unit system.  Needed to convert SI
        ///    convention vaporised oil composition values (Sm3/Sm3),
        ///    formation volume factors (Rm3/Sm3), and viscosity values
        ///    (Pa*s) to declared conventions of the run specification.
        ///
        /// \param[in] pvtg Collection of PVTG tables for all PVT regions.
        ///
        /// \return Linearised and padded 'TAB' vector values for output gas
        ///    PVT tables.  A unit-converted copy of the input table \p pvtg
        ///    with added derivatives.
        std::vector<double>
        fromPVTG(const std::size_t                  numCompNodes,
                 const std::size_t                  numPressNodes,
                 const Opm::UnitSystem&             units,
                 const std::vector<Opm::PvtgTable>& pvtg)
        {
            // Columns [ Rv, 1/Bg, 1/(Bg*mu_g), derivatives ]
            const auto numTab  = pvtg.size();
            const auto numPrim = numPressNodes;
            const auto numRows = numCompNodes;
            const auto numDep  = std::size_t{2}; // 1/Bg, 1/(Bg*mu_g)

            // PVTG fill value = -2.0e20
            const auto fillVal = -2.0e20;

            return createPropfuncTable(numTab, numPrim, numRows, numDep, fillVal,
                [&units, &pvtg](const std::size_t           tableID,
                                const std::size_t           primID,
                                Opm::LinearisedOutputTable& linTable)
                -> std::size_t
            {
                auto numActRows = std::size_t{0};

                if (primID >= pvtg[tableID].size()) {
                    // Composition node outside current table's active set.
                    // No active rows in this subtable.
                    return numActRows;
                }

                const auto& t = pvtg[tableID].getUnderSaturatedTable(primID);

                // Column 0: Rv
                {
                    constexpr auto uRv = ::Opm::UnitSystem::measure::oil_gas_ratio;

                    const auto& Rv = t.getColumn(0);

                    numActRows = Rv.size();

                    std::transform(std::begin(Rv), std::end(Rv),
                                   linTable.column(tableID, primID, 0),
                        [&units](const double rv) -> double
                    {
                        return units.from_si(uRv, rv);
                    });
                }

                // Column 1: 1/Bg
                {
                    constexpr auto uRecipFVF = ::Opm::UnitSystem::measure::
                        gas_inverse_formation_volume_factor;

                    const auto& Bg = t.getColumn(1);

                    std::transform(std::begin(Bg), std::end(Bg),
                                   linTable.column(tableID, primID, 1),
                        [&units](const double B) -> double
                    {
                        return units.from_si(uRecipFVF, 1.0 / B);
                    });
                }

                // Column 2: 1/(Bg*mu_g)
                {
                    constexpr auto uRecipFVF = ::Opm::UnitSystem::measure::
                        gas_inverse_formation_volume_factor;

                    constexpr auto uVisc = ::Opm::UnitSystem::measure::viscosity;

                    const auto& Bg   = t.getColumn(1);
                    const auto& mu_g = t.getColumn(2);

                    std::transform(std::begin(Bg), std::end(Bg),
                                   std::begin(mu_g),
                                   linTable.column(tableID, primID, 2),
                        [&units](const double B, const double mu) -> double
                    {
                        return units.from_si(uRecipFVF, 1.0 / B)
                            /  units.from_si(uVisc    , mu);
                    });
                }

                // Inform createPropfuncTable() of number of active rows in
                // this table.  Needed to compute slopes of piecewise linear
                // interpolants.
                return numActRows;
            });
        }

        /// Create linearised and padded 'TAB' vector entries of normalised
        /// gas pressure nodes for all PVT function regions from PVTG (wet
        /// gas with volatile/vaporised oil) keyword data.
        ///
        /// \param[in] numPressNodes Number of pressure nodes to allocate in
        ///    the output vector for each gas PVT table.  Expected to be
        ///    equal to the number of declared pressure nodes in the
        ///    simulation run's TABDIMS keyword (NPPVT, Item 4).
        ///
        /// \param[in] units Active unit system.  Needed to convert SI
        ///    convention pressure values (Pascal) to declared conventions
        ///    of the run specification.
        ///
        /// \param[in] pvtg Collection of PVTG tables for all PVT regions.
        ///
        /// \return Linearised and padded 'TAB' vector values for output gas
        ///    PVT tables.  A unit-converted copy of the primary keys in
        ///    input table \p pvtg.
        std::vector<double>
        pressureNodes(const std::size_t                  numPressNodes,
                      const Opm::UnitSystem&             units,
                      const std::vector<Opm::PvtgTable>& pvtg)
        {
            // Columns [ Pg ]
            const auto numTab = pvtg.size();

            // One set of pressure nodes per table.
            const auto numPrim = std::size_t{1};

            const auto numRows = numPressNodes;

            // No dependent variables.
            const auto numDep = std::size_t{0};

            // Pressure node fill value = +2.0e20
            const auto fillVal = +2.0e20;

            return createPropfuncTable(numTab, numPrim, numRows, numDep, fillVal,
                [&units, &pvtg](const std::size_t           tableID,
                                const std::size_t           primID,
                                Opm::LinearisedOutputTable& linTable)
                -> std::size_t
            {
                constexpr auto uPress = ::Opm::UnitSystem::measure::pressure;

                const auto& t  = pvtg[tableID].getSaturatedTable();
                const auto& Pg = t.getColumn(0);

                const auto numActRows = Pg.size();

                std::transform(std::begin(Pg), std::end(Pg),
                               linTable.column(tableID, primID, 0),
                    [&units](const double p) -> double
                {
                    return units.from_si(uPress, p);
                });

                return numActRows;
            });
        }

        /// Extract maximum effective composition nodes in PVTG table data
        ///
        /// \param[in] pvtg Collection of PVTG tables for all PVT regions.
        ///
        /// \return Maximum number of active rows across all sub-tables of \p
        ///    pvtg.
        std::size_t maxNumCompNodes(const std::vector<Opm::PvtgTable>& pvtg)
        {
            auto max_nrpvt = std::size_t{0};

            for (const auto& table : pvtg) {
                const auto max_nrpvt_tab =
                    std::accumulate(table.begin(), table.end(), std::size_t{0},
                        [](const std::size_t n, const Opm::SimpleTable& t)
                    {
                        return std::max(n, t.numRows());
                    });

                max_nrpvt = std::max(max_nrpvt, max_nrpvt_tab);
            }

            return max_nrpvt;
        }

        /// Extract maximum effective pressure nodes in PVDG table data
        ///
        /// \param[in] pvdg Collection of PVDG tables for all PVT regions.
        ///
        /// \return Maximum number of table rows across all tables of \p
        ///    pvdo.
        std::size_t maxNumPressNodes(const Opm::TableContainer& pvdg)
        {
            using PVDG = ::Opm::PvdgTable;
            auto tabID = std::vector<std::size_t>(pvdg.size());

            std::iota(tabID.begin(), tabID.end(), std::size_t{0});

            return std::accumulate(tabID.begin(), tabID.end(), std::size_t{0},
                [&pvdg](const std::size_t n, const std::size_t table)
            {
                return std::max(n, pvdg.getTable<PVDG>(table).numRows());
            });
        }

        /// Extract maximum effective pressure nodes in PVTG table data
        ///
        /// \param[in] pvtg Collection of PVTG tables for all PVT regions.
        ///
        /// \return Maximum number of active keys across all tables of \p
        ///    pvtg.
        std::size_t maxNumPressNodes(const std::vector<Opm::PvtgTable>& pvtg)
        {
            if (pvtg.empty()) { return 0; }

            return std::accumulate(std::begin(pvtg), std::end(pvtg), std::size_t{0},
                [](const std::size_t n, const Opm::PvtgTable& t) -> std::size_t
            {
                return std::max(n, t.getSaturatedTable().numRows());
            });
        }
    } // Gas

    /// Functions to create linearised, padded, and normalised oil PVT
    /// output tables from various input oil PVT function keywords.
    namespace Oil {
        /// Create linearised and padded 'TAB' vector entries of normalised
        /// oil tables for all PVT function regions from PVCDO (dead oil
        /// with constant oil compressibility) keyword data.
        ///
        /// \param[in] numPressNodes Number of pressure nodes (rows) to
        ///    allocate in the output vector for each table.  Expected to be
        ///    equal to the number of declared pressure nodes in the
        ///    simulation run's TABDIMS keyword (NPPVT, Item 4).  Note that
        ///    only the first row of each table contains any actual PVT
        ///    data.
        ///
        /// \param[in] units Active unit system.  Needed to convert SI
        ///    convention pressure values (Pascal), formation volume factors
        ///    (Rm3/Sm3), compressibility values (1/Pascal), viscosity
        ///    values (Pa*s), and viscosibility values (1/Pascal) to
        ///    declared conventions of the run specification.
        ///
        /// \param[in] pvcdo Collection of PVCDO tables for all PVT regions.
        ///
        /// \return Linearised and padded 'TAB' vector values for output oil
        ///    PVT tables.  A unit-converted copy of the input table \p
        ///    pvcdo.  No derivative information added.
        std::vector<double>
        fromPVCDO(const std::size_t      numPressNodes,
                  const Opm::UnitSystem& units,
                  const Opm::PvcdoTable& pvcdo)
        {
            // Recall: PvcdoTable is essentially vector<PVCDORecord> with
            //
            //   struct PVCDORecord {
            //       double reference_pressure;
            //       double volume_factor;
            //       double compressibility;
            //       double viscosity;
            //       double viscosibility;
            //   };
            //
            using M = ::Opm::UnitSystem::measure;

            // Columns [ Po, Bo, Co, mu_o, Cv ]
            //
            // Single active row per table.  No derivatives.  Can't reuse
            // createPropfuncTable here, so implement the return value
            // directly in terms of LinearisedOutputTable.

            const auto numTab  = pvcdo.size();
            const auto numPrim = std::size_t{1};
            const auto numCols = std::size_t{5};

            // PVCDO fill value: -1.0e20
            const auto fillVal = -1.0e20;

            auto lintable = ::Opm::LinearisedOutputTable {
                numTab, numPrim, numPressNodes, numCols, fillVal
            };

            // Note unit hack for compressibility and viscosibility.  The
            // unit of measurement for these quantities is 1/pressure, but
            // the UnitSystem does not define this unit.  Work around the
            // missing conversion by using *to_si()* rather than *from_si()*
            // for those quantities.

            const auto uPress = M::pressure;
            const auto uBo    = M::oil_formation_volume_factor;
            const auto uVisc  = M::viscosity;

            // Single primary key, ID = 0.
            const auto primID = std::size_t{0};

            for (auto tabID = 0*numTab; tabID < numTab; ++tabID) {
                const auto& t = pvcdo[tabID];

                auto iPo   = lintable.column(tabID, primID, 0);
                auto iBo   = lintable.column(tabID, primID, 1);
                auto iCo   = lintable.column(tabID, primID, 2);
                auto imu_o = lintable.column(tabID, primID, 3);
                auto iCv   = lintable.column(tabID, primID, 4);

                *iPo = units.from_si(uPress, t.reference_pressure);
                *iBo = units.from_si(uBo, t.volume_factor);

                // Compressibility unit hack here (*to_si()*)
                *iCo = units.to_si(uPress, t.compressibility);

                *imu_o = units.from_si(uVisc, t.viscosity);

                // Viscosibility unit hack here (*to_si()*)
                *iCv = units.to_si(uPress, t.viscosibility);
            }

            return lintable.getDataDestructively();
        }

        /// Create linearised and padded 'TAB' vector entries of normalised
        /// gas tables for all PVT function regions from PVDO (dead oil)
        /// keyword data.
        ///
        /// \param[in] numPressNodes Number of pressure nodes (rows) to
        ///    allocate in the output vector for each table.  Expected to be
        ///    equal to the number of declared pressure nodes in the
        ///    simulation run's TABDIMS keyword (NPPVT, Item 4).
        ///
        /// \param[in] units Active unit system.  Needed to convert SI
        ///    convention pressure values (Pascal), formation volume factors
        ///    (Rm3/Sm3), and viscosity values (Pa*s) to declared
        ///    conventions of the run specification.
        ///
        /// \param[in] pvdo Collection of PVDO tables for all PVT regions.
        ///
        /// \return Linearised and padded 'TAB' vector values for output gas
        ///    PVT tables.  A unit-converted copy of the input table \p pvdo
        ///    with added derivatives.
        std::vector<double>
        fromPVDO(const std::size_t          numPressNodes,
                 const Opm::UnitSystem&     units,
                 const Opm::TableContainer& pvdo)
        {
            // Columns [ Po, 1/Bo, 1/(Bo*mu_o), derivatives ]
            using PVDO = ::Opm::PvdoTable;

            const auto numTab  = pvdo.size();
            const auto numPrim = std::size_t{1}; // No sub-tables
            const auto numRows = numPressNodes;  // One row per pressure node.
            const auto numDep  = std::size_t{2}; // 1/Bo, 1/(Bo*mu_o)

            // PVDO fill value = +2.0e20
            const auto fillVal = +2.0e20;

            return createPropfuncTable(numTab, numPrim, numRows, numDep, fillVal,
                [&units, &pvdo](const std::size_t           tableID,
                                const std::size_t           primID,
                                Opm::LinearisedOutputTable& linTable)
                -> std::size_t
            {
                const auto& t = pvdo.getTable<PVDO>(tableID);

                auto numActRows = std::size_t{0};

                // Column 0: Po
                {
                    constexpr auto uPress = ::Opm::UnitSystem::measure::pressure;

                    const auto& Po = t.getPressureColumn();

                    numActRows = Po.size();

                    std::transform(std::begin(Po), std::end(Po),
                                   linTable.column(tableID, primID, 0),
                        [&units](const double p) -> double
                    {
                        return units.from_si(uPress, p);
                    });
                }

                // Column 1: 1/Bo
                {
                    constexpr auto uRecipFVF = ::Opm::UnitSystem::measure::
                        oil_inverse_formation_volume_factor;

                    const auto& Bo = t.getFormationFactorColumn();
                    std::transform(std::begin(Bo), std::end(Bo),
                                   linTable.column(tableID, primID, 1),
                        [&units](const double B) -> double
                    {
                        return units.from_si(uRecipFVF, 1.0 / B);
                    });
                }

                // Column 2: 1/(Bo*mu_o)
                {
                    constexpr auto uRecipFVF = ::Opm::UnitSystem::measure::
                        oil_inverse_formation_volume_factor;

                    constexpr auto uVisc = ::Opm::UnitSystem::measure::viscosity;

                    const auto& Bo   = t.getFormationFactorColumn();
                    const auto& mu_o = t.getViscosityColumn();

                    std::transform(std::begin(Bo), std::end(Bo),
                                   std::begin(mu_o),
                                   linTable.column(tableID, primID, 2),
                        [&units](const double B, const double mu) -> double
                    {
                        return units.from_si(uRecipFVF, 1.0 / B)
                            /  units.from_si(uVisc    , mu);
                    });
                }

                // Inform createPropfuncTable() of number of active rows in
                // this table.  Needed to compute slopes of piecewise linear
                // interpolants.
                return numActRows;
            });
        }

        /// Create linearised and padded 'TAB' vector entries of normalised
        /// gas tables for all PVT function regions from PVTO (live oil with
        /// dissolved gas) keyword data.
        ///
        /// \param[in] numCompNodes Number of composition nodes (sub-tables)
        ///    to allocate in the output vector for each oil PVT table.
        ///    Expected to be equal to the number of declared composition
        ///    nodes in the simulation run's TABDIMS keyword (NRPVT, Item
        ///    6).
        ///
        /// \param[in] numPressNodes Number of pressure nodes (rows per
        ///    sub-table) to allocate in the output vector for each oil PVT
        ///    table.  Expected to be equal to the number of declared
        ///    pressure nodes in the simulation run's TABDIMS keyword
        ///    (NPPVT, Item 4).
        ///
        /// \param[in] units Active unit system.  Needed to convert SI
        ///    convention pressure values (Pascal), formation volume factors
        ///    (S-volume to R-volume), and viscosity values (Pa*s) to
        ///    declared conventions of the run specification.
        ///
        /// \param[in] pvto Collection of PVTO tables for all PVT regions.
        ///
        /// \return Linearised and padded 'TAB' vector values for output gas
        ///    PVT tables.  A unit-converted copy of the input table \p pvto
        ///    with added derivatives.
        std::vector<double>
        fromPVTO(const std::size_t                  numCompNodes,
                 const std::size_t                  numPressNodes,
                 const Opm::UnitSystem&             units,
                 const std::vector<Opm::PvtoTable>& pvto)
        {
            // Columns [ Po, 1/Bo, 1/(Bo*mu_o), derivatives ]
            const auto numTab  = pvto.size();
            const auto numPrim = numCompNodes;
            const auto numRows = numPressNodes;
            const auto numDep  = std::size_t{2}; // 1/Bo, 1/(Bo*mu_o)

            // PVTO fill value = +2.0e20
            const auto fillVal = +2.0e20;

            return createPropfuncTable(numTab, numPrim, numRows, numDep, fillVal,
                [&units, &pvto](const std::size_t           tableID,
                                const std::size_t           primID,
                                Opm::LinearisedOutputTable& linTable)
                -> std::size_t
            {
                auto numActRows = std::size_t{0};

                if (primID >= pvto[tableID].size()) {
                    // Composition node outside current table's active set.
                    // No active rows in this sub-table.
                    return numActRows;
                }

                const auto& t = pvto[tableID].getUnderSaturatedTable(primID);

                // Column 0: Po
                {
                    constexpr auto uPo = ::Opm::UnitSystem::measure::pressure;

                    const auto& Po = t.getColumn(0);

                    numActRows = Po.size();

                    std::transform(std::begin(Po), std::end(Po),
                                   linTable.column(tableID, primID, 0),
                        [&units](const double po) -> double
                    {
                        return units.from_si(uPo, po);
                    });
                }

                // Column 1: 1/Bo
                {
                    constexpr auto uRecipFVF = ::Opm::UnitSystem::measure::
                        oil_inverse_formation_volume_factor;

                    const auto& Bo = t.getColumn(1);

                    std::transform(std::begin(Bo), std::end(Bo),
                                   linTable.column(tableID, primID, 1),
                        [&units](const double B) -> double
                    {
                        return units.from_si(uRecipFVF, 1.0 / B);
                    });
                }

                // Column 2: 1/(Bo*mu_o)
                {
                    constexpr auto uRecipFVF = ::Opm::UnitSystem::measure::
                        oil_inverse_formation_volume_factor;

                    constexpr auto uVisc = ::Opm::UnitSystem::measure::viscosity;

                    const auto& Bo   = t.getColumn(1);
                    const auto& mu_o = t.getColumn(2);

                    std::transform(std::begin(Bo), std::end(Bo),
                                   std::begin(mu_o),
                                   linTable.column(tableID, primID, 2),
                        [&units](const double B, const double mu) -> double
                    {
                        return units.from_si(uRecipFVF, 1.0 / B)
                            /  units.from_si(uVisc    , mu);
                    });
                }

                // Inform createPropfuncTable() of number of active rows in
                // this table.  Needed to compute slopes of piecewise linear
                // interpolants.
                return numActRows;
            });
        }

        /// Create linearised and padded 'TAB' vector entries of normalised
        /// composition nodes for all PVT function regions from PVTO (live
        /// oil with dissolved gas) keyword data.
        ///
        /// \param[in] numCompNodes Number of composition nodes to allocate
        ///    in the output vector for each oil PVT table.  Expected to be
        ///    equal to the number of declared composition nodes in the
        ///    simulation run's TABDIMS keyword (NRPVT, Item 6).
        ///
        /// \param[in] units Active unit system.  Needed to convert SI
        ///    convention dissolved gas composition values (Sm3/Sm3) to
        ///    declared conventions of the run specification.
        ///
        /// \param[in] pvto Collection of PVTO tables for all PVT regions.
        ///
        /// \return Linearised and padded 'TAB' vector values for output oil
        ///    PVT tables.  A unit-converted copy of the primary keys in
        ///    input table \p pvto.
        std::vector<double>
        compositionNodes(const std::size_t                  numCompNodes,
                         const Opm::UnitSystem&             units,
                         const std::vector<Opm::PvtoTable>& pvto)
        {
            // Columns [ Rs ]
            const auto numTab = pvto.size();

            // One set of composition nodes per table.
            const auto numPrim = std::size_t{1};
            const auto numRows = numCompNodes;

            // No dependent variables.
            const auto numDep = std::size_t{0};

            // PVTO fill value = +2.0e20
            const auto fillVal = +2.0e20;

            return createPropfuncTable(numTab, numPrim, numRows, numDep, fillVal,
                [&units, &pvto](const std::size_t           tableID,
                                const std::size_t           primID,
                                Opm::LinearisedOutputTable& linTable)
                -> std::size_t
            {
                constexpr auto uRs = ::Opm::UnitSystem::measure::gas_oil_ratio;

                const auto& t  = pvto[tableID].getSaturatedTable();
                const auto& Rs = t.getColumn(0);

                const auto numActRows = Rs.size();

                std::transform(std::begin(Rs), std::end(Rs),
                               linTable.column(tableID, primID, 0),
                    [&units](const double rs) -> double
                {
                    return units.from_si(uRs, rs);
                });

                return numActRows;
            });
        }

        /// Extract maximum effective composition nodes in PVTO table data
        ///
        /// \param[in] pvto Collection of PVTO tables for all PVT regions.
        ///
        /// \return Maximum number of active keys across all tables of \p
        ///    pvto.
        std::size_t maxNumCompNodes(const std::vector<Opm::PvtoTable>& pvto)
        {
            if (pvto.empty()) { return 0; }

            return std::accumulate(std::begin(pvto), std::end(pvto), std::size_t{0},
                [](const std::size_t n, const Opm::PvtoTable& t) -> std::size_t
            {
                return std::max(n, t.getSaturatedTable().numRows());
            });
        }

        /// Extract maximum effective pressure nodes in PVDO table data
        ///
        /// \param[in] pvdo Collection of PVDO tables for all PVT regions.
        ///
        /// \return Maximum number of table rows across all tables of \p
        ///    pvdo.
        std::size_t maxNumPressNodes(const Opm::TableContainer& pvdo)
        {
            using PVDO = ::Opm::PvdoTable;
            auto tabID = std::vector<std::size_t>(pvdo.size());

            std::iota(tabID.begin(), tabID.end(), std::size_t{0});

            return std::accumulate(tabID.begin(), tabID.end(), std::size_t{0},
                [&pvdo](const std::size_t n, const std::size_t table)
            {
                return std::max(n, pvdo.getTable<PVDO>(table).numRows());
            });
        }

        /// Extract maximum effective pressure nodes in PVTO table data
        ///
        /// \param[in] pvto Collection of PVTO tables for all PVT regions.
        ///
        /// \return Maximum number of active pressure rows across all
        ///    sub-tables of \p pvto.
        std::size_t maxNumPressNodes(const std::vector<Opm::PvtoTable>& pvto)
        {
            auto max_nppvt = std::size_t{0};

            for (const auto& table : pvto) {
                const auto max_nppvt_tab =
                    std::accumulate(table.begin(), table.end(), std::size_t{0},
                        [](const std::size_t n, const Opm::SimpleTable& t)
                    {
                        return std::max(n, t.numRows());
                    });

                max_nppvt = std::max(max_nppvt, max_nppvt_tab);
            }

            return max_nppvt;
        }
    } // Oil

    /// Functions to create linearised, padded, and normalised water PVT
    /// output tables from input water PVT function keyword.
    namespace Water {
        /// Create linearised 'TAB' vector entries of normalised water PVT
        /// tables for all PVT function regions from PVTW keyword data.
        ///
        /// \param[in] units Active unit system.  Needed to convert SI
        ///    convention pressure values (Pascal), formation volume factors
        ///    (Rm3/Sm3), compressibility values (1/Pascal), viscosity
        ///    values (Pa*s), and viscosibility values (1/Pascal) to
        ///    declared conventions of the run specification.
        ///
        /// \param[in] pvtw Collection of PVTW tables for all PVT regions.
        ///
        /// \return Linearised 'TAB' vector values for output water PVT
        ///    tables.  A unit-converted, transformed version of the input
        ///    table \p pvtw.  No derivative information added.
        std::vector<double>
        fromPVTW(const Opm::UnitSystem& units,
                 const Opm::PvtwTable&  pvtw)
        {
            // Recall: PvtwTable is essentially vector<PVTWRecord> in which
            //
            //   struct PVTWRecord {
            //       double reference_pressure;
            //       double volume_factor;
            //       double compressibility;
            //       double viscosity;
            //       double viscosibility;
            //   };

            using M = ::Opm::UnitSystem::measure;

            // Columns [ Pw, 1/Bw, Cw, 1/(Bw*mu_w), Cw - Cv ]
            //
            // Single row per table.  No derivatives.  Can't reuse
            // createPropfuncTable here, so implement the return value
            // directly in terms of LinearisedOutputTable.

            const auto numTab  = pvtw.size();
            const auto numPrim = std::size_t{1};
            const auto numRows = std::size_t{1};
            const auto numCols = std::size_t{5};

            auto lintable = ::Opm::LinearisedOutputTable {
                numTab, numPrim, numRows, numCols
            };

            // Note unit hack for compressibility and viscosibility.  The
            // unit of measurement for these quantities is 1/pressure, but
            // the UnitSystem does not define this unit.  Work around the
            // missing conversion by using *to_si()* rather than *from_si()*
            // for those quantities.

            const auto uPress    = M::pressure;
            const auto uRecipFVF = M::water_inverse_formation_volume_factor;
            const auto uVisc     = M::viscosity;

            const auto primID = std::size_t{0};
            for (auto tabID = 0*numTab; tabID < numTab; ++tabID) {
                const auto& t = pvtw[tabID];

                auto iPw           = lintable.column(tabID, primID, 0);
                auto irecipFvf     = lintable.column(tabID, primID, 1);
                auto iCw           = lintable.column(tabID, primID, 2);
                auto irecipFvfVisc = lintable.column(tabID, primID, 3);
                auto idiffCwCv     = lintable.column(tabID, primID, 4);

                *iPw       = units.from_si(uPress   , t.reference_pressure);
                *irecipFvf = units.from_si(uRecipFVF, 1.0 / t.volume_factor);

                // Compressibility unit hack here (*to_si()*)
                *iCw = units.to_si(uPress, t.compressibility);

                *irecipFvfVisc =
                    units.from_si(uRecipFVF, 1.0 / t.volume_factor)
                    / units.from_si(uVisc, t.viscosity);

                // Viscosibility unit hack here (*to_si()*)
                *idiffCwCv =
                    units.to_si(uPress, t.compressibility - t.viscosibility);
            }

            return lintable.getDataDestructively();
        }
    } // Water

}} // Anonymous::PVTFunc

namespace Opm {

    Tables::Tables(const UnitSystem& units)
        : units_  (units)
        , tabdims_(Ix::TabDimsNumElems, 0)
    {
        // Initialize subset of base pointers and dimensions to 1 to honour
        // requirements of TABDIMS protocol.
        std::fill_n(std::begin(this->tabdims_), 59, 1);
    }

    void Tables::addDensity(const DensityTable& density)
    {
        if (density.size() == 0) { return; }

        const auto default_value = -2.0e20;
        const auto nreg = density.size();
        const auto nph  = density[0].size;

        std::vector<double> densityData(nreg * nph, default_value);
        {
            const auto urho = UnitSystem::measure::density;

            auto* rho_o = &densityData[0*nreg + 0]; // Oil   <-> Column 0
            auto* rho_w = &densityData[1*nreg + 0]; // Water <-> Column 1
            auto* rho_g = &densityData[2*nreg + 0]; // Gas   <-> Column 2

            for (const auto& record : density) {
                *rho_o++ = this->units_.from_si(urho, record.oil);
                *rho_w++ = this->units_.from_si(urho, record.water);
                *rho_g++ = this->units_.from_si(urho, record.gas);
            }
        }

        this->addData(Ix::DensityTableStart, densityData);

        this->tabdims_[Ix::DensityNumTables] = nreg;
    }

    void Tables::addPVTTables(const EclipseState& es)
    {
        const auto& phases = es.runspec().phases();

        if (phases.active(Phase::GAS))   { this->addGasPVTTables  (es); }
        if (phases.active(Phase::OIL))   { this->addOilPVTTables  (es); }
        if (phases.active(Phase::WATER)) { this->addWaterPVTTables(es); }
    }

    void Tables::addSatFunc(const EclipseState& es)
    {
        const auto active = SatFunc::ActivePhases { es };
        const auto family = SatFunc::Family       { es, active };

        if (family.ambiguous()) {
            // No known saturation function family or at least two of the
            // families I, II, or III.  Can't have that.

            return; // Logging here?
        }

        const auto nssfun = es.runspec()
            .tabdims().getNumSatNodes();

        const auto tolcrit = es.runspec()
            .saturationFunctionControls()
            .minimumRelpermMobilityThreshold();

        if (active.gas) {
            this->addSatFuncGas(SatFunc::Gas::normalisedTable(family, nssfun,
                                                              tolcrit, this->units_,
                                                              es.getTableManager()));
        }

        if (active.oil) {
            this->addSatFuncOil(SatFunc::Oil::normalisedTable(family, active,
                                                              nssfun, tolcrit,
                                                              es.getTableManager()));
        }

        if (active.wat) {
            this->addSatFuncWater(SatFunc::Water::normalisedTable(family, nssfun,
                                                                  tolcrit, this->units_,
                                                                  es.getTableManager()));
        }
    }

    const std::vector<int>& Tables::tabdims() const
    {
        return this->tabdims_;
    }

    const std::vector<double>& Tables::tab() const
    {
        return this->data_;
    }

    // -----------------------------------------------------------------------
    // Private implementation functions below separator.
    // -----------------------------------------------------------------------

    void Tables::addData(const std::size_t          offset_index,
                         const std::vector<double>& new_data)
    {
        // Add one for INIT file's one-based indexing protocol.
        this->tabdims_[offset_index] = this->data_.size() + 1;

        this->data_.insert(this->data_.end(), new_data.begin(), new_data.end());

        // Track number of elements in TAB vector.
        this->tabdims_[Ix::TabSize] = this->data_.size();
    }

    void Tables::addSatFuncGas(const std::tuple<int, int, std::vector<double>>& sgfn)
    {
        this->addData(Ix::SgfnTableStart, std::get<std::vector<double>>(sgfn));

        this->tabdims_[Ix::SgfnNumSatNodes] = std::get<0>(sgfn);
        this->tabdims_[Ix::SgfnNumTables]   = std::get<1>(sgfn);
    }

    void Tables::addSatFuncOil(const std::tuple<int, int, std::vector<double>>& sofn)
    {
        this->addData(Ix::SofnTableStart, std::get<std::vector<double>>(sofn));

        this->tabdims_[Ix::SofnNumSatNodes] = std::get<0>(sofn);
        this->tabdims_[Ix::SofnNumTables]   = std::get<1>(sofn);
    }

    void Tables::addSatFuncWater(const std::tuple<int, int, std::vector<double>>& swfn)
    {
        this->addData(Ix::SwfnTableStart, std::get<std::vector<double>>(swfn));

        this->tabdims_[Ix::SwfnNumSatNodes] = std::get<0>(swfn);
        this->tabdims_[Ix::SwfnNumTables]   = std::get<1>(swfn);
    }

    void Tables::addGasPVTTables(const EclipseState& es)
    {
        const auto& tabMgr = es.getTableManager();
        const auto& tabd   = es.runspec().tabdims();

        const auto numPressNodes = tabd.getNumPressureNodes();

        const auto hasPVTG = !tabMgr.getPvtgTables().empty();
        const auto hasPVDG =  tabMgr.hasTables("PVDG");

        if (hasPVTG + hasPVDG != 1) {
            // Inconsistent table specification.  Maybe throw here?
            return;
        }

        if (hasPVTG) {
            // Wet gas with vaporised/volatile oil.
            const auto& pvtg = tabMgr.getPvtgTables();

            const auto numCompNodes =
                std::max(tabd.getNumRSNodes(), PVTFunc::Gas::maxNumCompNodes(pvtg));

            const auto numPrimary =
                std::max(numPressNodes, PVTFunc::Gas::maxNumPressNodes(pvtg));

            const auto tableData = PVTFunc::Gas::
                fromPVTG(numCompNodes, numPrimary, this->units_, pvtg);

            const auto pressData = PVTFunc::Gas::
                pressureNodes(numPrimary, this->units_, pvtg);

            this->addData(Ix::PvtgMainStart, tableData);
            this->addData(Ix::PvtgPressStart, pressData);

            this->tabdims_[Ix::NumPvtgPressNodes] = numPrimary;
            this->tabdims_[Ix::NumPvtgCompNodes]  = numCompNodes;
            this->tabdims_[Ix::NumPvtgTables]     = pvtg.size();
        }
        else {
            // Dry gas, pressure dependent compressibility.
            const auto& pvdg = tabMgr.getPvdgTables();

            const auto numRows =
                std::max(numPressNodes, PVTFunc::Gas::maxNumPressNodes(pvdg));

            const auto tableData = PVTFunc::Gas::fromPVDG(numRows, this->units_, pvdg);

            this->addData(Ix::PvtgMainStart, tableData);

            this->tabdims_[Ix::NumPvtgPressNodes] = numRows;
            this->tabdims_[Ix::NumPvtgTables]     = pvdg.size();
        }
    }

    void Tables::addOilPVTTables(const EclipseState& es)
    {
        const auto& tabMgr = es.getTableManager();
        const auto& tabd   = es.runspec().tabdims();

        const auto numPressNodes = tabd.getNumPressureNodes();

        const auto hasPVTO  = !tabMgr.getPvtoTables().empty();
        const auto hasPVDO  =  tabMgr.hasTables("PVDO");
        const auto hasPVCDO = !tabMgr.getPvcdoTable().empty();

        if (hasPVTO + hasPVDO + hasPVCDO != 1) {
            // Inconsistent table specification.  Maybe throw here?
            return;
        }

        if (hasPVTO) {
            // Live oil with dissolved gas.
            const auto& pvto = tabMgr.getPvtoTables();

            const auto numCompNodes =
                std::max(tabd.getNumRSNodes(), PVTFunc::Oil::maxNumCompNodes(pvto));

            const auto numRows =
                std::max(numPressNodes, PVTFunc::Oil::maxNumPressNodes(pvto));

            const auto tableData = PVTFunc::Oil::
                fromPVTO(numCompNodes, numRows, this->units_, pvto);

            const auto rsData = PVTFunc::Oil::
                compositionNodes(numCompNodes, this->units_, pvto);

            this->addData(Ix::PvtoMainStart, tableData);
            this->addData(Ix::PvtoCompStart, rsData);

            this->tabdims_[Ix::NumPvtoPressNodes] = numRows;
            this->tabdims_[Ix::NumPvtoCompNodes]  = numCompNodes;
            this->tabdims_[Ix::NumPvtoTables]     = pvto.size();
        }
        else if (hasPVDO) {
            // Dead oil, pressure dependent compressibility.
            const auto& pvdo = tabMgr.getPvdoTables();

            const auto numRows =
                std::max(numPressNodes, PVTFunc::Oil::maxNumPressNodes(pvdo));

            const auto tableData = PVTFunc::Oil::fromPVDO(numRows, this->units_, pvdo);

            this->addData(Ix::PvtoMainStart, tableData);

            this->tabdims_[Ix::NumPvtoPressNodes] = numRows;
            this->tabdims_[Ix::NumPvtoTables]     = pvdo.size();
        }
        else {
            // Dead oil, constant compressibility.
            const auto& pvcdo = tabMgr.getPvcdoTable();

            const auto numRows = std::max(numPressNodes, pvcdo.size());

            const auto tableData = PVTFunc::Oil::fromPVCDO(numRows, this->units_, pvcdo);

            this->addData(Ix::PvtoMainStart, tableData);

            this->tabdims_[Ix::NumPvtoPressNodes] = numRows;
            this->tabdims_[Ix::NumPvtoTables]     = pvcdo.size();
        }
    }

    void Tables::addWaterPVTTables(const EclipseState& es)
    {
        const auto& tabMgr = es.getTableManager();

        const auto& pvtw = tabMgr.getPvtwTable();

        if (pvtw.empty()) {
            return;
        }

        const auto tableData = PVTFunc::Water::fromPVTW(this->units_, pvtw);

        this->addData(Ix::PvtwStart, tableData);

        this->tabdims_[Ix::NumPvtwTables] = pvtw.size();
    }

} // namespace Opm
