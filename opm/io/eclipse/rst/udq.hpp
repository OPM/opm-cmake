/*
  Copyright 2021 Equinor ASA.

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

#ifndef RST_UDQ
#define RST_UDQ

#include <opm/input/eclipse/EclipseState/Phase.hpp>

#include <opm/common/utility/CSRGraphFromCoordinates.hpp>

#include <opm/input/eclipse/Schedule/UDQ/UDQEnums.hpp>

#include <cstddef>
#include <iterator>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

namespace Opm::RestartIO {

/// Container for a single user defined quantity (UDQ) reconstituted from
/// restart file information.
///
/// The producing side is expected to construct an RstUDQ object, and to
/// signal if the UDQ represents a quantity with a defining expression (four
/// parameter constructor) or without a defining expression (two parameter
/// constructor), the latter typically representing an assigned quantity
/// only.  Moreover, the producing side should call prepareValues() prior to
/// incorporating numeric values.  Then to assign the numeric values loaded
/// from a restart file using one of the functions assignScalarValue() or
/// addValue().  The former is intended for scalar UDQs, e.g., those at
/// field level, while the latter is intended for UDQ sets such as those
/// pertaining to wells, groups, or well segments.  Mixing
/// assignScalarValue() and addValue() for a single UDQ object will make
/// class RstUDQ throw an exception.  Once all values have been added, the
/// producing side is expected to call commitValues().
///
/// The consuming side, typically client code which forms \c UDQConfig and
/// \c UDQState objects, is expected to query the object for how to
/// interpret the values and then to call \code operator[]() \endcode to
/// retrieve a range of numeric values for the UDQ pertaining to a sing
/// entity.
class RstUDQ
{
public:
    /// Sequence of sub-entity IDs and values pertaining to a single entity.
    class ValueRange
    {
    public:
        /// Simple forward iterator over a UDQ set's values pertaining to a
        /// single entity (e.g., a well or a group).
        class Iterator
        {
        public:
            /// Iterator's category (forward iterator)
            using iterator_category = std::forward_iterator_tag;

            /// Iterator's value type
            ///
            /// Pair's \code .first \endcode element is the zero-based
            /// sub-entity index, while \code .second \endcode is the
            /// sub-entity's value.
            using value_type = std::pair<std::size_t, double>;

            /// Iterator's difference type.
            using difference_type = int;

            /// Iterator's pointer type (return type from operator->())
            using pointer = const value_type*;

            /// Iterator's reference type (return type from operator*())
            using reference = const value_type&;

            /// Pre-increment operator.
            ///
            /// \return *this.
            Iterator& operator++()
            {
                ++this->ix_;

                return *this;
            }

            /// Post-increment operator.
            ///
            /// \return Iterator pointing to element prior to increment result.
            Iterator operator++(int)
            {
                auto iter = *this;

                ++(*this);

                return iter;
            }

            /// Dereference operator
            ///
            /// \return Element at current position.
            reference operator*()
            {
                this->setValue();
                return this->deref_value_;
            }

            /// Indirection operator
            ///
            /// \return Pointer to element at current position.
            pointer operator->()
            {
                this->setValue();
                return &this->deref_value_;
            }

            /// Equality predicate
            ///
            /// \param[in] that Object to which \c *this will be compared for equality.
            ///
            /// \return Whether or not \c *this equals \p that.
            bool operator==(const Iterator& that) const
            {
                return (this->ix_    == that.ix_)
                    && (this->i_     == that.i_)
                    && (this->value_ == that.value_);
            }

            /// Inequality predicate
            ///
            /// \param[in] that Object to which \c *this will be compared for inequality.
            ///
            /// \return \code ! (*this == that)
            bool operator!=(const Iterator& that) const
            {
                return ! (*this == that);
            }

            friend class ValueRange;

        private:
            /// Constructor
            ///
            /// Accessible to RegionIndexRange only.
            ///
            /// \param[in] index range element value.
            Iterator(std::size_t ix,
                     const int*  i,
                     const std::variant<double, const double*>& value)
                : ix_{ix}, i_{i}, value_{value}
            {}

            /// Index range element value
            std::size_t ix_;
            const int* i_;
            std::variant<double, const double*> value_;

            value_type deref_value_{};

            void setValue();
        };

        /// Start of value range
        Iterator begin() const
        { return this->makeIterator(this->begin_); }

        /// End of value range
        Iterator end() const
        { return this->makeIterator(this->end_); }

        friend class RstUDQ;

    private:
        /// Constructor
        ///
        /// Scalar version.
        ///
        /// \param[in] begin_arg Index of start of range.
        ///
        /// \param[in] end_arg Index of element one past the end of the
        /// range.
        ///
        /// \param[in] i Range of sub-entity indices.
        ///
        /// \param[in] value Constant value applicable to every element in
        /// the range.
        ValueRange(std::size_t  begin_arg,
                   std::size_t  end_arg,
                   const int*   i,
                   const double value)
            : begin_{begin_arg}, end_{end_arg}, i_{i}, value_{value}
        {}

        /// Constructor
        ///
        /// Array version.
        ///
        /// \param[in] begin_arg Index of start of range.
        ///
        /// \param[in] end_arg Index of element one past the end of the
        /// range.
        ///
        /// \param[in] i Range of sub-entity indices.
        ///
        /// \param[in] value Range of values, one scalar value for each
        /// sub-entity in the value range.
        ValueRange(std::size_t   begin_arg,
                   std::size_t   end_arg,
                   const int*    i,
                   const double* value)
            : begin_{begin_arg}, end_{end_arg}, i_{i}, value_{value}
        {}

        /// Index of start of the range.
        std::size_t begin_{};

        /// Index of element one past the end of the range.
        std::size_t end_{};

        /// Sub-entities.
        const int* i_{};

        /// Values pertaining to each sub-entity of the range.
        ///
        /// Holds a \c double in the scalar version, and a \code const
        /// double* \end in the array version.
        std::variant<double, const double*> value_{};

        /// Form an iterator from a sequence index
        ///
        /// \param[in] ix Sequence index
        ///
        /// \return Range iterator pointing to position \p ix.
        Iterator makeIterator(const std::size_t ix) const
        {
            return { ix, this->i_, this->value_ };
        }
    };

    /// UDQ name.
    ///
    /// First argument to a DEFINE or an ASSIGN statement.
    std::string name{};

    /// UDQ's unit string.
    std::string unit{};

    /// UDQ's category, i.e., which level this UDQ applies to.
    ///
    /// Examples include, the FIELD (FU*), group (GU*), well (WU*), or well
    /// segments (SU*).
    UDQVarType category{UDQVarType::NONE};

    /// Constructor
    ///
    /// Quantity given by DEFINE statement.
    ///
    /// \param{in] name_arg UDQ name.
    ///
    /// \param{in] unit_arg UDQ's unit string.
    ///
    /// \param{in] define_arg UDQ's defining expression, previously entered
    /// in a DEFINE statement.
    ///
    /// \param{in] status_arg UDQ's update status, i.e., when to calculate
    /// new values for the quantity.
    RstUDQ(const std::string& name_arg,
           const std::string& unit_arg,
           const std::string& define_arg,
           UDQUpdate          status_arg);

    /// Constructor
    ///
    /// Quantity given by ASSIGN statement.
    ///
    /// \param{in] name_arg UDQ name.
    ///
    /// \param{in] unit_arg UDQ's unit string.
    RstUDQ(const std::string& name_arg,
           const std::string& unit_arg);

    /// Prepare inclusion of entities and values.
    ///
    /// Used by producing code, i.e., the layer which reads UDQ information
    /// directly from file.
    void prepareValues();

    /// Assign numeric UDQ value for a entity/sub-entity pair.
    ///
    /// Conflicts with assignScalarValue().  If assignScalarValue() has been
    /// called previously, then this function will throw an exception of
    /// type \code std::logic_error \endcode.
    ///
    /// \param[in] entity Numeric entity index.  Assumed to be zero-based
    /// linear index of an entity, such as a well or a group.
    ///
    /// \param[in] subEntity Numeric sub-entity index.  Assumed to be a
    /// zero-based linear index of a sub-entity such as a well segment.  For
    /// well or group level UDQs, pass zero for the sub-entity.
    ///
    /// \param[in] value Numeric UDQ value for this entity/sub-entity pair.
    void addValue(const int entity, const int subEntity, const double value);

    /// End value accumulation.
    ///
    /// Builds entity to sub-entity mapping, and sorts the corresponding
    /// numeric UDQ values.  Should normally be the final member function
    /// call on the producing side.
    void commitValues();

    /// Assign a scalar value for the UDQ.
    ///
    /// Conflicts with addValue().  If addValue() has been called
    /// previously, then this function will throw an exception of type \code
    /// std::logic_error \endcode.
    ///
    /// \param[in] value Single numeric value pertaining to all entities in
    /// this UDQ set.
    void assignScalarValue(const double value);

    /// Add a name for the last (or next) entity used in a call to addValue()
    ///
    /// \param[in] wgname Entity, typically well or group, name.
    void addEntityName(std::string_view wgname);

    /// Retrieve UDQ's scalar value.
    ///
    /// Throws a \code std::logic_error \endcode unless isScalar() returns
    /// \c true.
    ///
    /// Used by consuming code, i.e., the layer which forms \c UDQConfig and
    /// \c UDQState objects from restart file information.
    double scalarValue() const;

    /// Retrieve number of of entities known to this UDQ.
    std::size_t numEntities() const;

    /// Get read-only access to sub-entities and values associated to a
    /// single top-level entity.
    ///
    /// \param[in] i Entity index.  Must be in the range [0..numEntities()).
    ///
    /// \return Sequence of values pertaining to entity \p i.
    ValueRange operator[](const std::size_t i) const;

    /// Get sequence of UDQ's entity names.
    ///
    /// \return Sequence of entity names entered in earlier addEntityName()
    /// calls.  Order not guaranteed.
    const std::vector<std::string>& entityNames() const
    {
        return this->wgnames_;
    }

    /// Index map for entity names.
    ///
    /// The entity name of entity 'i' is
    /// \code
    ///   entityNames()[ nameIndex()[i] ]
    /// \endcode
    /// provided named entities are meaningful for this UDQ--i.e., if it
    /// pertains to the well, group, connection, or segment levels.
    const std::vector<int>& nameIndex() const;

    /// UDQ's defining expression
    ///
    /// Empty character sequence if this UDQ does not have a defining
    /// expression--i.e., if it was entered using ASSIGN statements.
    std::string_view definingExpression() const;

    /// UDQ's current update status.
    UDQUpdate currentUpdateStatus() const;

    /// Whether or not this UDQ has a defining expression
    ///
    /// Needed to properly reconstitute the UDQConfig object from restart
    /// file information.
    bool isDefine() const
    {
        return this->definition_.has_value();
    }

    /// Predicate for whether or not this UDQ is a scalar
    /// quantity--typically a FIELD level value.
    bool isScalar() const
    {
        return std::holds_alternative<double>(this->sa_);
    }

private:
    /// Entity mapping type.
    ///
    /// VertexID = int
    /// TrackCompressedIdx = true (need SA mapping)
    /// PermitSelfConnections = true (MS well 5 may have segment number 5).
    using Graph = utility::CSRGraphFromCoordinates<int, true, true>;

    /// Wrapper for a DEFINE expression
    struct Definition
    {
        /// Constructor
        ///
        /// Mostly to enable using optional<>::emplace().
        ///
        /// \param[in] expression_arg UDQ's defining expression
        ///
        /// \param[in] status_arg UDQ's update status.
        Definition(const std::string& expression_arg,
                   const UDQUpdate    status_arg)
            : expression { expression_arg }
            , status     { status_arg }
        {}

        /// UDQ's defining expression.
        std::string expression{};

        /// UDQ's update status.
        UDQUpdate status{};
    };

    /// Map entities to range of sub-entities.
    ///
    /// Examples include
    ///
    ///   -# Well and its segments (SU*)
    ///   -# Well and its connections (CU*)
    ///
    /// For field, group, and well level UDQs, the sub-entities are usually
    /// trivial (i.e., one sub-entity of numeric index zero), although one
    /// *might* choose to represent aquifer and block level UDQs as
    /// sub-entities of the FIELD.
    Graph entityMap_{};

    /// UDQ values.
    ///
    /// Single scalar if all entities in this UDQ share the same constant
    /// value--typically only for field level UDQs.  Range of scalars (\code
    /// vector<double> \endcode) if each entity/sub-entity potentially has a
    /// different value.  Monostate unless any values have been assigned.
    std::variant<std::monostate, double, std::vector<double>> sa_;

    /// Entity names.
    ///
    /// Typically well or group names.
    std::vector<std::string> wgnames_{};

    /// Largest entity index seen in all addValue() calls so far.
    ///
    /// Nullopt before first call to addValue().
    std::optional<int> maxEntityIdx_{};

    /// Map entity indices to entity names.
    mutable std::optional<std::vector<int>> wgNameIdx_{};

    /// UDQ's definition.
    ///
    /// Nullopt unless this UDQ has a defining expression (four parameter
    /// constructor used to form the object).
    std::optional<Definition> definition_{};

    /// Whether or not this UDQ represents a non-scalar UDQ set.
    bool isUDQSet() const
    {
        return std::holds_alternative<std::vector<double>>(this->sa_);
    }

    /// Construct a value range for a set of sub-entities that share a
    /// common, constant numeric value (\c sa_ is \c double).
    ///
    /// \param[in] i Entity index.
    ///
    /// \return Value range for entity \p i.
    ValueRange scalarRange(const std::size_t i) const;

    /// Construct a value range for a set of entities whose numeric values
    /// may differ between sub-entities (\c sa_ is \code vector<double>
    /// \endcode).
    ///
    /// \param[in] i Entity index.
    ///
    /// \return Value range for entity \p i.
    ValueRange udqSetRange(const std::size_t i) const;

    void ensureValidNameIndex() const;
};

/// Collection of UDAs loaded from restart file
struct RstUDQActive
{
    /// One single UDA
    struct RstRecord
    {
        /// Constructor.
        ///
        /// Convenience only as this enables using vector<>::emplace_back().
        ///
        /// \param[in] c Control keyword and item
        ///
        /// \param[in] i Input index.  Zero-based order in which the UDQ was
        /// entered.
        ///
        /// \param[in] numIuap Number of IUAP elements associated to this
        /// UDA.
        ///
        /// \param[in] u1 Use count.  Number of times this UDA is used in
        /// this particular way (same keyword and control item, e.g., for
        /// different wells or groups).
        ///
        /// \param[in] u2 IUAP start offset.
        RstRecord(UDAControl  c,
                  std::size_t i,
                  std::size_t numIuap,
                  std::size_t u1,
                  std::size_t u2);

        /// Control keyword and associate item for this UDA.
        UDAControl control;

        /// Input index.  Zero-based order in which the UDQ was entered.
        std::size_t input_index;

        /// Use count.  Number of times this UDA is used in this particular
        /// way (same keyword and control item, e.g., for different wells or
        /// groups).
        std::size_t use_count;

        /// IUAP start offset.
        std::size_t wg_offset;

        /// Number of IUAP elements.
        std::size_t num_wg_elems;
    };

    /// Default constructor.
    ///
    /// Represents an empty collection.
    RstUDQActive() = default;

    /// Constructor.
    ///
    /// Forms UDA collection from restart file information.
    ///
    /// \param[in] iuad.  Restart file IUAD array.
    ///
    /// \param[in] iuap.  Restart file IUAP array.  Wells/groups affected by
    /// each UDA.
    ///
    /// \param[in] igph.  Restart file IGPH array.  Injection phases for
    /// groups.
    RstUDQActive(const std::vector<int>& iuad,
                 const std::vector<int>& iuap,
                 const std::vector<int>& igph);

    /// Wells/groups affected by each UDA.
    std::vector<int> wg_index{};

    /// Exploded items of each UDA.
    std::vector<RstRecord> iuad{};

    /// Injection phases for groups.
    std::vector<Phase> ig_phase{};
};

} // namespace Opm::RestartIO

#endif // RST_UDQ
