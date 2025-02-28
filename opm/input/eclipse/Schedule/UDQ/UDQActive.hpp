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

#ifndef UDQ_USAGE_HPP
#define UDQ_USAGE_HPP

#include <opm/input/eclipse/EclipseState/Phase.hpp>

#include <opm/input/eclipse/Schedule/UDQ/UDQEnums.hpp>

#include <opm/input/eclipse/Deck/UDAValue.hpp>

#include <cstddef>
#include <cstdlib>
#include <optional>
#include <string>
#include <vector>

namespace Opm {

class UDQConfig;
class UnitSystem;

} // namespace Opm

namespace Opm::RestartIO {
    struct RstState;
} // namespace Opm::RestartIO

namespace Opm {

/// Internalised representation of all UDAs in a simulation run
class UDQActive
{
public:
    /// Single UDA created from restart file information
    struct RstRecord
    {
        /// Constructor
        ///
        /// Creates a general UDA from restart file information.
        ///
        /// \param[in] control_arg Which item/limit of which constraint
        /// keyword (e.g., WCONPROD or GCONPROD) for which this UDA supplies
        /// the numeric value.
        ///
        /// \param[in] value_arg UDA value loaded from restart file.
        /// Typically a UDQ name with associate unit conversion operators.
        ///
        /// \param[in] wgname_arg Well or group name affected by this UDA.
        RstRecord(const UDAControl   control_arg,
                  const UDAValue     value_arg,
                  const std::string& wgname_arg)
            : control { control_arg }
            , value   { value_arg }
            , wgname  { wgname_arg }
        {}

        /// Constructor
        ///
        /// Creates a group level UDA for an injection limit.
        ///
        /// \param[in] control_arg Which item/limit of GCONINJE for which
        /// this UDA supplies the numeric value.
        ///
        /// \param[in] value_arg UDA value loaded from restart file.
        /// Typically a UDQ name with associate unit conversion operators.
        ///
        /// \param[in] wgname_arg Group name affected by this UDA.
        ///
        /// \param[in] phase Injected phase.
        RstRecord(const UDAControl   control_arg,
                  const UDAValue     value_arg,
                  const std::string& wgname_arg,
                  const Phase        phase)
            : RstRecord { control_arg, value_arg, wgname_arg }
        {
            this->ig_phase = phase;
        }

        /// Item/limit of constraint keyword for which this UDA supplies the
        /// numeric value.
        UDAControl control;

        /// UDA value.
        ///
        /// Typically a UDQ name and associate unit conversion operators.
        UDAValue value;

        /// Name of well/group affected by this UDA.
        std::string wgname;

        /// Injected phase in group level injection.
        ///
        /// Nullopt unless the control is a GCONINJE item.
        std::optional<Phase> ig_phase{};
    };

    /// Single UDA with use counts and IUAP start offsets for restart file
    /// output purposes.
    ///
    /// This information is intended to go mostly unaltered into the IUAD
    /// restart file array.
    class OutputRecord
    {
    public:
        /// Default constructor.
        ///
        /// Creates an invalid OutputRecord which is mostly usable as a
        /// target for a deserialisation operation.
        OutputRecord();

        /// Constructor
        ///
        /// \param[in] udq_arg Name of UDQ from which this UDA derives its
        /// numeric value.
        ///
        /// \param[in] input_index_arg Zero-based index in order of
        /// appearance of the UDQ use for this UDA.
        ///
        /// \param[in] use_index_arg IUAP array start offset.
        ///
        /// \param[in] wgname_arg Name of well or group affected by this
        /// UDA.
        ///
        /// \param[in] control_arg Constraint keyword and item/limit for
        /// which this UDA supplies the numeric value.
        OutputRecord(const std::string& udq_arg,
                     const std::size_t  input_index_arg,
                     const std::string& wgname_arg,
                     const UDAControl   control_arg);

        /// Equality predicate.
        ///
        /// \param[in] other Object against which \code *this \endcode will be
        /// tested for equality.
        ///
        /// \return Whether or not \code *this \endcode is the same as \p
        /// other.
        bool operator==(const OutputRecord& other) const;

        /// Inequality predicate.
        ///
        /// \param[in] other Object against which \code *this \endcode will be
        /// tested for inequality.
        ///
        /// \return Whether or not \code *this \endcode is different from \p
        /// other.
        bool operator!=(const OutputRecord& other) const
        {
            return ! (*this == other);
        }

        /// Convert between byte array and object representation.
        ///
        /// \tparam Serializer Byte array conversion protocol.
        ///
        /// \param[in,out] serializer Byte array conversion object.
        template<class Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(udq);
            serializer(input_index);
            serializer(wgname);
            serializer(control);
            serializer(uda_code);
            serializer(use_count);
        }

        /// Name of UDQ from which this UDA derives its numeric value.
        std::string udq;

        /// Zero-based index in order of appearance of the UDQ use for this
        /// UDA.
        std::size_t input_index;

        /// Constraint keyword and item/limit for which this UDA supplies
        /// the numeric value.
        UDAControl control;

        /// Restart file integer representation of \c control.
        int uda_code{};

        /// Name of well/group affected by this UDA.
        ///
        /// Misleading if the use count is greater than one.
        const std::string& wg_name() const
        {
            return this->wgname;
        }

        /// Number of times this UDA is mentioned in this particular
        /// combination of constraint keyword and item/limit.  Effectively,
        /// how many wells/groups use this UDA for the same purpose.
        std::size_t use_count{};

    private:
        /// Name of well/group affected by this UDA.
        ///
        /// Misleading if the use count is greater than one.
        std::string wgname{};
    };

    /// Internalised representation of a UDA from the input file
    class InputRecord
    {
    public:
        /// Default constructor.
        ///
        /// Creates an invalid input record that is mostly useful as a
        /// target for a deserialisation operation.
        InputRecord();

        /// Constructor.
        ///
        /// \param[in] input_index_arg Zero-based index in order of
        /// appearance of the UDQ used for this UDA.  Needed for restart
        /// file output purposes.
        ///
        /// \param[in] udq_arg Name of UDQ used for this UDA.
        ///
        /// \param[in] wgname_arg Name of well/group affected by this UDA.
        ///
        /// \param[in] control_arg Keyword and item/limit for which this UDA
        /// supplies the numeric value.
        InputRecord(const std::size_t  input_index_arg,
                    const std::string& udq_arg,
                    const std::string& wgname_arg,
                    const UDAControl   control_arg);

        /// Equality predicate.
        ///
        /// \param[in] other Object against which \code *this \endcode will be
        /// tested for equality.
        ///
        /// \return Whether or not \code *this \endcode is the same as \p
        /// other.
        bool operator==(const InputRecord& other) const;

        /// Convert between byte array and object representation.
        ///
        /// \tparam Serializer Byte array conversion protocol.
        ///
        /// \param[in,out] serializer Byte array conversion object.
        template<class Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(input_index);
            serializer(udq);
            serializer(wgname);
            serializer(control);
        }

        /// Zero-based index in order of appearance of the UDQ use for this
        /// UDA.  Needed for restart file output purposes.
        std::size_t input_index;

        /// Name of the UDQ used in this UDA.
        std::string udq;

        /// Well or group affected by this UDA.
        std::string wgname;

        /// Constraint keyword and item/limit for which this UDA supplies
        /// the numeric value.
        UDAControl control;
    };

    /// Create a serialisation test object.
    static UDQActive serializationTestObject();

    /// Default constructor.
    ///
    /// Creates an empty collection of UDAs.  Resulting collection is usable
    /// as a target for a deserialisation operation, or as a container of
    /// new UDAs through calls to the update() member function.
    UDQActive() = default;

    /// Load UDAs from restart file.
    ///
    /// \param[in] units
    ///
    /// \param[in] udq_config Simulation run's collection of user defined
    /// quantities.
    ///
    /// \param[in] rst_state Restart file information.
    ///
    /// \param[in] well_names Run's wells at restart point.
    ///
    /// \param[in] group_names Run's groups at restart point.
    static std::vector<RstRecord>
    load_rst(const UnitSystem&               units,
             const UDQConfig&                udq_config,
             const RestartIO::RstState&      rst_state,
             const std::vector<std::string>& well_names,
             const std::vector<std::string>& group_names);

    /// Amend collection of input UDAs to account for a new entry
    ///
    /// Does nothing if the UDA is numeric.  Adds a new record if none
    /// exists for the particular combination of (UDA, well/group name, and
    /// constraint keyword item).  Removes a record if the UDA value is
    /// numeric and previously used a UDQ specification.  Replaces a record
    /// if a different UDA was used for the same combination of well/group
    /// name and keyword/item.
    ///
    /// \param[in] udq_config Simulation run's collection of user defined
    /// quantities.
    ///
    /// \param[in] uda Numeric or string UDA value entered for a single
    /// limit/item in a constraint keyword.
    ///
    /// \param[in] wgname Well/group name affected by \p uda.
    ///
    /// \param[in] control Constraint keyword and associate item/limit for
    /// which \p uda supplies the numeric value.
    ///
    /// \return Whether or not internal data structures were altered.  One
    /// (1) if changes were made, and zero (0) otherwise.
    int update(const UDQConfig&   udq_config,
               const UDAValue&    uda,
               const std::string& wgname,
               const UDAControl   control);

    /// UDA existence predicate
    ///
    /// \return Whether or not there are any UDAs registered in this
    /// collection.
    explicit operator bool() const;

    /// Retrieve current set of UDAs, condensed by use counts and IUAP
    /// offsets.
    ///
    /// Intended for restart file output purposes only.
    const std::vector<OutputRecord>& iuad() const;

    /// Retrieve current set of UDAs from which to form IUAP restart file
    /// array.
    ///
    /// Intendend for restart file output purposes only.
    ///
    /// \note This function's role could possibly be served by iuad() as
    /// well.  If so, that's a future performance benefit since we won't
    /// have to form a new vector<> on every call to the function.
    std::vector<InputRecord> iuap() const;

    /// Equality predicate.
    ///
    /// \param[in] data Object against which \code *this \endcode will be
    /// tested for equality.
    ///
    /// \return Whether or not \code *this \endcode is the same as \p data.
    bool operator==(const UDQActive& data) const;

    /// Convert between byte array and object representation.
    ///
    /// \tparam Serializer Byte array conversion protocol.
    ///
    /// \param[in,out] serializer Byte array conversion object.
    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(input_data);
        serializer(output_data);
    }

private:
    /// Current set of UDAs entered in the input source.
    std::vector<InputRecord> input_data{};

    /// Current set of UDAs condensed by use counts and IUAP start pointers.
    ///
    /// Intended for restart file output as the IUAD vector.  Cleared if
    /// input_data changes.  Formed in constructOutputRecords().
    ///
    /// Marked 'mutable' because we may need to construct this from the
    /// 'const' iuad() member function.
    mutable std::vector<OutputRecord> output_data{};

    /// Form output_data structure from input_data.
    void constructOutputRecords() const;
};

} // namespace Opm

#endif // UDQ_USAGE_HPP
