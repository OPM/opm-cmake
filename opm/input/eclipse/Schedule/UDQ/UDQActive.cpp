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

#include <opm/input/eclipse/Schedule/UDQ/UDQActive.hpp>

#include <opm/io/eclipse/rst/state.hpp>
#include <opm/io/eclipse/rst/udq.hpp>

#include <opm/input/eclipse/EclipseState/Phase.hpp>

#include <opm/input/eclipse/Schedule/UDQ/UDQConfig.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQEnums.hpp>

#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <string>
#include <vector>

#include <fmt/format.h>

namespace {
    bool isFieldUDA(const std::size_t num_wg_elems,
                    const int* const  localWgIndex)
    {
        // Note: all_of(localWgIndex == 0) corresponds to all_of(IUAP == 1)
        // in the restart file.  See RstUDQActive constructor.
        return (num_wg_elems > 1)
            && std::all_of(localWgIndex, localWgIndex + num_wg_elems,
                           [](const int wgIdx) { return wgIdx == 0; });
    }

    template <typename RstUDARecord>
    void loadRstWellUDAs(const RstUDARecord&                     record,
                         const std::vector<int>&                 wgIndex,
                         const Opm::UDAValue&                    uda,
                         const std::vector<std::string>&         wellNames,
                         std::vector<Opm::UDQActive::RstRecord>& udaRecords)
    {
        const auto* localWgIndex = &wgIndex[record.wg_offset + 0];

        for (auto i = 0*record.use_count; i < record.use_count; ++i) {
            udaRecords.emplace_back(record.control, uda, wellNames[localWgIndex[i]]);
        }
    }

    template <typename RstUDARecord, typename WGIndexOp>
    void loadGroupRstUDA(const RstUDARecord&     record,
                         const std::vector<int>& wgIndex,
                         WGIndexOp&&             wgIndexOp)
    {
        const auto* localWgIndex = &wgIndex[record.wg_offset + 0];

        for (auto i = 0*record.use_count; i < record.use_count; ++i) {
            wgIndexOp(record.num_wg_elems,
                      &localWgIndex[i*record.num_wg_elems + 0]);
        }
    }

    template <typename RstUDARecord>
    void loadRstGroupProdUDAs(const RstUDARecord&                     record,
                              const std::vector<int>&                 wgIndex,
                              const Opm::UDAValue&                    uda,
                              const std::vector<std::string>&         groupNames,
                              std::vector<Opm::UDQActive::RstRecord>& udaRecords)
    {
        loadGroupRstUDA(record, wgIndex,
                        [control = record.control,
                         &groupNames, &uda, &udaRecords]
                        (const std::size_t num_wg_elems,
                         const int* const  localWgIndex)
        {
            if (isFieldUDA(num_wg_elems, localWgIndex)) {
                udaRecords.emplace_back(control, uda, "FIELD");
            }
            else {
                udaRecords.emplace_back(control, uda, groupNames[localWgIndex[0] + 1]);
            }
        });
    }

    template <typename RstUDARecord>
    void loadRstGroupInjUDAs(const RstUDARecord&                     record,
                             const std::vector<int>&                 wgIndex,
                             const std::vector<Opm::Phase>&          igPhase,
                             const Opm::UDAValue&                    uda,
                             const std::vector<std::string>&         groupNames,
                             std::vector<Opm::UDQActive::RstRecord>& udaRecords)
    {
        loadGroupRstUDA(record, wgIndex,
                        [control = record.control,
                         &groupNames, &igPhase,
                         &uda, &udaRecords]
                        (const std::size_t num_wg_elems,
                         const int* const  localWgIndex)
        {
            const auto& group = isFieldUDA(num_wg_elems, localWgIndex)
                ? std::string { "FIELD" }
                : groupNames[localWgIndex[0] + 1];

            const auto phase = (group == "FIELD")
                ? igPhase.back()
                : igPhase[localWgIndex[0]];

            udaRecords.emplace_back(control, uda, group, phase);
        });
    }
} // Anonymous namespace

// ===========================================================================

Opm::UDQActive::OutputRecord::OutputRecord()
    : input_index { 0 }
    , control     { UDAControl::WCONPROD_ORAT }
    , uda_code    { 0 }
    , use_count   { 1 }
{}

Opm::UDQActive::OutputRecord::OutputRecord(const std::string& udq_arg,
                                           const std::size_t  input_index_arg,
                                           const std::string& wgname_arg,
                                           const UDAControl   control_arg)
    : udq         { udq_arg }
    , input_index { input_index_arg }
    , control     { control_arg }
    , uda_code    { UDQ::udaCode(control_arg) }
    , use_count   { 1 }
    , wgname      { wgname_arg }
{}

bool Opm::UDQActive::OutputRecord::operator==(const OutputRecord& other) const
{
    return (this->udq == other.udq)
        && (this->input_index == other.input_index)
        && (this->wgname == other.wgname)
        && (this->control == other.control)
        && (this->uda_code == other.uda_code)
        && (this->use_count == other.use_count)
        ;
}

// ---------------------------------------------------------------------------

Opm::UDQActive::InputRecord::InputRecord()
    : input_index { 0 }
    , control     { UDAControl::WCONPROD_ORAT }
{}

Opm::UDQActive::InputRecord::InputRecord(const std::size_t  input_index_arg,
                                         const std::string& udq_arg,
                                         const std::string& wgname_arg,
                                         const UDAControl   control_arg)
    : input_index { input_index_arg }
    , udq         { udq_arg }
    , wgname      { wgname_arg }
    , control     { control_arg }
{}

bool Opm::UDQActive::InputRecord::operator==(const InputRecord& other) const
{
    return (this->input_index == other.input_index)
        && (this->udq == other.udq)
        && (this->wgname == other.wgname)
        && (this->control == other.control)
        ;
}

// ---------------------------------------------------------------------------

Opm::UDQActive Opm::UDQActive::serializationTestObject()
{
    UDQActive result;

    result.input_data = {{1, "test1", "test2", UDAControl::WCONPROD_ORAT}};
    result.output_data = {{"test1", 1, "test2", UDAControl::WCONPROD_ORAT}};

    return result;
}

std::vector<Opm::UDQActive::RstRecord>
Opm::UDQActive::load_rst(const UnitSystem&               units,
                         const UDQConfig&                udq_config,
                         const RestartIO::RstState&      rst_state,
                         const std::vector<std::string>& well_names,
                         const std::vector<std::string>& group_names)
{
    auto udaRecords = std::vector<RstRecord> {};

    const auto& rst_active = rst_state.udq_active;
    const auto& wgIndex = rst_active.wg_index;

    for (const auto& record : rst_active.iuad) {
        const auto uda = UDAValue {
            udq_config[record.input_index].keyword(),
            units.uda_dim(record.control)
        };

        if (UDQ::well_control(record.control)) {
            loadRstWellUDAs(record, wgIndex, uda, well_names, udaRecords);
        }
        else if (UDQ::is_group_production_control(record.control)) {
            loadRstGroupProdUDAs(record, wgIndex, uda, group_names, udaRecords);
        }
        else {
            loadRstGroupInjUDAs(record, wgIndex, rst_active.ig_phase,
                                uda, group_names, udaRecords);
        }
    }

    return udaRecords;
}

Opm::UDQActive::operator bool() const
{
    return ! this->input_data.empty();
}

// We go through the current list of input records and compare with the
// supplied arguments (uda, wgname, control).  There are seven possible
// outcomes:
//
//   0. The uda variable is not defined (defaulted target/limit): fast return.
//
//   1. The uda variable is a double and no uda usage has been registered so far:
//      fast return.
//
//   2. The uda variable is a double, and the (wgname,control) combination
//      is found in the input data; this implies that uda has been used
//      previously for this particular (wgname, control) combination: We
//      remove that record from the input_data.
//
//   3. The uda variable is a string and we find that particular (udq,
//      wgname, control) combination in the input data: No changes
//
//   4. The uda variable is a string; but another udq was used for this
//      (wgname, control) combination: We replace the previous entry.
//
//   5. The uda variable is a string and we do not find this (wgname,
//      control) combination in the previous records: We add a new record.
//
//   6. The uda variable is a double, and the (wgname, control) combination
//      has not been encountered before: return 0
int Opm::UDQActive::update(const UDQConfig&   udq_config,
                           const UDAValue&    uda,
                           const std::string& wgname,
                           const UDAControl   control)
{
    // Alternative 0
    if (!uda.is_defined()) {
        return 0;
    }

    const auto isString = uda.is<std::string>();
    const auto quantity = isString ? uda.get<std::string>() : std::string{};

    // Alternative 1
    if (!isString && this->input_data.empty()) {
        return 0;
    }

    if (isString && !udq_config.has_keyword(quantity)) {
        throw std::logic_error {
            fmt::format("User defined quantity {} is not known and cannot "
                        "be used as a user defined argument in {} for {}\n"
                        "Missing ASSIGN or DEFINE for this UDQ?",
                        quantity, UDQ::controlName(control), wgname)
        };
    }

    for (auto record = this->input_data.begin();
         record != this->input_data.end(); ++record)
    {
        if ((record->wgname != wgname) || (record->control != control)) {
            continue;
        }

        if (isString && (record->udq == quantity)) {
            // Alternative 3
            return 0;
        }

        this->input_data.erase(record);
        this->output_data.clear();

        if (isString) {
            // Alternative 4
            break;
        }

        // Alternative 2
        return 1;
    }

    // Alternatives 4 & 5
    if (isString) {
        const auto udq_index = udq_config[quantity].index.insert_index;

        this->input_data.emplace_back(udq_index, quantity, wgname, control);
        this->output_data.clear();

        return 1;
    }

    // Alternative 6
    return 0;
}

const std::vector<Opm::UDQActive::OutputRecord>& Opm::UDQActive::iuad() const
{
    if (this->output_data.empty()) {
        this->constructOutputRecords();
    }

    return this->output_data;
}

std::vector<Opm::UDQActive::InputRecord> Opm::UDQActive::iuap() const
{
    auto iuap_data = std::vector<UDQActive::InputRecord>{};

    auto input_rcpy = this->input_data;
    while (! input_rcpy.empty()) {
        // Store next active control (new control).
        //
        // Recall: We must not take a reference to the new element here
        // since push_back() below may reallocate.  In that case, all
        // element references are invalidated.
        const auto cur_rec =
            iuap_data.emplace_back(input_rcpy.front());

        // Find and store active controls with same control and UDQ.
        auto it = input_rcpy.erase(input_rcpy.begin());
        while (it != input_rcpy.end()) {
            if ((it->control == cur_rec.control) && (it->udq == cur_rec.udq)) {
                iuap_data.push_back(*it);
                it = input_rcpy.erase(it);
            }
            else {
                ++it;
            }
        }
    }

    return iuap_data;
}

bool Opm::UDQActive::operator==(const UDQActive& data) const
{
    return (this->input_data  == data.input_data)
        && (this->output_data == data.output_data)
        ;
}

// ---------------------------------------------------------------------------
// Private member functions of class Opm::UDAQactive
// ---------------------------------------------------------------------------

void Opm::UDQActive::constructOutputRecords() const
{
    for (const auto& input_record : this->input_data) {
        auto it = std::find_if(this->output_data.begin(), this->output_data.end(),
                               [&input_record](const auto& output_record)
                               {
                                   return (output_record.udq     == input_record.udq)
                                       && (output_record.control == input_record.control);
                               });

        if (it != this->output_data.end()) {
            ++it->use_count;
        }
        else {
            // Recall: Constructor gives use_count = 1 in this case.
            this->output_data.emplace_back(input_record.udq,
                                           input_record.input_index,
                                           input_record.wgname,
                                           input_record.control);
        }
    }
}
