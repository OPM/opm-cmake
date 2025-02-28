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

#include <opm/io/eclipse/rst/udq.hpp>

#include <opm/common/utility/Visitor.hpp>

#include <opm/output/eclipse/UDQDims.hpp>
#include <opm/output/eclipse/VectorItems/udq.hpp>

#include <algorithm>
#include <cstddef>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include <fmt/format.h>

void Opm::RestartIO::RstUDQ::ValueRange::Iterator::setValue()
{
    this->deref_value_.first = this->i_[this->ix_];

    this->deref_value_.second =
        std::visit(VisitorOverloadSet {
            []              (const double  v) { return v; },
            [ix = this->ix_](const double* v) { return v[ix]; }
        }, this->value_);
}

// ---------------------------------------------------------------------------

Opm::RestartIO::RstUDQ::RstUDQ(const std::string& name_arg,
                               const std::string& unit_arg,
                               const std::string& define_arg,
                               const UDQUpdate    update_arg)
    : name       { name_arg }
    , unit       { unit_arg }
    , category   { UDQ::varType(name_arg) }
    , definition_{ std::in_place, define_arg, update_arg }
{}

Opm::RestartIO::RstUDQ::RstUDQ(const std::string& name_arg,
                               const std::string& unit_arg)
    : name     { name_arg }
    , unit     { unit_arg }
    , category { UDQ::varType(name_arg) }
{}

void Opm::RestartIO::RstUDQ::prepareValues()
{
    this->entityMap_.clear();
}

void Opm::RestartIO::RstUDQ::
addValue(const int entity, const int subEntity, const double value)
{
    if (this->isScalar()) {
        throw std::logic_error {
            fmt::format("UDQ {} cannot be defined as a scalar "
                        "and then used as UDQ set at restart time",
                        this->name)
        };
    }

    if (std::holds_alternative<std::monostate>(this->sa_)) {
        this->sa_.emplace<std::vector<double>>();
    }

    this->entityMap_.addConnection(entity, subEntity);
    std::get<std::vector<double>>(this->sa_).push_back(value);

    this->maxEntityIdx_ = ! this->maxEntityIdx_.has_value()
        ? entity
        : std::max(entity, *this->maxEntityIdx_);
}

void Opm::RestartIO::RstUDQ::commitValues()
{
    if (! this->isUDQSet()) {
        // Scalar or monostate.  Nothing to do.
        return;
    }

    // If we get here this is a UDQ set.  Form compressed UDQ mapping and
    // reorder the values to match this compressed mapping.

    this->entityMap_.compress(this->maxEntityIdx_.value() + 1);

    auto newSA = std::vector<double>(this->entityMap_.numEdges(), 0.0);
    {
        const auto& currSA = std::get<std::vector<double>>(this->sa_);

        const auto& ix = this->entityMap_.compressedIndexMap();

        for (auto i = 0*currSA.size(); i < currSA.size(); ++i) {
            newSA[ix[i]] += currSA[i];
        }
    }

    std::get<std::vector<double>>(this->sa_).swap(newSA);
}

void Opm::RestartIO::RstUDQ::assignScalarValue(const double value)
{
    if (this->isUDQSet()) {
        throw std::logic_error {
            fmt::format("UDQ {} cannot be defined as a UDQ set "
                        "and then used as a scalar at restart time",
                        this->name)
        };
    }

    this->sa_ = value;
}

void Opm::RestartIO::RstUDQ::addEntityName(std::string_view wgname)
{
    this->wgnames_.emplace_back(wgname);
}

double Opm::RestartIO::RstUDQ::scalarValue() const
{
    if (std::holds_alternative<std::monostate>(this->sa_)) {
        throw std::logic_error {
            fmt::format("Cannot request scalar value from UDQ "
                        "{} when no values have been assigned",
                        this->name)
        };
    }

    if (! this->isScalar()) {
        throw std::logic_error {
            fmt::format("Cannot request scalar value "
                        "from non-scalar UDQ set {}",
                        this->name)
        };
    }

    return std::get<double>(this->sa_);
}

std::size_t Opm::RestartIO::RstUDQ::numEntities() const
{
    return this->maxEntityIdx_.has_value()
        ? static_cast<std::size_t>(1 + *this->maxEntityIdx_)
        : this->wgnames_.size();
}

Opm::RestartIO::RstUDQ::ValueRange
Opm::RestartIO::RstUDQ::operator[](const std::size_t i) const
{
    if (std::holds_alternative<std::monostate>(this->sa_)) {
        throw std::logic_error {
            fmt::format("Cannot request values for "
                        "entity {} from UDQ {} when "
                        "no values have been assigned",
                        i, this->name)
        };
    }

    if (this->isScalar()) {
        return this->scalarRange(i);
    }

    return this->udqSetRange(i);
}

Opm::UDQUpdate Opm::RestartIO::RstUDQ::currentUpdateStatus() const
{
    return this->isDefine()
        ? this->definition_->status
        : UDQUpdate::OFF;
}

const std::vector<int>& Opm::RestartIO::RstUDQ::nameIndex() const
{
    this->ensureValidNameIndex();

    return *this->wgNameIdx_;
}

std::string_view Opm::RestartIO::RstUDQ::definingExpression() const
{
    return this->isDefine()
        ? std::string_view { this->definition_->expression }
        : std::string_view {};
}

Opm::RestartIO::RstUDQ::ValueRange
Opm::RestartIO::RstUDQ::scalarRange(const std::size_t i) const
{
    this->ensureValidNameIndex();

    return ValueRange {
        i, i + 1, this->wgNameIdx_->data(),
        std::get<double>(this->sa_)
    };
}

Opm::RestartIO::RstUDQ::ValueRange
Opm::RestartIO::RstUDQ::udqSetRange(const std::size_t i) const
{
    const auto& start = this->entityMap_.startPointers();
    const auto& cols  = this->entityMap_.columnIndices();
    const auto& sa    = std::get<std::vector<double>>(this->sa_);

    return ValueRange {
        start[i + 0], start[i + 1], cols.data(), sa.data()
    };
}

void Opm::RestartIO::RstUDQ::ensureValidNameIndex() const
{
    if (! this->wgNameIdx_.has_value()) {
        this->wgNameIdx_.emplace(this->wgnames_.size());

        std::iota(this->wgNameIdx_->begin(), this->wgNameIdx_->end(),
                  std::vector<std::string>::size_type{0});
    }
}

// ---------------------------------------------------------------------------

Opm::RestartIO::RstUDQActive::
RstRecord::RstRecord(const UDAControl  c,
                     const std::size_t i,
                     const std::size_t numIuap,
                     const std::size_t u1,
                     const std::size_t u2)
    : control     (c)
    , input_index (i)
    , use_count   (u1)
    , wg_offset   (u2)
    , num_wg_elems(numIuap)
{}

Opm::RestartIO::RstUDQActive::
RstUDQActive(const std::vector<int>& iuad_arg,
             const std::vector<int>& iuap,
             const std::vector<int>& igph)
    : wg_index { iuap }
{
    using Ix = Opm::RestartIO::Helpers::VectorItems::IUad::index;

    this->iuad.reserve(iuad_arg.size() / UDQDims::entriesPerIUAD());

    for (auto offset = 0*iuad_arg.size();
         offset < iuad_arg.size();
         offset += UDQDims::entriesPerIUAD())
    {
        const auto* uda = &iuad_arg[offset];

        this->iuad.emplace_back(UDQ::udaControl(uda[Ix::UDACode]),
                                uda[Ix::UDQIndex] - 1,
                                uda[Ix::NumIuapElm],
                                uda[Ix::UseCount],
                                uda[Ix::Offset] - 1);
    }

    std::transform(this->wg_index.begin(), this->wg_index.end(),
                   this->wg_index.begin(),
                   [](const int wgIdx) { return wgIdx - 1; });

    this->ig_phase.assign(igph.size(), Phase::OIL);
    std::transform(igph.begin(), igph.end(), this->ig_phase.begin(),
                   [](const int phase)
                   {
                       if (phase == 1) { return Phase::OIL;   }
                       if (phase == 2) { return Phase::WATER; }
                       if (phase == 3) { return Phase::GAS;   }

                       return Phase::OIL;
                   });
}
