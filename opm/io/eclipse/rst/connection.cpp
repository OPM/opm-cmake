/*
  Copyright 2020 Equinor ASA.

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

#include <opm/io/eclipse/rst/connection.hpp>

#include <opm/output/eclipse/VectorItems/connection.hpp>

#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>

namespace VI = ::Opm::RestartIO::Helpers::VectorItems;

namespace {

template <typename T>
T from_int(int);

template<> Opm::Connection::State from_int(int int_state)
{
    return (int_state == 1)
        ? Opm::Connection::State::OPEN
        : Opm::Connection::State::SHUT;
}

template<> Opm::Connection::Direction from_int(int int_dir)
{
    switch (int_dir) {
    case 1: return Opm::Connection::Direction::X;
    case 2: return Opm::Connection::Direction::Y;
    case 3: return Opm::Connection::Direction::Z;
    }

    throw std::invalid_argument {
        "Unable to convert direction value: " +
        std::to_string(int_dir) + " to Direction category"
    };
}

// Note: CTFKind originates in SCON and is indeed a float.
Opm::Connection::CTFKind from_float(float float_kind)
{
    return (float_kind == 0.0f)
        ? Opm::Connection::CTFKind::Defaulted
        : Opm::Connection::CTFKind::DeckValue;
}

float as_float(const double x)
{
    return static_cast<float>(x);
}

float staticDFactorCorrCoeff(const Opm::UnitSystem& usys,
                             const float            coeff)
{
    using M = ::Opm::UnitSystem::measure;

    // Coefficient's units are [D] * [viscosity]
    return as_float(usys.to_si(M::viscosity, usys.to_si(M::dfactor, coeff)));
}

double pressEquivRadius(const float denom,
                        const float skin,
                        const float rw)
{
    // Recall: denom = log(r0 / rw) + skin
    return rw * std::exp(denom - skin);
}

double pressEquivRadius(const Opm::UnitSystem& usys,
                        const std::size_t      nsconz,
                        const float*           scon,
                        const float            denom,
                        const float            skin,
                        const float            rw)
{
    return (nsconz > VI::SConn::PressEquivRad)
        ? usys.to_si(Opm::UnitSystem::measure::length, scon[VI::SConn::PressEquivRad])
        : pressEquivRadius(denom, skin, rw);
}

} // Anonymous namespace

double Opm::RestartIO::RstConnection::inverse_peaceman(double cf, double kh, double rw, double skin)
{
    auto alpha = 3.14159265 * 2 * kh / cf - skin;
    return rw * std::exp(alpha);
}

using M = ::Opm::UnitSystem::measure;

Opm::RestartIO::RstConnection::RstConnection(const UnitSystem& unit_system,
                                             const std::size_t rst_index_,
                                             const int         nsconz,
                                             const int*        icon,
                                             const float*      scon,
                                             const double*     xcon)
    : rst_index { rst_index_ }
      // -----------------------------------------------------------------
      // Integer values (ICON)
    , ijk             { icon[VI::IConn::CellI] - 1 ,
                        icon[VI::IConn::CellJ] - 1 ,
                        icon[VI::IConn::CellK] - 1 }
    , state           { from_int<Connection::State>(icon[VI::IConn::ConnStat]) }
    , drain_sat_table { icon[VI::IConn::Drainage] }
    , imb_sat_table   { icon[VI::IConn::Imbibition] }
    , completion      { icon[VI::IConn::ComplNum] }
    , dir             { from_int<Connection::Direction>(icon[VI::IConn::ConnDir]) }
    , segment         { icon[VI::IConn::Segment] }
      // -----------------------------------------------------------------
      // Float values (SCON)
    , cf_kind                { from_float(1.0f) }
    , skin_factor            { scon[VI::SConn::SkinFactor] }
    , cf                     { as_float(unit_system.to_si(M::transmissibility, scon[VI::SConn::ConnTrans])) }
    , depth                  { as_float(unit_system.to_si(M::length, scon[VI::SConn::Depth])) }
    , diameter               { as_float(unit_system.to_si(M::length, scon[VI::SConn::Diameter])) }
    , kh                     { as_float(unit_system.to_si(M::effective_Kh, scon[VI::SConn::EffectiveKH])) }
    , denom                  { scon[VI::SConn::CFDenom] }
    , length                 { as_float(unit_system.to_si(M::length, scon[VI::SConn::EffectiveLength])) }
    , static_dfac_corr_coeff { staticDFactorCorrCoeff(unit_system, scon[VI::SConn::StaticDFacCorrCoeff]) }
    , segdist_end            { as_float(unit_system.to_si(M::length, scon[VI::SConn::SegDistEnd])) }
    , segdist_start          { as_float(unit_system.to_si(M::length, scon[VI::SConn::SegDistStart])) }
      // -----------------------------------------------------------------
      // Double values (XCON)
    , oil_rate   { unit_system.to_si(M::liquid_surface_rate, xcon[VI::XConn::OilRate]) }
    , water_rate { unit_system.to_si(M::liquid_surface_rate, xcon[VI::XConn::WaterRate]) }
    , gas_rate   { unit_system.to_si(M::gas_surface_rate, xcon[VI::XConn::GasRate]) }
    , pressure   { unit_system.to_si(M::pressure, xcon[VI::XConn::Pressure]) }
    , resv_rate  { unit_system.to_si(M::rate, xcon[VI::XConn::ResVRate]) }
      // -----------------------------------------------------------------
      // Derived quantities
    , r0 { pressEquivRadius(unit_system, nsconz, scon, this->denom, this->skin_factor, this->diameter / 2) }
{
    if (static_cast<std::size_t>(nsconz) > VI::SConn::CFInDeck) {
        this->cf_kind = from_float(scon[VI::SConn::CFInDeck]);
    }
}
