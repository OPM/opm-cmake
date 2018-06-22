/*
  Copyright (c) 2018 Statoil ASA

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

#include <opm/output/eclipse/WriteRestartHelpers.hpp>
#include <ert/ecl_well/well_const.h> // containts ICON_XXX_INDEX
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/Units/UnitSystem.hpp>
#include <vector>

// ----------------------------------------------------------------------------
std::vector<double>
Opm::RestartIO::Helpers::
serialize_SCON(int lookup_step,
               int ncwmax,
               int nsconz,
               const std::vector<const Well*>& sched_wells,
               const UnitSystem& units)
// ----------------------------------------------------------------------------
{
    const size_t well_field_size = ncwmax * nsconz;
    std::vector<double> data(sched_wells.size() * well_field_size, 0);
    size_t well_offset = 0;
    for (const Opm::Well* well : sched_wells) {
        const auto& connections = well->getConnections( lookup_step );
        size_t connection_offset = 0;
        bool explicit_ctf_not_found = false;
        for (const auto& connection : connections) {
            const size_t offset = well_offset + connection_offset;
            const auto& ctf = connection.getConnectionTransmissibilityFactorAsValueObject();
            if (ctf.hasValue()) {
                // CTF explicitly set in deck, overrides calculation
                // from Peaceman model.  We should also give the Kh
                // factor, we output an explicitly invalid value
                // instead. This is acceptable since it will not be
                // used (the explicit CTF factor is used instead).
                const double ctf_SI = ctf.getValue();
                const double ctf_output = units.from_si(UnitSystem::measure::transmissibility, ctf_SI);
                data[ offset + SCON_CF_INDEX ] = ctf_output;
                data[ offset + SCON_KH_INDEX ] = UNIMPLEMENTED_VALUE;
            } else {
                // CTF not set in deck, Peaceman formula used to
                // compute it.  Here we should store the data for the
                // connection required to recalculate the CTF (the Kh
                // factor), as well as the actual CTF used by the
                // simulator, but that requires access to more data
                // from the simulator. As an interim measure we write
                // invalid values and give a warning.
                data[ offset + SCON_CF_INDEX ] = UNIMPLEMENTED_VALUE;
                data[ offset + SCON_KH_INDEX ] = UNIMPLEMENTED_VALUE;
                explicit_ctf_not_found = true;
            }
            connection_offset += nsconz;
        }
        if (explicit_ctf_not_found) {
            OpmLog::warning("restart output connection data missing",
                            "Explicit connection transmissibility factors for well " + well->name() + " missing, writing dummy values to restart file.");
        }
        well_offset += well_field_size;
    }
    return data; 
}

// ----------------------------------------------------------------------------
std::vector<int>
Opm::RestartIO::Helpers::
serialize_ICON(int lookup_step,
               int ncwmax,
               int niconz,
               const std::vector<const Opm::Well*>& sched_wells)
// ----------------------------------------------------------------------------  
{
    const size_t well_field_size = ncwmax * niconz;
    std::vector<int> data(sched_wells.size() * well_field_size, 0);
    size_t well_offset = 0;
    for (const Opm::Well* well : sched_wells) {
        const auto& connections = well->getConnections( lookup_step );
        size_t connection_offset = 0;
        for (const auto& connection : connections) {
            const size_t offset = well_offset + connection_offset;

            data[ offset + ICON_IC_INDEX ] = connection.complnum;
            data[ offset + ICON_I_INDEX ] = connection.getI() + 1;
            data[ offset + ICON_J_INDEX ] = connection.getJ() + 1;
            data[ offset + ICON_K_INDEX ] = connection.getK() + 1;
            data[ offset + ICON_DIRECTION_INDEX ] = connection.getDirection();
            data[ offset + ICON_STATUS_INDEX ] =
                (connection.state == WellCompletion::StateEnum::OPEN) ?
                1 : -1000;
            data[ offset + ICON_SEGMENT_INDEX ] =
                connection.attachedToSegment() ?
                connection.getSegmentNumber() : 0;
            connection_offset += niconz;
        }

        well_offset += well_field_size;
    }

    return data;
}
