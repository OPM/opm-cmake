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

#ifndef OPM_WRITE_RESTART_HELPERS_HPP
#define OPM_WRITE_RESTART_HELPERS_HPP

#include <vector>

// Forward declarations

namespace Opm {

    class EclipseGrid;
    class EclipseState;
    class Schedule;
    class Well;

} // Opm

namespace Opm { namespace RestartIO { namespace Helpers {

    std::vector<double>
    createDoubHead(const EclipseState& es,
                   const Schedule&     sched,
                   const std::size_t   lookup_step,
                   const double        simTime);



    std::vector<int>
    createInteHead(const EclipseState& es,
                   const EclipseGrid&  grid,
                   const Schedule&     sched,
                   const double        simTime,
                   const int           num_solver_steps,
                   const int           lookup_step,   // The integer index used to look up dynamic properties, e.g. the number of well.
                   const int           report_step);  // The integer number this INTEHEAD keyword will be saved to, typically report_step = lookup_step + 1.

    std::vector<bool>
    createLogiHead(const EclipseState& es);

    std::vector<int> serialize_ICON(int report_step,
                                    int ncwmax,
                                    int niconz, // should be entry 32 from createInteHead
                                    const std::vector<const Well*>& sched_wells);

}}} // Opm::RestartIO::Helpers

#endif  // OPM_WRITE_RESTART_HELPERS_HPP
