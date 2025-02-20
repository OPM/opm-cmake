/*
  Copyright 2023 SINTEF Digital

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
#ifndef OPM_TIMINGMACROS_HPP
#define OPM_TIMINGMACROS_HPP

// This file defines macros
// OPM_TIMEBLOCK - time block of main part of codes which do not effect performance
// OPM_TIMEFUNCTION - time block of main part of codes which do not effect performance with name from function
// OPM_TIMEBLOCK_LOCAL - detailed timing which may effect performance
// OPM_TIMEFUNCTION_LOCAL - detailed timing which may effect performance with name from function
#if USE_TRACY
#if USE_TRACY_LOCAL
#define DETAILED_PROFILING 1
#endif
#endif

#ifndef DETAILED_PROFILING
#define DETAILED_PROFILING 0 // set to 1 to enable invasive profiling
#endif


#if USE_TRACY
#define TRACY_ENABLE 1
#include <tracy/Tracy.hpp>
#define OPM_TIMEBLOCK(blockname) ZoneNamedN(blockname, #blockname, true)
#define OPM_TIMEFUNCTION() ZoneNamedN(myname, __func__, true)
#if DETAILED_PROFILING
#define OPM_TIMEBLOCK_LOCAL(blockname) ZoneNamedN(blockname, #blockname, true)
#define OPM_TIMEFUNCTION_LOCAL() ZoneNamedN(myname, __func__, true)
#endif
#endif

#ifndef OPM_TIMEBLOCK
#define OPM_TIMEBLOCK(x)\
    do { /* nothing */ } while (false)
#endif

// detailed timing which may effect performance
#ifndef OPM_TIMEBLOCK_LOCAL
#define OPM_TIMEBLOCK_LOCAL(x)\
    do { /* nothing */ } while (false)
#endif

#ifndef OPM_TIMEFUNCTION
#define OPM_TIMEFUNCTION()\
    do { /* nothing */ } while (false)
#endif

#ifndef OPM_TIMEFUNCTION_LOCAL
#define OPM_TIMEFUNCTION_LOCAL()\
    do { /* nothing */ } while (false)
#endif

#endif // OPM_TIMINGMACROS_HPP
