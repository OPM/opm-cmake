/*
  Copyright (c) 2024 Equinor ASA

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

#ifndef OPM_OUTPUT_ECLIPSE_VECTOR_UDQ_HPP
#define OPM_OUTPUT_ECLIPSE_VECTOR_UDQ_HPP

#include <vector>

namespace Opm::RestartIO::Helpers::VectorItems {

    namespace IUad {
        enum index : std::vector<int>::size_type {
            UDACode    = 0, // Integer code for keyword/item combination.
            UDQIndex   = 1, // UDQ insertion index (one-based)

            NumIuapElm = 2, // Number of elements in IUAP for this UDA.  One
                            // element for a regular group or well level
                            // UDA, two elements for a FIELD level UDA.

            UseCount   = 3, // Number of times this UDA is used in this
                            // way--i.e., number of wells (or groups) for
                            // which this UDA supplies this particular limit
                            // in this particular keyword.

            Offset     = 4, // One-based start offset in IUAP for this UDA.
        };

        namespace Value {
            enum IuapElems : int {
                Regular = 1, // UDA applies to a well or a non-field group.
                Field   = 2, // UDA applies to the field group.
            };
        } // namespace Value
    } // namespace IUad

    namespace ZUdn {
        enum index : std::vector<const char*>::size_type {
            Keyword = 0, // UDQ name
            Unit    = 1, // UDQ Unit
        };
    } // namespace ZUdn

} // Opm::RestartIO::Helpers::VectorItems

#endif // OPM_OUTPUT_ECLIPSE_VECTOR_UDQ_HPP
