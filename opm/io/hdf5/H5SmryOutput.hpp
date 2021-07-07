/*
   Copyright 2019 Statoil ASA.

   This file is part of the Open Porous Media project (OPM).

   OPM is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   OPM is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with OPM.  If not, see <http://www.gnu.org/licenses/>.
   */

#ifndef OPM_IO_H5SMRYOUTPUT_HPP
#define OPM_IO_H5SMRYOUTPUT_HPP

#include <string>

#include <opm/io/hdf5/Hdf5Util.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>


namespace Opm {

class EclipseState;

}

namespace Opm { namespace Hdf5IO {

class H5SmryOutput
{
public:
    H5SmryOutput(const std::vector<std::string>& valueKeys, const std::vector<std::string>& valueUnits,
                 const EclipseState& es, const time_t start_time);

    void write(const std::vector<float>& ts_data, int report_step);

    ~H5SmryOutput();


private:

    hid_t m_file_id;
    int m_nTimeSteps;
    int m_maxTimeSteps;

    int m_nVect;

    std::array<int, 3> ijk_from_global_index(const GridDims& dims, int globInd) const;
    std::vector<std::string> make_modified_keys(const std::vector<std::string> valueKeys, const GridDims& dims);
};


}} // namespace Opm::EclIO

#endif // OPM_IO_H5SMRYOUTPUT_HPP
