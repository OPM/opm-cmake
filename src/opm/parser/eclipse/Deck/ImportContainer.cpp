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

  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <unordered_set>

#include <fmt/format.h>

#include <opm/common/OpmLog/LogUtil.hpp>
#include <opm/io/eclipse/EclFile.hpp>
#include <opm/parser/eclipse/Deck/ImportContainer.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>

namespace Opm {

ImportContainer::ImportContainer(const Parser& parser, const UnitSystem& unit_system, const std::string& fname, bool formatted, std::size_t deck_size) {
    const std::unordered_set<std::string> float_kw = {"COORD", "MULTPV", "NTG", "PERMX", "PERMY", "PERMZ", "PORO", "SWATINIT", "ZCORN"};
    const std::unordered_set<std::string> int_kw = {"ACTNUM", "EQLNUM", "FIPNUM" , "MULTNUM", "FLUXNUM", "OPERNUM", "SATNUM"};

    EclIO::EclFile ecl_file(fname, EclIO::EclFile::Formatted{formatted});
    const auto& header = ecl_file.getList();
    for (std::size_t kw_index = 0; kw_index < header.size(); kw_index++) {
        const auto& [name, data_type, _] = header[kw_index];
        (void)_;

        if (float_kw.count(name)) {
            const auto& parser_kw = parser.getKeyword(name);
            if (data_type == EclIO::REAL) {
                auto& float_data = ecl_file.get<float>(kw_index);
                std::vector<double> double_data;
                std::copy(float_data.begin(), float_data.end(), std::back_inserter(double_data));
                this->keywords.emplace_back(parser_kw, double_data, unit_system, unit_system);
            } else if (data_type == EclIO::DOUB) {
                auto& double_data = ecl_file.get<double>(kw_index);
                this->keywords.emplace_back(parser_kw, double_data, unit_system, unit_system);
            }
            deck_size += 1;
            auto msg = fmt::format("{:5} Loading {:<8} from IMPORT file {}", deck_size, name, fname);
            OpmLog::info(msg);
            continue;
        }

        if (int_kw.count(name)) {
            const auto& parser_kw = parser.getKeyword(name);
            const auto& data = ecl_file.get<int>(kw_index);
            this->keywords.emplace_back(parser_kw, data);
            deck_size += 1;
            auto msg = fmt::format("{:5} Loading {:<8} from IMPORT file {}", deck_size, name, fname);
            continue;
        }

        OpmLog::info(fmt::format("{:<5} Skipping {} from IMPORT file {}", "", name, fname));
    }
}



}
