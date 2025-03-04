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

#include <opm/input/eclipse/Schedule/ScheduleGrid.hpp>

#include <opm/input/eclipse/Schedule/CompletedCells.hpp>

#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FieldPropsManager.hpp>

#include <cstddef>
#include <string>
#include <fmt/format.h>

Opm::ScheduleGrid::ScheduleGrid(const Opm::EclipseGrid& ecl_grid,
                                const Opm::FieldPropsManager& fpm,
                                Opm::CompletedCells& completed_cells,
                                std::vector<CompletedCells>& completed_cells_lgr)
    : grid  { &ecl_grid }
    , fp    { &fpm }
    , cells { completed_cells}
    , cells_lgr{ completed_cells_lgr}
{
    init_lgr_grid();
}

Opm::ScheduleGrid::ScheduleGrid(Opm::CompletedCells& completed_cells, std::vector<CompletedCells>& completed_cells_lgr)
    : cells { completed_cells}
    , cells_lgr{completed_cells_lgr}
{
    init_lgr_grid();
}

namespace {
    double try_get_value(const Opm::FieldPropsManager& fp,
                         const std::string& kw,
                         const std::size_t active_index)
    {
        if (! fp.has_double(kw)) {
            throw std::logic_error(fmt::format("FieldPropsManager is missing keyword '{}'", kw));
        }

        return fp.try_get<double>(kw)->at(active_index);
    }

    double try_get_ntg_value(const Opm::FieldPropsManager& fp,
                             const std::string& kw,
                             const std::size_t active_index)
    {
        return fp.has_double(kw)
            ? fp.try_get<double>(kw)->at(active_index)
            : 1.0;
    }
}

const Opm::CompletedCells::Cell&
Opm::ScheduleGrid::get_cell(std::size_t i, std::size_t j, std::size_t k) const
{
    if (this->grid == nullptr) {
        return this->cells.get(i, j, k);
    }

    auto [valid, cellRef] = this->cells.try_get(i, j, k);

    if (!valid) {
        cellRef.depth = this->grid->getCellDepth(i, j, k);
        cellRef.dimensions = this->grid->getCellDimensions(i, j, k);

        if (this->grid->cellActive(i, j, k)) {
            const auto active_index = this->grid->getActiveIndex(i, j, k);
            const double porv = try_get_value(*this->fp, "PORV", active_index);
            if (this->grid->cellActiveAfterMINPV(i, j, k, porv)) {
                auto& props = cellRef.props.emplace(CompletedCells::Cell::Props{});
                props.active_index = active_index;
                props.permx = try_get_value(*this->fp, "PERMX", props.active_index);
                props.permy = try_get_value(*this->fp, "PERMY", props.active_index);
                props.permz = try_get_value(*this->fp, "PERMZ", props.active_index);
                props.poro = try_get_value(*this->fp, "PORO", props.active_index);
                props.satnum = this->fp->get_int("SATNUM").at(props.active_index);
                props.pvtnum = this->fp->get_int("PVTNUM").at(props.active_index);
                props.ntg = try_get_ntg_value(*this->fp, "NTG", props.active_index);
            }
        }
    }

    return cellRef;
}

const Opm::CompletedCells::Cell&
Opm::ScheduleGrid::get_cell(std::size_t i, std::size_t j, std::size_t k, std::optional<std::string> tag) const
{
    if (!tag.has_value()) 
    {
        return get_cell(i, j, k);
    }
    else
    {
        return get_cell_lgr(i, j, k,tag.value());
    }
}

const Opm::CompletedCells::Cell&
Opm::ScheduleGrid::get_cell_lgr(std::size_t i, std::size_t j, std::size_t k, std::string tag) const
{
    std::size_t tag_index = label_to_index.at(tag);
    if (tag_index == 0)
    {
        return get_cell(i, j, k);
    }   
    tag_index--;
    const EclipseGridLGR& lgr_grid =  this->grid->getLGRCell(tag);
    auto [valid, cellRef] = this->cells_lgr[tag_index].try_get(i, j, k);
    int father_global_id =  this->grid->getLGR_global_father(cellRef.global_index,tag);
    auto [fi, fj, fk] = this->grid->getIJK(father_global_id);
    if (!valid) {
        // this part relies on the ZCORN and COORDS of the host cells that have not been parsed yet.
        // the following implementations are a path that compute depths and dimensions of the LGR cells
        // based on the host cells.
        cellRef.depth = this->grid->getCellDepthLGR(i, j, k, tag);
        cellRef.dimensions = this->grid->getCellDimensionsLGR(fi, fj, fk, tag);
        
        if (this->grid->cellActive(fi, fj, fk) && lgr_grid.cellActive(i, j, k)) 
        {
            const auto father_active_index = this->grid->getActiveIndex(fi, fj, fk);
            const auto active_index = lgr_grid.getActiveIndex(i, j, k);
            const double porv = try_get_value(*this->fp, "PORV", father_active_index);
            if (this->grid->cellActiveAfterMINPV(i, j, k, porv)) {
                auto& props = cellRef.props.emplace(CompletedCells::Cell::Props{});
                props.active_index = active_index;
                props.permx = try_get_value(*this->fp, "PERMX", father_active_index);
                props.permy = try_get_value(*this->fp, "PERMY", father_active_index);
                props.permz = try_get_value(*this->fp, "PERMZ", father_active_index);
                props.poro = try_get_value(*this->fp, "PORO", father_active_index);
                props.satnum = this->fp->get_int("SATNUM").at(father_active_index);
                props.pvtnum = this->fp->get_int("PVTNUM").at(father_active_index);
                props.ntg = try_get_ntg_value(*this->fp, "NTG", father_active_index);
            }
        }
    }

    return cellRef;
 }



const Opm::EclipseGrid* Opm::ScheduleGrid::get_grid() const
{
    return this->grid;
}

void Opm::ScheduleGrid::init_lgr_grid()
{
    std::size_t index = 0;
    for (const std::string& label : this->grid->get_all_labels())
    {
        label_to_index[label] = index;
        index++;
    }
}
