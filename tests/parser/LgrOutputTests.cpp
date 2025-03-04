/*
  Copyright (C) 2023 Equinor

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

#define BOOST_TEST_MODULE LgrOutputTests

#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>

#include <opm/io/eclipse/EGrid.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Grid/Carfin.hpp>
#include <opm/input/eclipse/EclipseState/Grid/LgrCollection.hpp>

#include <opm/input/eclipse/Deck/DeckSection.hpp>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>

#include <array>
#include <cstddef>
#include <filesystem>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using namespace Opm;

namespace {

LgrCollection read_lgr(const std::string& deck_string,
                       const std::size_t nx,
                       const std::size_t ny,
                       const std::size_t nz)
{
    Opm::Parser parser;
    Opm::EclipseGrid eclipse_grid(GridDims(nx,ny,nz));    
    Opm::Deck deck = parser.parseString(deck_string);
    const Opm::GRIDSection gridSection ( deck );
    return LgrCollection(gridSection, eclipse_grid);
}

std::pair<std::vector<double>, std::vector<double>>
read_cpg_from_egrid(const std::string& file_path,
                    const std::string& lgr_label)
{
    Opm::EclIO::EGrid egrid_global(file_path, lgr_label);
    egrid_global.load_grid_data();

    const auto& global_coord = egrid_global.get_coord();
    const auto& global_zcorn = egrid_global.get_zcorn();

    return {
        std::piecewise_construct,
        std::forward_as_tuple(global_coord.begin(), global_coord.end()),
        std::forward_as_tuple(global_zcorn.begin(), global_zcorn.end())
    };
}

} // Anonymous namespace

BOOST_AUTO_TEST_CASE(TestLgrOutputBasicLGR) { 
    const std::string deck_string = R"(
RUNSPEC

DIMENS
  3 3 1 /

GRID

CARFIN
-- NAME I1-I2 J1-J2 K1-K2 NX NY NZ
'LGR1'  2  2  2  2  1  1  3  3   /
ENDFIN

DX 
  9*1000 /
DY
	9*1000 /
DZ
	9*20 /

TOPS
	9*8325 /

PORO
  9*0.15 /

PERMX
  9*1 /

COPY
  PERMX PERMZ /
  PERMX PERMY /
/

EDIT

OIL
GAS

TITLE
Opm-the-best-open-source-simulator

START
16 JUN 1988 /

PROPS

REGIONS

SOLUTION

SCHEDULE
)";

    Opm::UnitSystem units(1);
    std::vector<Opm::NNCdata> vecNNC;
    std::array<int,3> global_grid_dim = {3,3,1};
    std::vector<double> coord_g, zcorn_g, coord_l, zcorn_l, coord_g_opm, zcorn_g_opm, coord_l_opm, zcorn_l_opm;
    // Intialize LgrCollection from string.
    LgrCollection lgr_col = read_lgr(deck_string,global_grid_dim[0],global_grid_dim[1],global_grid_dim[2]);
    //  Read global COORD and ZCORN from reference simulator output.
    std::tie(coord_g, zcorn_g) = read_cpg_from_egrid("CARFIN5.EGRID", "global");
    //  Read LGR CELL COORD and ZCORN from reference simulator output.
    std::tie(coord_l, zcorn_l) = read_cpg_from_egrid("CARFIN5.EGRID", "LGR1");
    //  Eclipse Grid is intialzied with COORD and ZCORN.
    Opm::EclipseGrid eclipse_grid_file(global_grid_dim, coord_g, zcorn_g);    
    //  LgrCollection is used to initalize LGR Cells in the Eclipse Grid.
    eclipse_grid_file.init_lgr_cells(lgr_col);
    // LGR COORD and ZCORN is parsed to EclipseGridLGR children cell. (Simulates the process of recieving the LGR refinement.)   
    eclipse_grid_file.set_lgr_refinement("LGR1",coord_l,zcorn_l);
    // Intialize host_cell numbering.
    eclipse_grid_file.init_children_host_cells();
    // Save EclipseGrid.
    eclipse_grid_file.save("OPMCARFIN5.EGRID",false,vecNNC,units);
    // Once the new EGRID is saved, another EclipseGrid Object is created for the sake of comparison.
    std::tie(coord_g_opm, zcorn_g_opm)  = read_cpg_from_egrid("OPMCARFIN5.EGRID", "global");
    //  Read LGR CELL COORD and ZCORN from reference simulator output.
    std::tie(coord_l_opm, zcorn_l_opm) = read_cpg_from_egrid("OPMCARFIN5.EGRID", "LGR1");
    //  Eclipse Grid is intialzied with COORD and ZCORN. 
    Opm::EclipseGrid eclipse_grid_OPM(global_grid_dim, coord_g_opm, zcorn_g_opm);    
//  LgrCollection is used to initalize LGR Cells in the Eclipse Grid.
    eclipse_grid_OPM.init_lgr_cells(lgr_col);
    // LGR COORD and ZCORN is parsed to EclipseGridLGR children cell. (Simulates the process of recieving the LGR refinement.)   
    eclipse_grid_OPM.set_lgr_refinement("LGR1",coord_l_opm,zcorn_l_opm);
    // Intialize host_cell numbering.
    eclipse_grid_OPM.init_children_host_cells();
    BOOST_CHECK_EQUAL( coord_g_opm.size() , coord_g.size());
    BOOST_CHECK_EQUAL( zcorn_g_opm.size() , zcorn_g.size());
    BOOST_CHECK_EQUAL( coord_l_opm.size() , coord_l.size());
    BOOST_CHECK_EQUAL( zcorn_l_opm.size() , zcorn_l.size());
    std::size_t index ;
    for (index = 0; index < coord_g.size(); index++) {
      BOOST_CHECK_EQUAL( coord_g_opm[index] , coord_g[index]);
    }
    for (index = 0; index < zcorn_g.size(); index++) {
      BOOST_CHECK_EQUAL( zcorn_g_opm[index] , zcorn_g[index]);
    }
  }

BOOST_AUTO_TEST_CASE(TestLgrOutputColumnLGR) { 
    const std::string deck_string = R"(
RUNSPEC

DIMENS
  3 3 1 /

GRID

CARFIN
-- NAME I1-I2 J1-J2 K1-K2 NX NY NZ
'LGR1'  1  1  1  2  1  1  2  4   1/
ENDFIN

DX 
  9*1000 /
DY
	9*1000 /
DZ
	9*20 /

TOPS
	9*8325 /

PORO
  9*0.15 /

PERMX
  9*1 /

COPY
  PERMX PERMZ /
  PERMX PERMY /
/

EDIT

OIL
GAS

TITLE
The title

START
16 JUN 1988 /

PROPS

REGIONS

SOLUTION

SCHEDULE
)";

    Opm::UnitSystem units(1);
    std::vector<Opm::NNCdata> vecNNC;
    std::array<int,3> global_grid_dim = {3,3,1};
    std::vector<double> coord_g, zcorn_g, coord_l, zcorn_l, coord_g_opm, zcorn_g_opm, coord_l_opm, zcorn_l_opm;
    // Intialize LgrCollection from string.
    LgrCollection lgr_col = read_lgr(deck_string,global_grid_dim[0],global_grid_dim[1],global_grid_dim[2]);
    //  Read global COORD and ZCORN from reference simulator output.
    std::tie(coord_g, zcorn_g) = read_cpg_from_egrid("CARFIN-COLUMN.EGRID", "global");
    //  Read LGR CELL COORD and ZCORN from reference simulator output.
    std::tie(coord_l, zcorn_l) = read_cpg_from_egrid("CARFIN-COLUMN.EGRID", "LGR1");
    //  Eclipse Grid is intialzied with COORD and ZCORN.
    Opm::EclipseGrid eclipse_grid_file(global_grid_dim, coord_g, zcorn_g);    
    //  LgrCollection is used to initalize LGR Cells in the Eclipse Grid.
    eclipse_grid_file.init_lgr_cells(lgr_col);
    // LGR COORD and ZCORN is parsed to EclipseGridLGR children cell. (Simulates the process of recieving the LGR refinement.)   
    eclipse_grid_file.set_lgr_refinement("LGR1",coord_l,zcorn_l);
    // Intialize host_cell numbering.
    eclipse_grid_file.init_children_host_cells();
    // Save EclipseGrid.
    eclipse_grid_file.save("OPMCARFIN-COLUMN.EGRID",false,vecNNC,units);
    // Once the new EGRID is saved, another EclipseGrid Object is created for the sake of comparison.
    std::tie(coord_g_opm, zcorn_g_opm)  = read_cpg_from_egrid("OPMCARFIN-COLUMN.EGRID", "global");
    //  Read LGR CELL COORD and ZCORN from reference simulator output.
    std::tie(coord_l_opm, zcorn_l_opm) = read_cpg_from_egrid("OPMCARFIN-COLUMN.EGRID", "LGR1");
    //  Eclipse Grid is intialzied with COORD and ZCORN. 
    Opm::EclipseGrid eclipse_grid_OPM(global_grid_dim, coord_g_opm, zcorn_g_opm);    
//  LgrCollection is used to initalize LGR Cells in the Eclipse Grid.
    eclipse_grid_OPM.init_lgr_cells(lgr_col);
    // LGR COORD and ZCORN is parsed to EclipseGridLGR children cell. (Simulates the process of recieving the LGR refinement.)   
    eclipse_grid_OPM.getLGRCell(0).set_lgr_refinement(coord_l_opm,zcorn_l_opm);
    // Intialize host_cell numbering.
    eclipse_grid_OPM.init_children_host_cells();
    BOOST_CHECK_EQUAL( coord_g_opm.size() , coord_g.size());
    BOOST_CHECK_EQUAL( zcorn_g_opm.size() , zcorn_g.size());
    BOOST_CHECK_EQUAL( coord_l_opm.size() , coord_l.size());
    BOOST_CHECK_EQUAL( zcorn_l_opm.size() , zcorn_l.size());
    std::size_t index ;
    for (index = 0; index < coord_g.size(); index++) {
      BOOST_CHECK_EQUAL( coord_g_opm[index] , coord_g[index]);
    }
    for (index = 0; index < zcorn_g.size(); index++) {
      BOOST_CHECK_EQUAL( zcorn_g_opm[index] , zcorn_g[index]);
    }
  }

BOOST_AUTO_TEST_CASE(TestLgrOutputDoubleLGR) { 
    const std::string deck_string = R"(
RUNSPEC

DIMENS
  3 3 1 /

GRID

CARFIN
-- NAME I1-I2 J1-J2 K1-K2 NX NY NZ
'LGR1'  2  2  2  2  1  1  2  2  1  /
ENDFIN

CARFIN
-- NAME I1-I2 J1-J2 K1-K2 NX NY NZ
'LGR2'  2  2  1  1  1  1  2  2  1 /
ENDFIN

DX 
  9*1000 /
DY
	9*1000 /
DZ
	9*20 /

TOPS
	9*8325 /

PORO
  9*0.15 /

PERMX
  9*1 /

COPY
  PERMX PERMZ /
  PERMX PERMY /
/

EDIT

OIL
GAS

TITLE
The title

START
16 JUN 1988 /

PROPS

REGIONS

SOLUTION

SCHEDULE
)";
    Opm::UnitSystem units(1);
    std::vector<Opm::NNCdata> vecNNC;
    std::array<int,3> global_grid_dim = {3,3,1};
    std::vector<double> coord_g, zcorn_g, coord_l1, zcorn_l1, coord_l2, zcorn_l2,
                        coord_g_opm, zcorn_g_opm, coord_l1_opm, zcorn_l1_opm, coord_l2_opm, zcorn_l2_opm;
    // Intialize LgrCollection from string.
    LgrCollection lgr_col = read_lgr(deck_string,global_grid_dim[0],global_grid_dim[1],global_grid_dim[2]);
    //  Read global COORD and ZCORN from reference simulator output.
    std::tie(coord_g, zcorn_g) = read_cpg_from_egrid("CARFIN-DOUBLE.EGRID", "global");
    //  Read LGR CELL COORD and ZCORN from reference simulator output.
    std::tie(coord_l1, zcorn_l1) = read_cpg_from_egrid("CARFIN-DOUBLE.EGRID", "LGR1");
    std::tie(coord_l2, zcorn_l2) = read_cpg_from_egrid("CARFIN-DOUBLE.EGRID", "LGR2");

    //  Eclipse Grid is intialzied with COORD and ZCORN.
    Opm::EclipseGrid eclipse_grid_file(global_grid_dim, coord_g, zcorn_g);    
    //  LgrCollection is used to initalize LGR Cells in the Eclipse Grid.
    eclipse_grid_file.init_lgr_cells(lgr_col);
    // LGR COORD and ZCORN is parsed to EclipseGridLGR children cell. (Simulates the process of recieving the LGR refinement.)   
    eclipse_grid_file.getLGRCell(1).set_lgr_refinement(coord_l1,zcorn_l1);
    eclipse_grid_file.getLGRCell(0).set_lgr_refinement(coord_l2,zcorn_l2);
    // Intialize host_cell numbering.
    eclipse_grid_file.init_children_host_cells();
    // Save EclipseGrid.
    eclipse_grid_file.save("OPMCARFIN-DOUBLE.EGRID",false,vecNNC,units);
    // Once the new EGRID is saved, another EclipseGrid Object is created for the sake of comparison.
    std::tie(coord_g_opm, zcorn_g_opm)  = read_cpg_from_egrid("OPMCARFIN-DOUBLE.EGRID", "global");
    //  Read LGR CELL COORD and ZCORN from reference simulator output.
    std::tie(coord_l1_opm, zcorn_l1_opm) = read_cpg_from_egrid("OPMCARFIN-DOUBLE.EGRID", "LGR1");
    std::tie(coord_l2_opm, zcorn_l2_opm) = read_cpg_from_egrid("OPMCARFIN-DOUBLE.EGRID", "LGR2");
    //  Eclipse Grid is intialzied with COORD and ZCORN. 
    Opm::EclipseGrid eclipse_grid_OPM(global_grid_dim, coord_g_opm, zcorn_g_opm);    
//  LgrCollection is used to initalize LGR Cells in the Eclipse Grid.
    eclipse_grid_OPM.init_lgr_cells(lgr_col);
    // LGR COORD and ZCORN is parsed to EclipseGridLGR children cell. (Simulates the process of recieving the LGR refinement.)   
    eclipse_grid_OPM.getLGRCell(1).set_lgr_refinement(coord_l1_opm,zcorn_l1_opm);
    eclipse_grid_OPM.getLGRCell(0).set_lgr_refinement(coord_l2_opm,zcorn_l2_opm);

    // Intialize host_cell numbering.
    eclipse_grid_OPM.init_children_host_cells();
    BOOST_CHECK_EQUAL( coord_g_opm.size() , coord_g.size());
    BOOST_CHECK_EQUAL( zcorn_g_opm.size() , zcorn_g.size());
    BOOST_CHECK_EQUAL( coord_l1_opm.size() , coord_l1.size());
    BOOST_CHECK_EQUAL( zcorn_l1_opm.size() , zcorn_l1.size());
    BOOST_CHECK_EQUAL( coord_l2_opm.size() , coord_l2.size());
    BOOST_CHECK_EQUAL( zcorn_l2_opm.size() , zcorn_l2.size());    
    std::size_t index ;
    for (index = 0; index < coord_g.size(); index++) {
      BOOST_CHECK_EQUAL( coord_g_opm[index] , coord_g[index]);
    }
    for (index = 0; index < zcorn_g.size(); index++) {
      BOOST_CHECK_EQUAL( zcorn_g_opm[index] , zcorn_g[index]);
    }
  }

BOOST_AUTO_TEST_CASE(TestLgrOutputNESTED) { 
    const std::string deck_string = R"(
RUNSPEC

DIMENS
  3 3 1 /

GRID

CARFIN
-- NAME I1-I2 J1-J2 K1-K2 NX NY NZ
LGR1  2  2  2  2  1  1  3  3   1 1*  GLOBAL/
ENDFIN

CARFIN
-- NAME I1-I2 J1-J2 K1-K2 NX NY NZ
LGR2  2  2  2  2  1  1  3  3   1 1*  LGR1/
ENDFIN

DX 
  9*1000 /
DY
	9*1000 /
DZ
	9*20 /

TOPS
	9*8325 /

PORO
  9*0.15 /

PERMX
  9*1 /

COPY
  PERMX PERMZ /
  PERMX PERMY /
/

EDIT

OIL
GAS

TITLE
The title

START
16 JUN 1988 /

PROPS

REGIONS

SOLUTION

SCHEDULE
)";
    Opm::UnitSystem units(1);
    std::vector<Opm::NNCdata> vecNNC;
    std::array<int,3> global_grid_dim = {3,3,1};
    std::vector<double> coord_g, zcorn_g, coord_l1, zcorn_l1, coord_l2, zcorn_l2;
    // Intialize LgrCollection from string.
    LgrCollection lgr_col = read_lgr(deck_string,global_grid_dim[0],global_grid_dim[1],global_grid_dim[2]);
    //  Read global COORD and ZCORN from reference simulator output.
    std::tie(coord_g, zcorn_g) = read_cpg_from_egrid("CARFIN-NESTED.EGRID", "global");
    //  Read LGR CELL COORD and ZCORN from reference simulator output.
    std::tie(coord_l1, zcorn_l1) = read_cpg_from_egrid("CARFIN-NESTED.EGRID", "LGR1");
    std::tie(coord_l2, zcorn_l2) = read_cpg_from_egrid("CARFIN-NESTED.EGRID", "LGR2");

    //  Eclipse Grid is intialzied with COORD and ZCORN.
     Opm::EclipseGrid eclipse_grid_file(global_grid_dim, coord_g, zcorn_g);    
    //  LgrCollection is used to initalize LGR Cells in the Eclipse Grid.
    eclipse_grid_file.init_lgr_cells(lgr_col);
    // LGR COORD and ZCORN is parsed to EclipseGridLGR children cell. (Simulates the process of recieving the LGR refinement.)   
    eclipse_grid_file.set_lgr_refinement("LGR1",coord_l1,zcorn_l1);
    eclipse_grid_file.set_lgr_refinement("LGR2",coord_l2,zcorn_l2);
    // Intialize host_cell numbering.
    eclipse_grid_file.init_children_host_cells();
    // Save EclipseGrid.
    eclipse_grid_file.save("OPMCARFIN-NESTED.EGRID",false,vecNNC,units);
  }

  BOOST_AUTO_TEST_CASE(Test_lgr_host_cells_logical) { 
    const std::string deck_string = R"(
RUNSPEC

DIMENS
  3 1 1 /

GRID

CARFIN
-- NAME I1-I2 J1-J2 K1-K2 NX NY NZ
LGR1  1  3  1  1  1  1  12  1   1 1*  GLOBAL/
ENDFIN


DX 
  3*1000 /
DY
	3*1000 /
DZ
	3*20 /

TOPS
	3*8325 /

PORO
  3*0.15 /

PERMX
  3*1 /

COPY
  PERMX PERMZ /
  PERMX PERMY /
/

EDIT

OIL
GAS

TITLE
The title

START
16 JUN 1988 /

PROPS

REGIONS

SOLUTION

SCHEDULE
)";
    
    Opm::Parser parser;
    Opm::Deck deck = parser.parseString(deck_string);
    Opm::EclipseState state(deck);
    Opm::EclipseGrid eclipse_grid = state.getInputGrid();
    Opm::LgrCollection lgrs = state.getLgrs();
    // Intialize host_cell numbering.
    eclipse_grid.init_children_host_cells();
    std::vector<int> host_vec_sol = { 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3};
    std::vector<int> host_vec = eclipse_grid.getLGRCell(0).save_hostnum();
    BOOST_CHECK_EQUAL_COLLECTIONS(host_vec.begin(), host_vec.end(), host_vec_sol.begin(), host_vec_sol.end());  
 
  }
