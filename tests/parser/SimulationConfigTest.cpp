/*
  Copyright 2015 Statoil ASA.

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

#define BOOST_TEST_MODULE SimulationConfigTests

#include <boost/test/unit_test.hpp>

#include <opm/input/eclipse/EclipseState/SimulationConfig/RockConfig.hpp>

#include <opm/common/utility/OpmInputError.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FieldPropsManager.hpp>
#include <opm/input/eclipse/EclipseState/Runspec.hpp>
#include <opm/input/eclipse/EclipseState/SimulationConfig/SimulationConfig.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>

#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/C.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Deck/DeckSection.hpp>

#include <cstddef>
#include <string>

using namespace Opm;

const std::string& inputStr = "RUNSPEC\n"
                              "EQLOPTS\n"
                              "THPRES /\n"
                              "DIMENS\n"
                              "10 3 4 /\n"
                              "\n"
                              "GRID\n"
                              "REGIONS\n"
                              "EQLNUM\n"
                              "10*1 10*2 100*3 /\n "
                              "\n"
                              "SOLUTION\n"
                              "THPRES\n"
                              "1 2 12.0/\n"
                              "1 3 5.0/\n"
                              "2 3 7.0/\n"
                              "/\n"
                              "\n";


const std::string& inputStr_noTHPRES = "RUNSPEC\n"
                                       "DIMENS\n"
                                       "10 3 4 /\n"
                                       "\n"
                                       "GRID\n"
                                       "REGIONS\n"
                                       "EQLNUM\n"
                                       "10*1 10*2 100*3 /\n "
                                       "\n"
                                       "SOLUTION\n"
                                       "\n";

const std::string& inputStr_cpr = "RUNSPEC\n"
    "DIMENS\n"
    "10 3 4 /\n"
    "CPR\n"
    "/\n"
    "SUMMARY\n";


const std::string& inputStr_INVALID = "RUNSPEC\n"
    "DIMENS\n"
    "10 3 4 /\n"
    "CPR\n"
    "WEll 10 10 17/"
    "/\n"
    "SUMMARY\n";



const std::string& inputStr_cpr_in_SUMMARY = "RUNSPEC\n"
    "DIMENS\n"
    "10 3 4 /\n"
    "SUMMARY\n"
        "CPR\n"
        "well1 10 27 10/\n/\n";

const std::string& inputStr_cpr_BOTH = "RUNSPEC\n"
    "DIMENS\n"
    "10 3 4 /\n"
    "CPR\n"
    "/\n"
    "SUMMARY\n"
    "CPR\n"
    "well1 10 20 30/\n/\n";

const std::string& inputStr_nonnc = "RUNSPEC\n"
    "NONNC\n"
    "DIMENS\n"
    "10 3 4 /\n";

const std::string& inputStr_vap_dis = R"(
RUNSPEC
VAPOIL
DISGAS
VAPWAT
DIMENS
 10 3 4 /
GRID
REGIONS
)";

const std::string& inputStr_bc = R"(
RUNSPEC
DIMENS
 10 3 4 /
GRID
SOLUTION
BC
1 1 1 1 1 10 X- FREE GAS/
/
REGIONS
)";

const std::string& inputStr_bccon = R"(
RUNSPEC
DIMENS
 10 3 4 /
GRID
BCCON
1 1 1 1 1 1 1 X- /
2 10 10 4* X /
/
REGIONS
)";

namespace {
    std::string simDeckStringTEMP()
    {
        return R"(
RUNSPEC

TEMP

DIMENS
  10 3 4 /
)";
    }

    std::string simDeckStringTHERMAL()
    {
        return R"(
RUNSPEC

THERMAL

DIMENS
  10 3 4 /
)";
    }
}

static Deck createDeck(const std::string& input) {
    Opm::Parser parser;
    return parser.parseString(input);
}

BOOST_AUTO_TEST_CASE(SimulationConfigBCCON) {
    auto deck = createDeck(inputStr_bccon);
    BCConfig bcconfig(deck);
    std::vector<int> i1s = { 1, 10};
    std::vector<int> k2s = { 1, 4};
    std::vector<Opm::FaceDir::DirEnum> dirs = { Opm::FaceDir::XMinus , Opm::FaceDir::XPlus};
    for (const auto bc : bcconfig) {
        BOOST_CHECK(bc.i1 == (i1s[bc.index - 1 ]) - 1);
        BOOST_CHECK(bc.k2 == (k2s[bc.index - 1 ]) - 1);
        BOOST_CHECK(bc.dir == dirs[bc.index -1 ]);
    }

}
BOOST_AUTO_TEST_CASE(SimulationConfigBC) {
    BOOST_CHECK_THROW( BCConfig(createDeck(inputStr_bc)), OpmInputError );
}


BOOST_AUTO_TEST_CASE(SimulationConfigGetThresholdPressureTableTest) {
    auto deck = createDeck(inputStr);
    TableManager tm(deck);
    EclipseGrid eg(10, 3, 4);
    FieldPropsManager fp(deck, Phases{true, true, true}, eg, tm);
    BOOST_CHECK_NO_THROW( SimulationConfig(false, deck, fp) );
}


BOOST_AUTO_TEST_CASE(SimulationConfigNOTHPRES) {
    auto deck = createDeck(inputStr_noTHPRES);
    TableManager tm(deck);
    EclipseGrid eg(10, 3, 4);
    FieldPropsManager fp(deck, Phases{true, true, true}, eg, tm);
    SimulationConfig simulationConfig(false, deck, fp);
    BOOST_CHECK( !simulationConfig.useThresholdPressure() );
}

BOOST_AUTO_TEST_CASE(SimulationConfigCPRNotUsed) {
    auto deck = createDeck(inputStr_noTHPRES);
    TableManager tm(deck);
    EclipseGrid eg(10, 3, 4);
    FieldPropsManager fp(deck, Phases{true, true, true}, eg, tm);
    SimulationConfig simulationConfig(false, deck, fp);
    BOOST_CHECK( ! simulationConfig.useCPR());
}

BOOST_AUTO_TEST_CASE(SimulationConfigCPRUsed) {
    auto deck = createDeck(inputStr_cpr);
    TableManager tm(deck);
    EclipseGrid eg(10, 3, 4);
    SUMMARYSection summary(deck);
    FieldPropsManager fp(deck, Phases{true, true, true}, eg, tm);
    SimulationConfig simulationConfig(false, deck, fp);
    BOOST_CHECK(     simulationConfig.useCPR() );
    BOOST_CHECK(  !  summary.hasKeyword("CPR") );
}


BOOST_AUTO_TEST_CASE(SimulationConfigCPRInSUMMARYSection) {
    auto deck = createDeck(inputStr_cpr_in_SUMMARY);
    TableManager tm(deck);
    EclipseGrid eg(10, 3, 4);
    SUMMARYSection summary(deck);
    FieldPropsManager fp(deck, Phases{true, true, true}, eg, tm);
    SimulationConfig simulationConfig(false, deck, fp);
    BOOST_CHECK( ! simulationConfig.useCPR());
    BOOST_CHECK(   summary.hasKeyword("CPR"));
}


BOOST_AUTO_TEST_CASE(SimulationConfigCPRBoth) {
    auto deck = createDeck(inputStr_cpr_BOTH);
    TableManager tm(deck);
    EclipseGrid eg(10, 3, 4);
    SUMMARYSection summary(deck);
    FieldPropsManager fp(deck, Phases{true, true, true}, eg, tm);
    SimulationConfig simulationConfig(false, deck, fp);
    BOOST_CHECK(  simulationConfig.useCPR());
    BOOST_CHECK(  summary.hasKeyword("CPR"));

    const auto& cpr = summary.get<ParserKeywords::CPR>().back();
    const auto& record = cpr.getRecord(0);
    BOOST_CHECK_EQUAL( 1U , cpr.size());
    BOOST_CHECK_EQUAL( record.getItem<ParserKeywords::CPR::WELL>().get< std::string >(0) , "well1");
    BOOST_CHECK_EQUAL( record.getItem<ParserKeywords::CPR::I>().get< int >(0) , 10);
    BOOST_CHECK_EQUAL( record.getItem<ParserKeywords::CPR::J>().get< int >(0) , 20);
    BOOST_CHECK_EQUAL( record.getItem<ParserKeywords::CPR::K>().get< int >(0) , 30);
}


BOOST_AUTO_TEST_CASE(SimulationConfigCPRRUnspecWithData) {
    BOOST_CHECK_THROW( createDeck(inputStr_INVALID) , Opm::OpmInputError );
}

BOOST_AUTO_TEST_CASE(SimulationConfig_NONNC) {

    EclipseGrid eg(10, 3, 4);

    auto deck = createDeck(inputStr);
    TableManager tm(deck);
    FieldPropsManager fp(deck, Phases{true, true, true}, eg, tm);
    SimulationConfig simulationConfig(false, deck, fp);
    BOOST_CHECK_EQUAL( false , simulationConfig.useNONNC());

    auto deck_nonnc = createDeck(inputStr_nonnc);
    TableManager tm_nonnc(deck_nonnc);
    FieldPropsManager fp_nonnc(deck_nonnc, Phases{true, true, true}, eg, tm_nonnc);
    SimulationConfig simulationConfig_nonnc(false, deck_nonnc, fp_nonnc);
    BOOST_CHECK_EQUAL( true , simulationConfig_nonnc.useNONNC());
}

BOOST_AUTO_TEST_CASE(SimulationConfig_VAPOIL_DISGAS_VAPWAT) {

    EclipseGrid eg(10, 3, 4);

    auto deck = createDeck(inputStr);
    TableManager tm(deck);
    FieldPropsManager fp(deck, Phases{true, true, true}, eg, tm);
    SimulationConfig simulationConfig(false, deck, fp);
    BOOST_CHECK_EQUAL( false , simulationConfig.hasDISGAS());
    BOOST_CHECK_EQUAL( false , simulationConfig.hasVAPOIL());
    BOOST_CHECK_EQUAL( false , simulationConfig.hasVAPWAT());

    auto deck_vd = createDeck(inputStr_vap_dis);
    FieldPropsManager fp_vd(deck_vd, Phases{true, true, true}, eg, tm);
    SimulationConfig simulationConfig_vd(false, deck_vd, fp_vd);
    BOOST_CHECK_EQUAL( true , simulationConfig_vd.hasDISGAS());
    BOOST_CHECK_EQUAL( true , simulationConfig_vd.hasVAPOIL());
    BOOST_CHECK_EQUAL( true , simulationConfig_vd.hasVAPWAT());
    
}

BOOST_AUTO_TEST_CASE(SimulationConfig_TEMP_THERMAL)
{
    {
        const auto deck = createDeck(inputStr);
        const auto tm = TableManager(deck);
        auto eg = EclipseGrid(10, 3, 4);
        const auto fp = FieldPropsManager(deck, Phases{true, true, true}, eg, tm);
        const auto simulationConfig = Opm::SimulationConfig(false, deck, fp);

        BOOST_CHECK(! simulationConfig.isThermal());
    }

    {
        const auto deck = createDeck(simDeckStringTEMP());
        const auto tm = TableManager(deck);
        auto eg = EclipseGrid(10, 3, 4);
        const auto fp = FieldPropsManager(deck, Phases{true, true, true}, eg, tm);
        const auto simulationConfig = Opm::SimulationConfig(false, deck, fp);

        BOOST_CHECK(simulationConfig.isTemp());
    }

    {
        const auto deck = createDeck(simDeckStringTHERMAL());
        const auto tm = TableManager(deck);
        auto eg = EclipseGrid(10, 3, 4);
        const auto fp = FieldPropsManager(deck, Phases{true, true, true}, eg, tm);
        const auto simulationConfig = Opm::SimulationConfig(false, deck, fp);

        BOOST_CHECK(simulationConfig.isThermal());
    }
}


BOOST_AUTO_TEST_CASE(TESTRockConfig_Standard)
{
    const std::string deck_string = R"(
RUNSPEC

ROCKCOMP
/

TABDIMS
  * 3 /

PROPS

ROCK
   1  0.1 /
   2  0.2 /
   3  0.3 /

ROCKOPTS
   1* STORE SATNUM  /

)";
    Opm::Parser parser;
    const auto deck = parser.parseString(deck_string);
    EclipseGrid grid(10,10,10);
    FieldPropsManager fp(deck, Phases{true, true, true}, grid, TableManager());
    RockConfig rc(deck, fp);
    BOOST_CHECK_EQUAL(rc.rocknum_property(), "SATNUM");
    BOOST_CHECK_EQUAL(rc.store(), true);
    const auto& comp = rc.comp();
    BOOST_CHECK_EQUAL(comp.size(), 3U);
}

BOOST_AUTO_TEST_CASE(TESTRockConfig_Default)
{
    const auto deck = Parser{}.parseString(R"(
RUNSPEC

TABDIMS
  3  / -- NTSFUN = 3

PROPS

ROCKOPTS
  1* NOSTORE SATNUM /

ROCK
123.4 0.40E-05 /
/
271.8 1.61e-05 /

)");

    auto grid = EclipseGrid { 10, 10, 10 };
    const auto fp = FieldPropsManager {
        deck, Phases{true, true, true}, grid, TableManager()
    };

    const auto rc = RockConfig { deck, fp };

    BOOST_CHECK_EQUAL(rc.rocknum_property(), "SATNUM");
    BOOST_CHECK_EQUAL(rc.store(), false);

    const auto& comp = rc.comp();
    BOOST_REQUIRE_EQUAL(comp.size(), std::size_t{3});

    BOOST_CHECK_CLOSE(comp[0].pref, 123.4*1.0e5, 1.0e-8);
    BOOST_CHECK_CLOSE(comp[0].compressibility, 0.4e-5/1.0e5, 1.0e-8);

    BOOST_CHECK_CLOSE(comp[1].pref, 123.4*1.0e5, 1.0e-8);
    BOOST_CHECK_CLOSE(comp[1].compressibility, 0.4e-5/1.0e5, 1.0e-8);

    BOOST_CHECK_CLOSE(comp[2].pref, 271.8*1.0e5, 1.0e-8);
    BOOST_CHECK_CLOSE(comp[2].compressibility, 1.61e-5/1.0e5, 1.0e-8);
}

BOOST_AUTO_TEST_CASE(DatumDepth_Zero)
{
    const auto es = EclipseState { createDeck(R"(RUNSPEC
DIMENS
1 5 2 /
GRID
DXV
100.0 /
DYV
5*100.0 /
DZV
2*10.0 /
DEPTHZ
12*2000.0 /
END
)") };

    const auto& dd = es.getSimulationConfig().datumDepths();

    BOOST_CHECK_CLOSE(dd(0), 0.0, 1.0e-8);
}

BOOST_AUTO_TEST_CASE(DatumDepth_Equil)
{
    const auto es = EclipseState { createDeck(R"(RUNSPEC
DIMENS
1 5 2 /
EQLDIMS
/
GRID
DXV
100.0 /
DYV
5*100.0 /
DZV
2*10.0 /
DEPTHZ
12*2000.0 /
SOLUTION
EQUIL
  2005.0 123.4 2015.0 2.34 1995.0 0.0 /
END
)") };

    const auto& dd = es.getSimulationConfig().datumDepths();

    BOOST_CHECK_CLOSE(dd("FIPNUM", 0), 2005.0, 1.0e-8);
    BOOST_CHECK_CLOSE(dd("FIPABC", 42), 2005.0, 1.0e-8);
}
