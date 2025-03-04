/*
  Copyright 2014 Statoil ASA.

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
#include <opm/io/eclipse/PaddedOutputString.hpp>
#include <cstddef>
#include <cstdlib>
#include <optional>
#include <stdexcept>
#include <vector>
#include <string>
#define _USE_MATH_DEFINES

#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/EclipseState/Grid/LgrCollection.hpp>

#include <opm/input/eclipse/EclipseState/Grid/FieldProps.hpp>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/utility/numeric/calculateCellVol.hpp>
#include <opm/common/utility/numeric/VectorOps.hpp>
#include <opm/common/utility/String.hpp>

#include <opm/io/eclipse/EclFile.hpp>
#include <opm/io/eclipse/EclOutput.hpp>

#include <opm/input/eclipse/Units/Units.hpp>
#include <opm/input/eclipse/Deck/DeckSection.hpp>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Deck/DeckItem.hpp>
#include <opm/input/eclipse/Deck/DeckKeyword.hpp>
#include <opm/input/eclipse/Deck/DeckRecord.hpp>
#include <opm/input/eclipse/EclipseState/Grid/NNC.hpp>
#include <opm/common/utility/OpmInputError.hpp>
#include <opm/common/utility/numeric/VectorUtil.hpp>
#include <opm/common/utility/numeric/GeometryUtil.hpp>

#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <opm/input/eclipse/Parser/ParserKeywords/A.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/C.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/D.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/G.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/I.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/M.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/P.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/R.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/S.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/T.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/Z.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <cstdint>
#include <functional>
#include <initializer_list>
#include <numeric>
#include <tuple>
#include <stdexcept>
#include <vector>

#include <fmt/format.h>

namespace Opm {

namespace {


std::optional<UnitSystem> make_grid_units(const std::string& grid_unit) {
    if (grid_unit == "METRES")
        return UnitSystem(UnitSystem::UnitType::UNIT_TYPE_METRIC);

    if (grid_unit == "FEET")
        return UnitSystem(UnitSystem::UnitType::UNIT_TYPE_FIELD);

    if (grid_unit == "CM")
        return UnitSystem(UnitSystem::UnitType::UNIT_TYPE_LAB);

    return std::nullopt;
}

void apply_GRIDUNIT(const UnitSystem& deck_units, const UnitSystem& grid_units, std::vector<double>& data)
{
    double scale_factor = grid_units.getDimension(UnitSystem::measure::length).getSIScaling() / deck_units.getDimension(UnitSystem::measure::length).getSIScaling();
    std::transform(data.begin(), data.end(), data.begin(),
                   [scale_factor](const auto& v) { return v * scale_factor; });
}

}
EclipseGrid::EclipseGrid()
    : GridDims(),
      m_minpvMode(MinpvMode::Inactive),
      m_pinchoutMode(PinchMode::TOPBOT),
      m_multzMode(PinchMode::TOP),
      m_pinchGapMode(PinchMode::GAP),
      m_pinchMaxEmptyGap(ParserKeywords::PINCH::MAX_EMPTY_GAP::defaultValue)
{}

EclipseGrid::EclipseGrid(const std::array<int, 3>& dims ,
                         const std::vector<double>& coord ,
                         const std::vector<double>& zcorn ,
                         const int * actnum)
    : GridDims(dims),
      m_minpvMode(MinpvMode::Inactive),
      m_pinchoutMode(PinchMode::TOPBOT),
      m_multzMode(PinchMode::TOP),
      m_pinchGapMode(PinchMode::GAP),
      m_pinchMaxEmptyGap(ParserKeywords::PINCH::MAX_EMPTY_GAP::defaultValue)
{
    initCornerPointGrid( coord , zcorn , actnum );
}

/**
   Will create an EclipseGrid instance based on an existing
   GRID/EGRID file.
*/


EclipseGrid::EclipseGrid(const std::string& fileName )
    : GridDims(),
      m_minpvMode(MinpvMode::Inactive),
      m_pinchoutMode(PinchMode::TOPBOT),
      m_multzMode(PinchMode::TOP),
      m_pinchGapMode(PinchMode::GAP),
      m_pinchMaxEmptyGap(ParserKeywords::PINCH::MAX_EMPTY_GAP::defaultValue)
{

    Opm::EclIO::EclFile egridfile(fileName);
    m_useActnumFromGdfile=true;
    initGridFromEGridFile(egridfile, fileName);
}


EclipseGrid::EclipseGrid(const GridDims& gd)
    : GridDims(gd),
      m_minpvMode(MinpvMode::Inactive),
      m_pinchoutMode(PinchMode::TOPBOT),
      m_multzMode(PinchMode::TOP),
      m_pinchGapMode(PinchMode::GAP),
      m_pinchMaxEmptyGap(ParserKeywords::PINCH::MAX_EMPTY_GAP::defaultValue)
{
    this->m_nactive = this->getCartesianSize();
    this->active_volume = std::nullopt;
    // Nothing else initialized. Leaving in particular as empty:
    // m_actnum,
    // m_global_to_active,
    // m_active_to_global.
    //
    // EclipseGrid will only be usable for constructing Box objects with
    // all cells active.
}


EclipseGrid::EclipseGrid(std::size_t nx, std::size_t ny , std::size_t nz,
                         double dx, double dy, double dz, double top)
    : GridDims(nx, ny, nz),
      m_minpvMode(MinpvMode::Inactive),
      m_pinchoutMode(PinchMode::TOPBOT),
      m_multzMode(PinchMode::TOP),
      m_pinchGapMode(PinchMode::GAP),
      m_pinchMaxEmptyGap(ParserKeywords::PINCH::MAX_EMPTY_GAP::defaultValue)
{

    m_coord.reserve((nx+1)*(ny+1)*6);

    for (std::size_t j = 0; j < (ny+1); j++) {
        for (std::size_t i = 0; i < (nx+1); i++) {
            m_coord.push_back(i*dx);
            m_coord.push_back(j*dy);
            m_coord.push_back(0.0);
            m_coord.push_back(i*dx);
            m_coord.push_back(j*dy);
            m_coord.push_back(nz*dz);
        }
    }

    m_zcorn.assign(nx*ny*nz*8, 0);

    for (std::size_t k = 0; k < nz ; k++) {
        for (std::size_t j = 0; j < ny ; j++) {
            for (std::size_t i = 0; i < nx ; i++) {

                // top face of cell
                int zind = i*2 + j*nx*4 + k*nx*ny*8;

                double zt = top + k*dz;
                double zb = top + (k+1)*dz;

                m_zcorn[zind] = zt;
                m_zcorn[zind + 1] = zt;

                zind = zind + nx*2;

                m_zcorn[zind] = zt;
                m_zcorn[zind + 1] = zt;

                // bottom face of cell
                zind = i*2 + j*nx*4 + k*nx*ny*8 + nx*ny*4;

                m_zcorn[zind] = zb;
                m_zcorn[zind + 1] = zb;

                zind = zind + nx*2;

                m_zcorn[zind] = zb;
                m_zcorn[zind + 1] = zb;
            }
        }
    }

    resetACTNUM();
}

EclipseGrid::EclipseGrid(const EclipseGrid& src, const double* zcorn, const std::vector<int>& actnum)
    : EclipseGrid(src)
{

    if (zcorn != nullptr) {
        std::size_t sizeZcorn = this->getCartesianSize()*8;

        for (std::size_t n=0; n < sizeZcorn; n++) {
            m_zcorn[n] = zcorn[n];
        }

        ZcornMapper mapper( getNX(), getNY(), getNZ());
        zcorn_fixed = mapper.fixupZCORN( m_zcorn );
    }

    resetACTNUM(actnum);
}

EclipseGrid::EclipseGrid(const EclipseGrid& src, const std::vector<int>& actnum)
        : EclipseGrid( src , nullptr , actnum )
{ }

/*
  This is the main EclipseGrid constructor, it will inspect the
  input Deck for grid keywords, either the corner point keywords
  COORD and ZCORN, or the various rectangular keywords like DX,DY
  and DZ.

  Actnum is treated specially:

    1. If an actnum pointer is passed in that should be a pointer
       to 0 and 1 values which will be used as actnum mask.

    2. If the actnum pointer is not present the constructor will
       look in the deck for an actnum keyword, and use that if it
       is found. This is a best effort which will work in many
       cases, but if the ACTNUM keyword is manipulated in the deck
       those manipulations will be silently lost; if the ACTNUM
       keyword has size different from nx*ny*nz it will also be
       silently ignored.

  With a mutable EclipseGrid instance you can later call the
  EclipseGrid::resetACTNUM() method when you have complete actnum
  information. The EclipseState based construction of EclipseGrid
  is a two-pass operation, which guarantees that actnum is handled
  correctly.
*/



EclipseGrid::EclipseGrid(const Deck& deck, const int * actnum)
    : GridDims(deck),
      m_minpvMode(MinpvMode::Inactive),
      m_pinchoutMode(PinchMode::TOPBOT),
      m_multzMode(PinchMode::TOP),
      m_pinchGapMode(PinchMode::GAP),
      m_pinchMaxEmptyGap(ParserKeywords::PINCH::MAX_EMPTY_GAP::defaultValue)
{
    if (deck.hasKeyword("GDFILE")){

        if (deck.hasKeyword("ACTNUM")){
            if (keywInputBeforeGdfile(deck, "ACTNUM"))  {
                m_useActnumFromGdfile = true;
            }
        } else {
            m_useActnumFromGdfile = true;
        }

    }

    updateNumericalAquiferCells(deck);

    initGrid(deck, actnum);

    if (deck.hasKeyword<ParserKeywords::MAPAXES>())
        this->m_mapaxes = std::make_optional<MapAxes>( deck );

    /*
      The GRIDUNIT handling is simplified compared to the full specification:

        1. The optional second item 'MAP' is ignored.

        2. In Eclipse the action of the GRIDUNIT keyword only applies to
           keywords in the same file as the GRIDUNIT keyword itself is in; in
           our implementation we apply the GRIDUNIT transformation if the deck
           contains the GRIDUNUT keyword - irrespective of exactly where it is
           located in the deck.
    */

    if (deck.hasKeyword<ParserKeywords::GRIDUNIT>()) {
        const auto& kw = deck.get<ParserKeywords::GRIDUNIT>().front();
        const auto& length_unit = trim_copy(kw.getRecord(0).getItem(0).get<std::string>(0));
        auto grid_units = make_grid_units(length_unit);
        if (!grid_units.has_value())
            throw OpmInputError(fmt::format("Invalid length specifier: [{}]", length_unit), kw.location());

        if (grid_units.value() != deck.getActiveUnitSystem()) {
            apply_GRIDUNIT(deck.getActiveUnitSystem(), grid_units.value(), this->m_zcorn);
            apply_GRIDUNIT(deck.getActiveUnitSystem(), grid_units.value(), this->m_coord);
            if (this->m_rv.has_value())
                apply_GRIDUNIT(deck.getActiveUnitSystem(), grid_units.value(), this->m_rv.value());

            if (this->m_input_coord.has_value()) {
                apply_GRIDUNIT(deck.getActiveUnitSystem(), grid_units.value(), this->m_input_coord.value());
            }

            if (this->m_input_zcorn.has_value()) {
                apply_GRIDUNIT(deck.getActiveUnitSystem(), grid_units.value(), this->m_input_zcorn.value());
            }
        }
    }
}


    bool EclipseGrid::circle( ) const{
        return this->m_circle;
    }


    void EclipseGrid::initGrid(const Deck& deck, const int* actnum)
    {
        enum GridType { COORD, DEPTHZ, TOPS, RADIAL, SPIDER, GDFILE, GT_SIZE };

        std::vector<int> found;
        if (hasCornerPointKeywords(deck)) {
            found.push_back(GridType::COORD);
        }
        if (hasDVDEPTHZKeywords(deck)) {
            found.push_back(GridType::DEPTHZ);
        }
        if (hasDTOPSKeywords(deck)) {
            found.push_back(GridType::TOPS);
        }
        if (hasRadialKeywords(deck)) {
            found.push_back(GridType::RADIAL);
        }
        if (hasSpiderKeywords(deck)) {
            found.push_back(GridType::SPIDER);
        }
        if (hasGDFILE(deck)) {
            found.push_back(GridType::GDFILE);
        }

        const std::string messages[GT_SIZE] {
            "COORD with ZCORN creates a corner-point grid",
            "DEPTHZ with DXV, DYV, DZV creates a cartesian grid",
            "TOPS with DX/DXV, DY/DYV, DZ/DZV creates a cartesian grid",
            "RADIAL with DR/DRV, DTHETA/DTHETAV, DZ/DZV and TOPS creates a cylindrical grid",
            "SPIDER with DR/DRV, DTHETA/DTHETAV, DZ/DZV and TOPS creates a spider grid",
            "GDFILE reads a grid from file"};

        if (found.empty()) {
            std::string message = "The grid must be specified using one of these options:";
            for (int grid_type = 0; grid_type < GridType::GT_SIZE; grid_type++) {
                message += "\n    " + messages[grid_type];
            }
            throw std::invalid_argument(message);
        }

        if (found.size() > 1) {
            std::string message = "The specification of the grid is ambiguous:";
            for (const int& grid_type : found) {
                message += "\n    " + messages[grid_type];
            }
            throw std::invalid_argument(message);
        }

        switch (found.front()) {
        case GridType::COORD:
            this->initCornerPointGrid(deck);
            break;
        case GridType::DEPTHZ:
        case GridType::TOPS:
            this->initCartesianGrid(deck);
            break;
        case GridType::RADIAL:
            this->initCylindricalGrid(deck);
            break;
        case GridType::SPIDER:
            this->initSpiderwebGrid(deck);
            break;
        case GridType::GDFILE:
            this->initBinaryGrid(deck);
            break;
        }

        if (deck.hasKeyword<ParserKeywords::PINCH>()) {
            const auto& record = deck.get<ParserKeywords::PINCH>( ).back().getRecord(0);
            const auto& item = record.getItem<ParserKeywords::PINCH::THRESHOLD_THICKNESS>( );
            m_pinch = item.getSIDouble(0);

            auto pinchoutString = record.getItem<ParserKeywords::PINCH::PINCHOUT_OPTION>().get< std::string >(0);
            m_pinchoutMode = PinchModeFromString(pinchoutString);

            auto multzString = record.getItem<ParserKeywords::PINCH::MULTZ_OPTION>().get< std::string >(0);
            m_multzMode = PinchModeFromString(multzString);
            auto pinchGapString = record.getItem<ParserKeywords::PINCH::CONTROL_OPTION>().get< std::string >(0);
            m_pinchGapMode = PinchModeFromString(pinchGapString);
            m_pinchMaxEmptyGap = record.getItem<ParserKeywords::PINCH::MAX_EMPTY_GAP>().getSIDouble(0);
        }

        m_minpvVector.resize(getCartesianSize(), 0.0);
        if (deck.hasKeyword<ParserKeywords::MINPV>()) {
            const auto& record = deck.get<ParserKeywords::MINPV>( ).back().getRecord(0);
            const auto& item = record.getItem<ParserKeywords::MINPV::VALUE>( );
            std::fill(m_minpvVector.begin(), m_minpvVector.end(), item.getSIDouble(0));
            m_minpvMode = MinpvMode::EclSTD;
        } else if (deck.hasKeyword<ParserKeywords::MINPORV>()) {
            const auto& record = deck.get<ParserKeywords::MINPORV>( ).back().getRecord(0);
            const auto& item = record.getItem<ParserKeywords::MINPORV::VALUE>( );
            std::fill(m_minpvVector.begin(), m_minpvVector.end(), item.getSIDouble(0));
            m_minpvMode = MinpvMode::EclSTD;
        }
        // Note that MINPVV is not handled here but in a second stage in EclipseState where we
        // call setMINPVV to also support BOX via grid properties

        if (actnum != nullptr) {
            this->resetACTNUM(actnum);
        }
        else if (! this->m_useActnumFromGdfile) {
            const auto fp = FieldProps {
                deck, EclipseGrid { static_cast<GridDims&>(*this) }
            };

            this->resetACTNUM(fp.actnumRaw());
        }
    }

    void EclipseGrid::initGridFromEGridFile(Opm::EclIO::EclFile& egridfile,
                                            const std::string&   fileName)
    {
        OpmLog::info(fmt::format("\nCreating grid from: {} ", fileName));

        if (!egridfile.hasKey("GRIDHEAD")) {
            throw std::invalid_argument("file: " + fileName + " is not a valid egrid file, GRIDHEAD not found");
        }

        if (!egridfile.hasKey("COORD")) {
            throw std::invalid_argument("file: " + fileName + " is not a valid egrid file, COORD not found");
        }

        if (!egridfile.hasKey("ZCORN")) {
            throw std::invalid_argument("file: " + fileName + " is not a valid egrid file, ZCORN not found");
        }

        if (!egridfile.hasKey("GRIDUNIT")) {
            throw std::invalid_argument("file: " + fileName + " is not a valid egrid file, GRIDUNIT not found");
        }

        {
            const std::vector<int>& gridhead = egridfile.get<int>("GRIDHEAD");

            this->m_nx = gridhead[1];
            this->m_ny = gridhead[2];
            this->m_nz = gridhead[3];
        }

        {
            const std::vector<float>& coord_f = egridfile.get<float>("COORD");
            const std::vector<float>& zcorn_f = egridfile.get<float>("ZCORN");

            m_coord.assign(coord_f.begin(), coord_f.end());
            m_zcorn.assign(zcorn_f.begin(), zcorn_f.end());
        }

        if (const auto& gridunit = egridfile.get<std::string>("GRIDUNIT");
            gridunit[0] != "METRES")
        {
            constexpr auto length = ::Opm::UnitSystem::measure::length;

            if (gridunit[0] == "FEET") {
                Opm::UnitSystem units(Opm::UnitSystem::UnitType::UNIT_TYPE_FIELD);
                units.to_si(length, m_coord);
                units.to_si(length, m_zcorn);
            }
            else if (gridunit[0] == "CM") {
                Opm::UnitSystem units(Opm::UnitSystem::UnitType::UNIT_TYPE_LAB);
                units.to_si(length, m_coord);
                units.to_si(length, m_zcorn);
            }
            else {
                std::string message = "gridunit '" + gridunit[0] + "' doesn't correspong to a valid unit system";
                throw std::invalid_argument(message);
            }
        }

        if (egridfile.hasKey("ACTNUM") && m_useActnumFromGdfile) {
            resetACTNUM(egridfile.get<int>("ACTNUM"));
        }
        else {
            this->resetACTNUM();
        }

        if (egridfile.hasKey("MAPAXES")) {
            this->m_mapaxes = std::make_optional<MapAxes>(egridfile);
        }

        ZcornMapper mapper(getNX(), getNY(), getNZ());
        zcorn_fixed = mapper.fixupZCORN(m_zcorn);
    }

    bool EclipseGrid::keywInputBeforeGdfile(const Deck& deck, const std::string& keyword) const {

        std::vector<std::string> keywordList;
        keywordList.reserve(deck.size());

        for (std::size_t n=0;n<deck.size();n++){
            DeckKeyword kw = deck[n];
            keywordList.push_back(kw.name());
        }

        int indKeyw = -1;
        int indGdfile = -1;

        for (std::size_t i = 0; i < keywordList.size(); i++){
            if (keywordList[i]=="GDFILE"){
                indGdfile=i;
            }

            if (keywordList[i]==keyword){
                indKeyw=i;
            }
        }

        if (indGdfile == -1) {
            throw std::runtime_error("keyword GDFILE not found in deck ");
        }

        if (indKeyw == -1) {
            throw std::runtime_error("keyword " + keyword + " not found in deck ");
        }

        return indKeyw < indGdfile;
    }


    std::size_t EclipseGrid::activeIndex(std::size_t i, std::size_t j, std::size_t k) const {

        return activeIndex( getGlobalIndex( i,j,k ));
    }

    std::size_t EclipseGrid::activeIndex(std::size_t globalIndex) const {
        if (m_global_to_active.empty()) {
            return globalIndex;
        }

        if (m_global_to_active[ globalIndex] == -1) {
            throw std::invalid_argument("Input argument does not correspond to an active cell");
        }

        return m_global_to_active[ globalIndex];
    }

    /**
       Observe: the input argument is assumed to be in the space
       [0,num_active).
    */

    std::size_t EclipseGrid::getGlobalIndex(std::size_t active_index) const {
        return m_active_to_global.at(active_index);
    }

    bool EclipseGrid::isPinchActive( ) const {
        return m_pinch.has_value();
    }

    double EclipseGrid::getPinchThresholdThickness( ) const {
        return m_pinch.value();
    }

    PinchMode EclipseGrid::getPinchOption( ) const {
        return m_pinchoutMode;
    }

    PinchMode EclipseGrid::getMultzOption( ) const {
        return m_multzMode;
    }

    MinpvMode EclipseGrid::getMinpvMode() const {
        return m_minpvMode;
    }

    PinchMode EclipseGrid::getPinchGapMode() const {
        return m_pinchGapMode;
    }

    double EclipseGrid::getPinchMaxEmptyGap() const {
        return m_pinchMaxEmptyGap;
    }


    const std::vector<double>& EclipseGrid::getMinpvVector( ) const {
        return m_minpvVector;
    }

    void EclipseGrid::initBinaryGrid(const Deck& deck)
    {
        const auto& gdfile = deck["GDFILE"].back().getRecord(0);

        const auto formatted = EclIO::EclFile::Formatted {
            gdfile.getItem("formatted").getTrimmedString(0).front() == 'F'
        };

        auto filename = deck
            .makeDeckPath(gdfile.getItem("filename").getTrimmedString(0));

        std::unique_ptr<EclIO::EclFile> egridfile{};

        // Some windows applications export .DATA files with relative GDFILE
        // keywords with a windows formatted path.  If open fails we give it
        // one more try with the replacement '\'' -> '/'.
        try {
            egridfile = std::make_unique<Opm::EclIO::EclFile>(filename, formatted);
        }
        catch (const std::runtime_error&) {
            std::replace(filename.begin(), filename.end(), '\\', '/');
            egridfile = std::make_unique<Opm::EclIO::EclFile>(filename, formatted);
        }

        this->initGridFromEGridFile(*egridfile, filename);
    }

    void EclipseGrid::initCartesianGrid(const Deck& deck) {

        if (hasDVDEPTHZKeywords( deck )) {
            initDVDEPTHZGrid(deck );
        } else if (hasDTOPSKeywords(deck)) {
            initDTOPSGrid( deck );
        } else {
            throw std::invalid_argument("Tried to initialize cartesian grid without all required keywords");
        }
    }

    void EclipseGrid::initDVDEPTHZGrid(const Deck& deck) {
      OpmLog::info(fmt::format("\nCreating grid from keywords DXV, DYV, DZV and DEPTHZ"));
        const std::vector<double>& DXV = deck.get<ParserKeywords::DXV>().back().getSIDoubleData();
        const std::vector<double>& DYV = deck.get<ParserKeywords::DYV>().back().getSIDoubleData();
        const std::vector<double>& DZV = deck.get<ParserKeywords::DZV>().back().getSIDoubleData();
        const std::vector<double>& DEPTHZ = deck.get<ParserKeywords::DEPTHZ>().back().getSIDoubleData();
        auto nx = this->getNX();
        auto ny = this->getNY();
        auto nz = this->getNZ();

        assertVectorSize( DEPTHZ , static_cast<size_t>( (nx + 1)*(ny +1 )) , "DEPTHZ");
        assertVectorSize( DXV    , static_cast<size_t>( nx ) , "DXV");
        assertVectorSize( DYV    , static_cast<size_t>( ny ) , "DYV");
        assertVectorSize( DZV    , static_cast<size_t>( nz ) , "DZV");

        m_coord = makeCoordDxvDyvDzvDepthz(DXV, DYV, DZV, DEPTHZ);
        m_zcorn = makeZcornDzvDepthz(DZV, DEPTHZ);

        ZcornMapper mapper( getNX(), getNY(), getNZ());
        zcorn_fixed = mapper.fixupZCORN( m_zcorn );
    }

    void EclipseGrid::initDTOPSGrid(const Deck& deck) {
      OpmLog::info(fmt::format("\nCreating grid from keywords DX, DY, DZ and TOPS"));

        std::vector<double> DX = EclipseGrid::createDVector(  this->getNXYZ(), 0 , "DX" , "DXV" , deck);
        std::vector<double> DY = EclipseGrid::createDVector(  this->getNXYZ(), 1 , "DY" , "DYV" , deck);
        std::vector<double> DZ = EclipseGrid::createDVector(  this->getNXYZ(), 2 , "DZ" , "DZV" , deck);
        std::vector<double> TOPS = EclipseGrid::createTOPSVector( this->getNXYZ(), DZ , deck );

        m_coord = makeCoordDxDyDzTops(DX, DY, DZ, TOPS);
        m_zcorn = makeZcornDzTops(DZ, TOPS);

        ZcornMapper mapper( getNX(), getNY(), getNZ());
        zcorn_fixed = mapper.fixupZCORN( m_zcorn );
    }


    void EclipseGrid::getCellCorners(const std::array<int, 3>& ijk, const std::array<int, 3>& dims,
                                     std::array<double,8>& X,
                                     std::array<double,8>& Y,
                                     std::array<double,8>& Z) const
    {

        std::array<int, 8> zind;
        std::array<int, 4> pind;

        // calculate indices for grid pillars in COORD arrray
        const std::size_t p_offset = ijk[1]*(dims[0]+1)*6 + ijk[0]*6;

        pind[0] = p_offset;
        pind[1] = p_offset + 6;
        pind[2] = p_offset + (dims[0]+1)*6;
        pind[3] = pind[2] + 6;

        // get depths from zcorn array in ZCORN array
        const std::size_t z_offset = ijk[2]*dims[0]*dims[1]*8 + ijk[1]*dims[0]*4 + ijk[0]*2;

        zind[0] = z_offset;
        zind[1] = z_offset + 1;
        zind[2] = z_offset + dims[0]*2;
        zind[3] = zind[2] + 1;

        for (int n = 0; n < 4; n++)
            zind[n+4] = zind[n] + dims[0]*dims[1]*4;


        for (int n = 0; n< 8; n++)
           Z[n] = m_zcorn[zind[n]];


        for (int  n=0; n<4; n++) {
            const double xt = m_coord[pind[n]];
            const double yt = m_coord[pind[n] + 1];
            const double zt = m_coord[pind[n] + 2];

            const double zb = m_coord[pind[n]+5];

            if (zt == zb) {
                X[n] = xt;
                X[n + 4] = xt;

                Y[n] = yt;
                Y[n + 4] = yt;
            } else {
                const double xb = m_coord[pind[n] + 3];
                const double yb = m_coord[pind[n] + 4];
                X[n] = xt + (xb-xt) / (zt-zb) * (zt - Z[n]);
                X[n+4] = xt + (xb-xt) / (zt-zb) * (zt-Z[n+4]);

                Y[n] = yt+(yb-yt)/(zt-zb)*(zt-Z[n]);
                Y[n+4] = yt+(yb-yt)/(zt-zb)*(zt-Z[n+4]);
            }
        }
    }


    void EclipseGrid::getCellCorners(const std::size_t globalIndex,
                                     std::array<double,8>& X,
                                     std::array<double,8>& Y,
                                     std::array<double,8>& Z) const
    {
        this->assertGlobalIndex(globalIndex);
        auto ijk = this->getIJK(globalIndex);
        this->getCellCorners(ijk, this->getNXYZ(), X, Y, Z);
    }


    std::vector<double> EclipseGrid::makeCoordDxvDyvDzvDepthz(const std::vector<double>& dxv, const std::vector<double>& dyv, const std::vector<double>& dzv, const std::vector<double>& depthz) const {
        auto nx = this->getNX();
        auto ny = this->getNY();
        auto nz = this->getNZ();

        std::vector<double> coord;
        coord.reserve((nx+1)*(ny+1)*6);

        std::vector<double> x(nx + 1, 0.0);
        std::partial_sum(dxv.begin(), dxv.end(), x.begin() + 1);

        std::vector<double> y(ny + 1, 0.0);
        std::partial_sum(dyv.begin(), dyv.end(), y.begin() + 1);

        std::vector<double> z(nz + 1, 0.0);
        std::partial_sum(dzv.begin(), dzv.end(), z.begin() + 1);


        for (std::size_t j = 0; j < ny + 1; j++) {
            for (std::size_t i = 0; i<nx + 1; i++) {

                const double x0 = x[i];
                const double y0 = y[j];

                std::size_t ind = i+j*(nx+1);

                double zb = depthz[ind] + z[nz];

                coord.push_back(x0);
                coord.push_back(y0);
                coord.push_back(depthz[ind]);
                coord.push_back(x0);
                coord.push_back(y0);
                coord.push_back(zb);
            }
        }

        return coord;
    }

    std::vector<double> EclipseGrid::makeZcornDzvDepthz(const std::vector<double>& dzv, const std::vector<double>& depthz) const {
        auto nx = this->getNX();
        auto ny = this->getNY();
        auto nz = this->getNZ();

        std::vector<double> zcorn;
        std::size_t sizeZcorn = nx*ny*nz*8;

        zcorn.assign (sizeZcorn, 0.0);

        std::vector<double> z(nz + 1, 0.0);
        std::partial_sum(dzv.begin(), dzv.end(), z.begin() + 1);


        for (std::size_t k = 0; k < nz; k++) {
            for (std::size_t j = 0; j < ny; j++) {
                for (std::size_t i = 0; i < nx; i++) {

                    const double z0 = z[k];

                    // top face of cell
                    auto zind = i*2 + j*nx*4 + k*nx*ny*8;

                    zcorn[zind] = depthz[i+j*(nx+1)] + z0;
                    zcorn[zind + 1] = depthz[i+j*(nx+1) +1 ] + z0;

                    zind = zind + nx*2;

                    zcorn[zind] = depthz[i+(j+1)*(nx+1)] + z0;
                    zcorn[zind + 1] = depthz[i+(j+1)*(nx+1) + 1] + z0;

                    // bottom face of cell
                    zind = i*2 + j*nx*4 + k*nx*ny*8 + nx*ny*4;

                    zcorn[zind] = depthz[i+j*(nx+1)] + z0 + dzv[k];
                    zcorn[zind + 1] = depthz[i+j*(nx+1) +1 ] + z0 + dzv[k];

                    zind = zind + nx*2;

                    zcorn[zind] = depthz[i+(j+1)*(nx+1)] + z0 + dzv[k];
                    zcorn[zind + 1] = depthz[i+(j+1)*(nx+1) + 1] + z0 + dzv[k];
                }
            }
        }

        return zcorn;
    }

    namespace
    {
        std::vector<double> makeSumIdirAtK(const int nx, const int ny, const int k, const std::vector<double>& dx)
        {
            std::vector<double> s(nx * ny, 0.0);
            for (int j = 0; j < ny; ++j) {
                double sum = 0.0;
                for (int i = 0; i < nx; ++i) {
                    sum += dx[i + j*nx + k*nx*ny];
                    s[i + j*nx] = sum;
                }
            }
            return s;
        }

        std::vector<double> makeSumJdirAtK(const int nx, const int ny, const int k, const std::vector<double>& dy)
        {
            std::vector<double> s(nx * ny, 0.0);
            for (int i = 0; i < nx; ++i) {
                double sum = 0.0;
                for (int j = 0; j < ny; ++j) {
                    sum += dy[i + j*nx + k*nx*ny];
                    s[i + j*nx] = sum;
                }
            }
            return s;
        }

        std::vector<double> makeSumKdir(const int nx, const int ny, const int nz, const std::vector<double>& dz)
        {
            std::vector<double> s(nx * ny, 0.0);
            for (int i = 0; i < nx; ++i) {
                for (int j = 0; j < ny; ++j) {
                    double sum = 0.0;
                    for (int k = 0; k < nz; ++k) {
                        sum += dz[i + j*nx + k*nx*ny];
                    }
                    s[i + j*nx] = sum;
                }
            }
            return s;
        }

    } // anonymous namespace

    std::vector<double> EclipseGrid::makeCoordDxDyDzTops(const std::vector<double>& dx, const std::vector<double>& dy, const std::vector<double>& dz, const std::vector<double>& tops) const {
        auto nx = this->getNX();
        auto ny = this->getNY();
        auto nz = this->getNZ();

        std::vector<double> coord;
        coord.reserve((nx+1)*(ny+1)*6);

        std::vector<double> sum_idir_top = makeSumIdirAtK(nx, ny, 0, dx);
        std::vector<double> sum_idir_bot = makeSumIdirAtK(nx, ny, nz - 1, dx);
        std::vector<double> sum_jdir_top = makeSumJdirAtK(nx, ny, 0, dy);
        std::vector<double> sum_jdir_bot = makeSumJdirAtK(nx, ny, nz - 1, dy);
        std::vector<double> sum_kdir = makeSumKdir(nx, ny, nz, dz);

        for (std::size_t j = 0; j < ny; j++) {

            double zt = tops[0];
            double zb = zt + sum_kdir[0 + 0*nx];

            if (j == 0) {
                double x0 = 0.0;
                double y0 = 0;

                coord.push_back(x0);
                coord.push_back(y0);
                coord.push_back(zt);
                coord.push_back(x0);
                coord.push_back(y0);
                coord.push_back(zb);

                for (std::size_t i = 0; i < nx; i++) {

                    std::size_t ind = i+j*nx+1;

                    if (i == (nx-1)) {
                        ind=ind-1;
                    }

                    zt = tops[ind];
                    zb = zt + sum_kdir[i + j*nx];

                    double xt = x0 + dx[i + j*nx];
                    double xb = sum_idir_bot[i + j*nx];

                    coord.push_back(xt);
                    coord.push_back(y0);
                    coord.push_back(zt);
                    coord.push_back(xb);
                    coord.push_back(y0);
                    coord.push_back(zb);

                    x0=xt;
                }
            }

            std::size_t ind = (j+1)*nx;

            if (j == (ny-1) ) {
                ind = j*nx;
            }

            double x0 = 0.0;

            double yt = sum_jdir_top[0 + j*nx];
            double yb = sum_jdir_bot[0 + j*nx];

            zt = tops[ind];
            zb = zt + sum_kdir[0 + j*nx];

            coord.push_back(x0);
            coord.push_back(yt);
            coord.push_back(zt);
            coord.push_back(x0);
            coord.push_back(yb);
            coord.push_back(zb);

            for (std::size_t i=0; i < nx; i++) {

                ind = i+(j+1)*nx+1;

                if (j == (ny-1) ) {
                    ind = i+j*nx+1;
                }

                if (i == (nx - 1) ) {
                    ind=ind-1;
                }

                zt = tops[ind];
                zb = zt + sum_kdir[i + j*nx];

                double xt=-999;
                double xb;

                if (j == (ny-1) ) {
                    xt = sum_idir_top[i + j*nx];
                    xb = sum_idir_bot[i + j*nx];
                } else {
                    xt = sum_idir_top[i + (j+1)*nx];
                    xb = sum_idir_bot[i + (j+1)*nx];
                }

                if (i == (nx - 1) ) {
                    yt = sum_jdir_top[i + j*nx];
                    yb = sum_jdir_bot[i + j*nx];
                } else {
                    yt = sum_jdir_top[(i + 1) + j*nx];
                    yb = sum_jdir_bot[(i + 1) + j*nx];
                }

                coord.push_back(xt);
                coord.push_back(yt);
                coord.push_back(zt);
                coord.push_back(xb);
                coord.push_back(yb);
                coord.push_back(zb);
            }
        }

        return coord;
    }

    std::vector<double> EclipseGrid::makeZcornDzTops(const std::vector<double>& dz, const std::vector<double>& tops) const {

        std::vector<double> zcorn;
        std::size_t sizeZcorn = this->getCartesianSize() * 8;
        auto nx = this->getNX();
        auto ny = this->getNY();
        auto nz = this->getNZ();
        zcorn.assign (sizeZcorn, 0.0);

        for (std::size_t j = 0; j < ny; j++) {
            for (std::size_t i = 0; i < nx; i++) {
                std::size_t ind = i + j*nx;
                double z = tops[ind];

                for (std::size_t k = 0; k < nz; k++) {

                    // top face of cell
                    std::size_t zind = i*2 + j*nx*4 + k*nx*ny*8;

                    zcorn[zind] = z;
                    zcorn[zind + 1] = z;

                    zind = zind + nx*2;

                    zcorn[zind] = z;
                    zcorn[zind + 1] = z;

                    z = z + dz[i + j*nx + k*nx*ny];

                    // bottom face of cell
                    zind = i*2 + j*nx*4 + k*nx*ny*8 + nx*ny*4;

                    zcorn[zind] = z;
                    zcorn[zind + 1] = z;

                    zind = zind + nx*2;

                    zcorn[zind] = z;
                    zcorn[zind + 1] = z;
                }
            }
        }

        return zcorn;
    }

    void EclipseGrid::initCylindricalGrid(const Deck& deck)
    {
        initSpiderwebOrCylindricalGrid(deck, true);
    }

    void EclipseGrid::initSpiderwebGrid(const Deck& deck)
    {
        initSpiderwebOrCylindricalGrid(deck, false);
    }

    /*
      Limited implementaton - requires keywords: DRV, DTHETAV, DZV or DZ, and TOPS.
    */

    void EclipseGrid::initSpiderwebOrCylindricalGrid(const Deck& deck, const bool is_cylindrical)
    {
        const std::string kind = is_cylindrical ? "cylindrical" : "spiderweb";

        // The hasCylindricalKeywords( ) checks according to the
        // eclipse specification for RADIAL grid. We currently do not support all
        // aspects of cylindrical grids, we therefor have an
        // additional test here, which checks if we have the keywords
        // required by the current implementation.
        if (!hasCylindricalKeywords(deck))
            throw std::invalid_argument("Not all keywords required for " + kind + " grids present");

        if (!deck.hasKeyword<ParserKeywords::DTHETAV>())
            throw std::logic_error("The current implementation *must* have theta values specified using the DTHETAV keyword");

        if (!deck.hasKeyword<ParserKeywords::DRV>())
            throw std::logic_error("The current implementation *must* have radial values specified using the DRV keyword");

        if (!(deck.hasKeyword<ParserKeywords::DZ>() || deck.hasKeyword<ParserKeywords::DZV>()) || !deck.hasKeyword<ParserKeywords::TOPS>())
            throw std::logic_error("The vertical cell size must be specified using the DZ or DZV, and the TOPS keywords");

        const std::vector<double>& drv     = deck.get<ParserKeywords::DRV>().back().getSIDoubleData();
        const std::vector<double>& dthetav = deck.get<ParserKeywords::DTHETAV>().back().getSIDoubleData();
        const std::vector<double>& tops    = deck.get<ParserKeywords::TOPS>().back().getSIDoubleData();
        OpmLog::info(fmt::format("\nCreating {} grid from keywords DRV, DTHETAV, DZV and TOPS", kind));

        if (drv.size() != this->getNX())
            throw std::invalid_argument("DRV keyword should have exactly " + std::to_string( this->getNX() ) + " elements");

        if (dthetav.size() != this->getNY())
            throw std::invalid_argument("DTHETAV keyword should have exactly " + std::to_string( this->getNY() ) + " elements");

        auto area = this->getNX() * this->getNY();
        std::size_t volume = this->getNX() * this->getNY() * this->getNZ();

        std::vector<double> dz(volume);
        if (deck.hasKeyword<ParserKeywords::DZ>()) {
            const std::vector<double>& dz_deck = deck.get<ParserKeywords::DZ>().back().getSIDoubleData();
            if (dz_deck.size() != volume)
                throw std::invalid_argument("DZ keyword should have exactly " + std::to_string( volume ) + " elements");
            dz = dz_deck;
        } else {
            const std::vector<double>& dzv = deck.get<ParserKeywords::DZV>().back().getSIDoubleData();
            if (dzv.size() != this->getNZ())
                throw std::invalid_argument("DZV keyword should have exactly " + std::to_string( this->getNZ() ) + " elements");
            for (std::size_t k= 0; k < this->getNZ(); k++)
                std::fill(dz.begin() + k*area, dz.begin() + (k+1)*area, dzv[k]);
        }

        if (tops.size() != (this->getNX() * this->getNY()))
            throw std::invalid_argument("TOPS keyword should have exactly " + std::to_string( this->getNX() * this->getNY() ) + " elements");

        {
            double total_angle = std::accumulate(dthetav.begin(), dthetav.end(), 0.0);

            if (std::abs( total_angle - 360 ) < 0.01)
                m_circle = deck.hasKeyword<ParserKeywords::CIRCLE>();
            else {
                if (total_angle > 360)
                    throw std::invalid_argument("More than 360 degrees rotation - cells will be double covered");
            }
        }

        /*
          Now the data has been validated, now we continue to create
          ZCORN and COORD vectors, and we are almost done.
        */
        {
            ZcornMapper zm( this->getNX(), this->getNY(), this->getNZ());
            CoordMapper cm(this->getNX(), this->getNY());
            std::vector<double> zcorn( zm.size() );
            std::vector<double> coord( cm.size() );
            {
                std::vector<double> depth(tops);
                for (std::size_t k = 0; k < this->getNZ(); k++) {
                    for (std::size_t j = 0; j < this->getNY(); j++) {
                        for (std::size_t i = 0; i < this->getNX(); i++) {
                            auto current_depth = depth[j * this->getNX() + i];
                            auto next_depth = current_depth + dz[k * area + j * this->getNX() + i];
                            for (std::size_t c=0; c < 4; c++) {
                                zcorn[ zm.index(i,j,k,c) ]     = current_depth;
                                zcorn[ zm.index(i,j,k,c + 4) ] = next_depth;
                            }
                            depth[j * this->getNX() + i] = next_depth;
                        }
                    }
                }
            }
            {
                std::vector<double> ri(this->getNX() + 1);
                std::vector<double> tj(this->getNY() + 1);
                double z1 = *std::min_element( zcorn.begin() , zcorn.end());
                double z2 = *std::max_element( zcorn.begin() , zcorn.end());
                ri[0] = deck.get<ParserKeywords::INRAD>().back().getRecord(0).getItem(0).getSIDouble( 0 );
                for (std::size_t i = 1; i <= this->getNX(); i++)
                    ri[i] = ri[i - 1] + drv[i - 1];

                tj[0] = 0;
                for (std::size_t j = 1; j <= this->getNY(); j++)
                    tj[j] = tj[j - 1] + dthetav[j - 1];


                for (std::size_t j = 0; j <= this->getNY(); j++) {
                    /*
                      The theta value is supposed to go counterclockwise, starting at 'twelve o clock'.
                    */
                    double t = M_PI * (90 - tj[j]) / 180;
                    double c = cos( t );
                    double s = sin( t );
                    for (std::size_t i = 0; i <= this->getNX(); i++) {
                        double r = ri[i];
                        double x = r*c;
                        double y = r*s;

                        coord[ cm.index(i,j,0,0) ] = x;
                        coord[ cm.index(i,j,1,0) ] = y;
                        coord[ cm.index(i,j,2,0) ] = z1;

                        coord[ cm.index(i,j,0,1) ] = x;
                        coord[ cm.index(i,j,1,1) ] = y;
                        coord[ cm.index(i,j,2,1) ] = z2;
                    }
                }

                // Save angles, used by the cylindrical grid to calculate volumes.
                if (is_cylindrical) {
                    m_rv = ri;
                    m_thetav = dthetav;
                }

            }
            initCornerPointGrid( coord, zcorn, nullptr);
        }
    }



    void EclipseGrid::initCornerPointGrid(const std::vector<double>& coord ,
                                          const std::vector<double>& zcorn ,
                                          const int * actnum)


    {
        m_coord = coord;
        m_zcorn = zcorn;

        m_input_coord = coord;
        m_input_zcorn = zcorn;

        ZcornMapper mapper( getNX(), getNY(), getNZ());
        zcorn_fixed = mapper.fixupZCORN( m_zcorn );
        this->resetACTNUM(actnum);
    }

    void EclipseGrid::initCornerPointGrid(const Deck& deck)
    {
        this->assertCornerPointKeywords(deck);

        OpmLog::info(fmt::format("\nCreating corner-point grid from "
                                 "keywords COORD, ZCORN and others"));

        const auto& coord = deck.get<ParserKeywords::COORD>().back();
        const auto& zcorn = deck.get<ParserKeywords::ZCORN>().back();

        this->initCornerPointGrid(coord.getSIDoubleData(),
                                  zcorn.getSIDoubleData(),
                                  nullptr);
    }

    bool EclipseGrid::hasCornerPointKeywords(const Deck& deck) {
        if (deck.hasKeyword<ParserKeywords::ZCORN>() && deck.hasKeyword<ParserKeywords::COORD>())
            return true;
        else
            return false;
    }


    void EclipseGrid::assertCornerPointKeywords(const Deck& deck)
    {
        const int nx = this->getNX();
        const int ny = this->getNY();
        const int nz = this->getNZ();
        {
            const auto& ZCORNKeyWord = deck.get<ParserKeywords::ZCORN>().back();

            if (ZCORNKeyWord.getDataSize() != static_cast<size_t>(8*nx*ny*nz)) {
                const std::string msg =
                    "Wrong size of the ZCORN keyword: Expected 8*nx*ny*nz = "
                    + std::to_string(static_cast<long long>(8*nx*ny*nz)) + " is "
                    + std::to_string(static_cast<long long>(ZCORNKeyWord.getDataSize()));
                OpmLog::error(msg);
                throw std::invalid_argument(msg);
            }
        }

        {
            const auto& COORDKeyWord = deck.get<ParserKeywords::COORD>().back();
            if (COORDKeyWord.getDataSize() != static_cast<size_t>(6*(nx + 1)*(ny + 1))) {
                const std::string msg =
                    "Wrong size of the COORD keyword: Expected 6*(nx + 1)*(ny + 1) = "
                    + std::to_string(static_cast<long long>(6*(nx + 1)*(ny + 1))) + " is "
                    + std::to_string(static_cast<long long>(COORDKeyWord.getDataSize()));
                OpmLog::error(msg);
                throw std::invalid_argument(msg);
            }
        }
    }


    bool EclipseGrid::hasGDFILE(const Deck& deck) {
        return deck.hasKeyword<ParserKeywords::GDFILE>();
    }


    bool EclipseGrid::hasCartesianKeywords(const Deck& deck) {
        if (hasDVDEPTHZKeywords( deck ))
            return true;
        else
            return hasDTOPSKeywords(deck);
    }


    bool EclipseGrid::hasRadialKeywords(const Deck& deck) {
        return deck.hasKeyword<ParserKeywords::RADIAL>() && hasCylindricalKeywords(deck);
    }


    bool EclipseGrid::hasSpiderKeywords(const Deck& deck) {
        return deck.hasKeyword<ParserKeywords::SPIDER>() && hasCylindricalKeywords(deck);
    }


    bool EclipseGrid::hasCylindricalKeywords(const Deck& deck) {
        if (deck.hasKeyword<ParserKeywords::INRAD>() &&
            deck.hasKeyword<ParserKeywords::TOPS>() &&
            (deck.hasKeyword<ParserKeywords::DZ>() || deck.hasKeyword<ParserKeywords::DZV>()) &&
            (deck.hasKeyword<ParserKeywords::DRV>() || deck.hasKeyword<ParserKeywords::DR>()) &&
            (deck.hasKeyword<ParserKeywords::DTHETA>() || deck.hasKeyword<ParserKeywords::DTHETAV>()))
            return true;
        else
            return false;
    }


    bool EclipseGrid::hasDVDEPTHZKeywords(const Deck& deck) {
        if (deck.hasKeyword<ParserKeywords::DXV>() &&
            deck.hasKeyword<ParserKeywords::DYV>() &&
            deck.hasKeyword<ParserKeywords::DZV>() &&
            deck.hasKeyword<ParserKeywords::DEPTHZ>())
            return true;
        else
            return false;
    }

    bool EclipseGrid::hasEqualDVDEPTHZ(const Deck& deck) {
        const std::vector<double>& DXV = deck.get<ParserKeywords::DXV>().back().getSIDoubleData();
        const std::vector<double>& DYV = deck.get<ParserKeywords::DYV>().back().getSIDoubleData();
        const std::vector<double>& DZV = deck.get<ParserKeywords::DZV>().back().getSIDoubleData();
        const std::vector<double>& DEPTHZ = deck.get<ParserKeywords::DEPTHZ>().back().getSIDoubleData();
        if (EclipseGrid::allEqual(DXV) &&
            EclipseGrid::allEqual(DYV)&&
            EclipseGrid::allEqual(DZV)&&
            EclipseGrid::allEqual(DEPTHZ))
            return true;
        else
            return false;
    }

    bool EclipseGrid::hasDTOPSKeywords(const Deck& deck) {
        if ((deck.hasKeyword<ParserKeywords::DX>() || deck.hasKeyword<ParserKeywords::DXV>()) &&
            (deck.hasKeyword<ParserKeywords::DY>() || deck.hasKeyword<ParserKeywords::DYV>()) &&
            (deck.hasKeyword<ParserKeywords::DZ>() || deck.hasKeyword<ParserKeywords::DZV>()) &&
            deck.hasKeyword<ParserKeywords::TOPS>())
            return true;
        else
            return false;
    }

    void EclipseGrid::assertVectorSize(const std::vector<double>& vector , std::size_t expectedSize , const std::string& vectorName) {
        if (vector.size() != expectedSize)
            throw std::invalid_argument("Wrong size for keyword: " + vectorName + ". Expected: " + std::to_string(expectedSize) + " got: " + std::to_string(vector.size()));
    }


    /*
      The body of the for loop in this method looks slightly
      peculiar. The situation is as follows:

        1. The EclipseGrid class will assemble the necessary keywords
           and create an ert ecl_grid instance.

        2. The ecl_grid instance will export ZCORN, COORD and ACTNUM
           data which will be used by the UnstructureGrid constructor
           in opm-core. If the ecl_grid is created with ZCORN as an
           input keyword that data is retained in the ecl_grid
           structure, otherwise the ZCORN data is created based on the
           internal cell geometries.

        3. When constructing the UnstructuredGrid structure strict
           numerical comparisons of ZCORN values are used to detect
           cells in contact, if all the elements the elements in the
           TOPS vector are specified[1] we will typically not get
           bitwise equality between the bottom of one cell and the top
           of the next.

       To remedy this we enforce bitwise equality with the
       construction:

          if (std::abs(nextValue - TOPS[targetIndex]) < z_tolerance)
             TOPS[targetIndex] = nextValue;

       [1]: This is of course assuming the intention is to construct a
            fully connected space covering grid - if that is indeed
            not the case the barriers must be thicker than 1e-6m to be
            retained.
    */

    std::vector<double> EclipseGrid::createTOPSVector(const std::array<int, 3>& dims,
            const std::vector<double>& DZ, const Deck& deck)
    {
        std::size_t volume = dims[0] * dims[1] * dims[2];
        std::size_t area = dims[0] * dims[1];
        const auto& TOPSKeyWord = deck.get<ParserKeywords::TOPS>().back();
        std::vector<double> TOPS = TOPSKeyWord.getSIDoubleData();

        if (TOPS.size() >= area) {
            std::size_t initialTOPSize = TOPS.size();
            TOPS.resize( volume );

            for (std::size_t targetIndex = area; targetIndex < volume; targetIndex++) {
                std::size_t sourceIndex = targetIndex - area;
                double nextValue = TOPS[sourceIndex] + DZ[sourceIndex];

                if (targetIndex >= initialTOPSize)
                    TOPS[targetIndex] = nextValue;
                else {
                    constexpr double z_tolerance = 1e-6;
                    if (std::abs(nextValue - TOPS[targetIndex]) < z_tolerance)
                        TOPS[targetIndex] = nextValue;
                }

            }
        }

        if (TOPS.size() != volume)
            throw std::invalid_argument("TOPS size mismatch");

        return TOPS;
    }

std::vector<double> EclipseGrid::createDVector(const std::array<int,3>& dims, std::size_t dim, const std::string& DKey,
            const std::string& DVKey, const Deck& deck)
    {
        std::size_t volume = dims[0] * dims[1] * dims[2];
        std::vector<double> D;
        if (deck.hasKeyword(DKey)) {
            D = deck[DKey].back().getSIDoubleData();

            const std::size_t area = dims[0] * dims[1];
            if (D.size() >= area && D.size() < volume) {
                /*
                  Only the top layer is required; for layers below the
                  top layer the value from the layer above is used.
                */
                std::size_t initialDSize = D.size();
                D.resize( volume );
                for (std::size_t targetIndex = initialDSize; targetIndex < volume; targetIndex++) {
                    std::size_t sourceIndex = targetIndex - area;
                    D[targetIndex] = D[sourceIndex];
                }
            }

            if (D.size() != volume)
                throw std::invalid_argument(DKey + " size mismatch");
        } else {
            const auto& DVKeyWord = deck[DVKey].back();
            const std::vector<double>& DV = DVKeyWord.getSIDoubleData();
            if (DV.size() != static_cast<std::size_t>(dims[dim]))
                throw std::invalid_argument(DVKey + " size mismatch");
            D.resize( volume );
            scatterDim( dims , dim , DV , D );
        }
        return D;
    }


    void EclipseGrid::scatterDim(const std::array<int, 3>& dims , std::size_t dim , const std::vector<double>& DV , std::vector<double>& D) {
        int index[3];
        for (index[2] = 0;  index[2] < dims[2]; index[2]++) {
            for (index[1] = 0; index[1] < dims[1]; index[1]++) {
                for (index[0] = 0;  index[0] < dims[0]; index[0]++) {
                    std::size_t globalIndex = index[2] * dims[1] * dims[0] + index[1] * dims[0] + index[0];
                    D[globalIndex] = DV[ index[dim] ];
                }
            }
        }
    }

    bool EclipseGrid::allEqual(const std::vector<double> &v) {
        auto comp = [](double x, double y) { return std::fabs(x - y) < 1e-12; };
        return std::equal(v.begin() + 1, v.end(), v.begin(), comp);
    }

    bool EclipseGrid::equal(const EclipseGrid& other) const {

        //double reltol = 1.0e-6;

        if (m_coord.size() != other.m_coord.size())
            return false;

        if (m_zcorn.size() != other.m_zcorn.size())
            return false;

        if (!(m_mapaxes == other.m_mapaxes))
            return false;

        if (m_actnum != other.m_actnum)
            return false;

        if (m_coord != other.m_coord)
            return false;

        if (m_zcorn != other.m_zcorn)
            return false;

        bool status = ((m_pinch == other.m_pinch)  && (m_minpvMode == other.getMinpvMode()));

        if (m_minpvMode != MinpvMode::Inactive) {
            status = status && (m_minpvVector == other.getMinpvVector());
        }

        return status;
    }

    std::size_t EclipseGrid::getNumActive( ) const {
        return this->m_nactive;
    }

    bool EclipseGrid::allActive( ) const {
        return (getNumActive() == getCartesianSize());
    }

    bool EclipseGrid::cellActive( std::size_t globalIndex ) const {
        assertGlobalIndex( globalIndex );
        if (m_actnum.empty()) {
            return true;
        } else {
            return m_actnum[globalIndex]>0;
        }
    }

    bool EclipseGrid::cellActive( std::size_t i , std::size_t j , std::size_t k ) const {
        assertIJK(i,j,k);

        std::size_t globalIndex = getGlobalIndex(i,j,k);
        return this->cellActive(globalIndex);
    }

    bool EclipseGrid::cellActiveAfterMINPV( std::size_t i , std::size_t j , std::size_t k, double cell_porv ) const {
        assertIJK(i,j,k);

        std::size_t globalIndex = getGlobalIndex(i,j,k);
        if (!this->cellActive(globalIndex))
            return false;

        return m_minpvMode == MinpvMode::Inactive || cell_porv >= m_minpvVector[globalIndex];
    }

    const std::vector<double>& EclipseGrid::activeVolume() const {
        if (!this->active_volume.has_value()) {
            std::vector<double> volume(this->m_nactive);

            #pragma omp parallel for schedule(static)
            for (std::int64_t active_index = 0; active_index < static_cast<std::int64_t>(this->m_active_to_global.size()); active_index++) {
                std::array<double,8> X;
                std::array<double,8> Y;
                std::array<double,8> Z;
                auto global_index = this->m_active_to_global[active_index];
                this->getCellCorners(global_index, X, Y, Z );
                if (m_rv && m_thetav) {
                    const auto[i,j,k] = this->getIJK(global_index);
                    const auto& r = *m_rv;
                    const auto& t = *m_thetav;
                    volume[active_index] = calculateCylindricalCellVol(r[i], r[i+1], t[j], Z[4] - Z[0]);
                } else
                    volume[active_index] = calculateCellVol(X, Y, Z);
            }

            this->active_volume = std::move(volume);
        }

        return this->active_volume.value();
    }


    double EclipseGrid::getCellVolume(std::size_t globalIndex) const {
        assertGlobalIndex( globalIndex );
        if (this->cellActive(globalIndex) && this->active_volume.has_value()) {
            auto active_index = this->activeIndex(globalIndex);
            return this->active_volume.value()[active_index];
        }

        std::array<double,8> X;
        std::array<double,8> Y;
        std::array<double,8> Z;
        this->getCellCorners(globalIndex, X, Y, Z );
        if (m_rv && m_thetav) {
            const auto[i,j,k] = this->getIJK(globalIndex);
            const auto& r = *m_rv;
            const auto& t = *m_thetav;
           return calculateCylindricalCellVol(r[i], r[i+1], t[j], Z[4] - Z[0]);
        } else {
            return calculateCellVol(X, Y, Z);
        }
    }

    double EclipseGrid::getCellVolume(std::size_t i , std::size_t j , std::size_t k) const {
        this->assertIJK(i,j,k);
        std::size_t globalIndex = getGlobalIndex(i,j,k);
        return this->getCellVolume(globalIndex);
    }

    double EclipseGrid::getCellThickness(std::size_t i , std::size_t j , std::size_t k) const {
        assertIJK(i,j,k);

        const std::array<int, 3> dims = getNXYZ();
        std::size_t globalIndex = i + j*dims[0] + k*dims[0]*dims[1];

        return getCellThickness(globalIndex);
    }

    double EclipseGrid::getCellThickness(std::size_t globalIndex) const {
        assertGlobalIndex( globalIndex );
        std::array<double,8> X;
        std::array<double,8> Y;
        std::array<double,8> Z;
        this->getCellCorners(globalIndex, X, Y, Z );

        double z2 = (Z[4]+Z[5]+Z[6]+Z[7])/4.0;
        double z1 = (Z[0]+Z[1]+Z[2]+Z[3])/4.0;
        double dz = z2-z1;
        return dz;
    }


    std::array<double, 3> EclipseGrid::getCellDims(std::size_t globalIndex) const {
        assertGlobalIndex( globalIndex );
        std::array<double,8> X;
        std::array<double,8> Y;
        std::array<double,8> Z;
        this->getCellCorners(globalIndex, X, Y, Z );

        // calculate dx
        double x1 = (X[0]+X[2]+X[4]+X[6])/4.0;
        double y1 = (Y[0]+Y[2]+Y[4]+Y[6])/4.0;
        double x2 = (X[1]+X[3]+X[5]+X[7])/4.0;
        double y2 = (Y[1]+Y[3]+Y[5]+Y[7])/4.0;
        double dx = sqrt(pow((x2-x1), 2.0) + pow((y2-y1), 2.0) );

        // calculate dy
        x1 = (X[0]+X[1]+X[4]+X[5])/4.0;
        y1 = (Y[0]+Y[1]+Y[4]+Y[5])/4.0;
        x2 = (X[2]+X[3]+X[6]+X[7])/4.0;
        y2 = (Y[2]+Y[3]+Y[6]+Y[7])/4.0;
        double dy = sqrt(pow((x2-x1), 2.0) + pow((y2-y1), 2.0));

        // calculate dz

        double z2 = (Z[4]+Z[5]+Z[6]+Z[7])/4.0;
        double z1 = (Z[0]+Z[1]+Z[2]+Z[3])/4.0;
        double dz = z2-z1;


        return std::array<double,3> {{dx, dy, dz}};
    }

    std::array<double, 3> EclipseGrid::getCellDims(std::size_t i , std::size_t j , std::size_t k) const {
        assertIJK(i,j,k);
        {
            std::size_t globalIndex = getGlobalIndex( i,j,k );

            return getCellDims(globalIndex);
        }
    }


    std::tuple<std::array<double, 3>,std::array<double, 3>,std::array<double, 3>>
    EclipseGrid::getCellAndBottomCenterNormal(std::size_t globalIndex) const {
        // bottom face is spanned by corners 4, ..., 7 and connects cells
        // with vertical index k and k+1
        assertGlobalIndex( globalIndex );
        std::array<double,8> X;
        std::array<double,8> Y;
        std::array<double,8> Z;
        this->getCellCorners(globalIndex, X, Y, Z );
        std::array<double,3> bottomCenter{{std::accumulate(X.begin() + 4, X.end(), 0.0) / 4.0,
                                           std::accumulate(Y.begin() + 4, Y.end(), 0.0) / 4.0,
                                           std::accumulate(Z.begin() + 4, Z.end(), 0.0) / 4.0}};

        // Calculate normal scaled with area via the triangles spanned by the center and two neighboring
        // corners.
        std::array<double,3> bottomFaceNormal{};
        // Reorder counter-clock.wise
        std::array<std::size_t,4> bottomIndices = {4, 5, 7, 6};
        std::array<double,3> oldCorner{ *(X.begin() + 6),
                                        *(Y.begin() + 6),
                                        *(Z.begin() + 6)};

        for(const auto cornerIndex: bottomIndices)
        {
            std::array<double,3> newCorner{ *(X.begin() + cornerIndex), *(Y.begin() + cornerIndex),
                                            *(Z.begin() + cornerIndex) };
            std::array<double,3> side1, side2, normalTriangle;
            std::transform(oldCorner.begin(), oldCorner.end(), bottomCenter.begin(),
                           side1.begin(), std::minus<double>());
            std::transform(newCorner.begin(), newCorner.end(), bottomCenter.begin(),
                           side2.begin(), std::minus<double>());

            cross(side1, side2, normalTriangle);
            // use std::plus to make normal point outwards (Z-axis points downwards)
            std::transform(bottomFaceNormal.begin(), bottomFaceNormal.end(),
                           normalTriangle.begin(), bottomFaceNormal.begin(),
                           std::plus<double>());
            oldCorner = newCorner;
        }

        std::transform(bottomFaceNormal.begin(), bottomFaceNormal.end(),
                       bottomFaceNormal.begin(), [](const double& a){ return 0.5 * a;});
        std::array<double,3> cellCenter = { std::accumulate(X.begin(), X.end(), 0.0) / 8.0,
                                             std::accumulate(Y.begin(), Y.end(), 0.0) / 8.0,
                                             std::accumulate(Z.begin(), Z.end(), 0.0) / 8.0 };
            return std::make_tuple(cellCenter,
                                   bottomCenter,
                                   bottomFaceNormal);
    }

    std::array<double, 3> EclipseGrid::getCellCenter(std::size_t globalIndex) const {
        assertGlobalIndex( globalIndex );
        std::array<double,8> X;
        std::array<double,8> Y;
        std::array<double,8> Z;
        this->getCellCorners(globalIndex, X, Y, Z );
        return std::array<double,3> { { std::accumulate(X.begin(), X.end(), 0.0) / 8.0,
                                        std::accumulate(Y.begin(), Y.end(), 0.0) / 8.0,
                                        std::accumulate(Z.begin(), Z.end(), 0.0) / 8.0 } };
    }


    std::array<double, 3> EclipseGrid::getCellCenter(std::size_t i,std::size_t j, std::size_t k) const {
        assertIJK(i,j,k);

        const std::array<int, 3> dims = getNXYZ();

        std::size_t globalIndex = i + j*dims[0] + k*dims[0]*dims[1];

        return getCellCenter(globalIndex);
    }

    /*
      This is the numbering of the corners in the cell.

       bottom                           j
         6---7                        /|\
         |   |                         |
         4---5                         |
                                       |
       top                             o---------->  i
         2---3
         |   |
         0---1

     */

    std::array<double, 3> EclipseGrid::getCornerPos(std::size_t i, std::size_t j, std::size_t k, std::size_t corner_index) const {
        assertIJK(i,j,k);
        if (corner_index >= 8)
            throw std::invalid_argument("Invalid corner position");
        {
            const std::array<int, 3> dims = getNXYZ();

            std::array<double,8> X = {0.0};
            std::array<double,8> Y = {0.0};
            std::array<double,8> Z = {0.0};

            std::array<int, 3> ijk;

            ijk[0] = i;
            ijk[1] = j;
            ijk[2] = k;

            getCellCorners(ijk, dims, X, Y, Z );

            return std::array<double, 3> {{X[corner_index], Y[corner_index], Z[corner_index]}};
        }
    }

    bool EclipseGrid::isValidCellGeomtry(const std::size_t globalIndex,
                                         const UnitSystem& usys) const
    {
        const auto threshold = usys.to_si(UnitSystem::measure::length, 1.0e+20f);

        auto is_finite = [threshold](const double coord_component)
        {
            return std::abs(coord_component) < threshold;
        };

        std::array<double,8> X, Y, Z;
        this->getCellCorners(globalIndex, X, Y, Z);

        const auto finite_coord = std::all_of(X.begin(), X.end(), is_finite)
            && std::all_of(Y.begin(), Y.end(), is_finite)
            && std::all_of(Z.begin(), Z.end(), is_finite);

        if (! finite_coord) {
            return false;
        }

        const auto max_pillar_point_distance = std::max({
            Z[0 + 4] - Z[0 + 0],
            Z[1 + 4] - Z[1 + 0],
            Z[2 + 4] - Z[2 + 0],
            Z[3 + 4] - Z[3 + 0],
        });

        // Define points as "well separated" if maximum distance exceeds
        // 1e-4 length units (e.g., 0.1 mm in METRIC units).  May consider
        // using a coarser tolerance/threshold here.
        return max_pillar_point_distance
            >  usys.to_si(UnitSystem::measure::length, 1.0e-4);
    }

    double EclipseGrid::getCellDepth(std::size_t globalIndex) const {
        assertGlobalIndex( globalIndex );

        auto it = this->m_aquifer_cell_depths.find(globalIndex);
        return it != this->m_aquifer_cell_depths.end() ? it->second : computeCellGeometricDepth(globalIndex);
    }

    double EclipseGrid::computeCellGeometricDepth(std::size_t globalIndex) const {
        std::array<double,8> X;
        std::array<double,8> Y;
        std::array<double,8> Z;
        this->getCellCorners(globalIndex, X, Y, Z );

        double z2 = (Z[4]+Z[5]+Z[6]+Z[7])/4.0;
        double z1 = (Z[0]+Z[1]+Z[2]+Z[3])/4.0;
        return (z1 + z2)/2.0;
    }

    double EclipseGrid::getCellDepth(std::size_t i, std::size_t j, std::size_t k) const {
        this->assertIJK(i,j,k);
        std::size_t globalIndex = getGlobalIndex(i,j,k);
        return this->getCellDepth(globalIndex);
    }

    const std::map<std::size_t, std::array<int,2>>& EclipseGrid::getAquiferCellTabnums() const {
        return m_aquifer_cell_tabnums;
    }

    const std::vector<int>& EclipseGrid::getACTNUM( ) const {

        return m_actnum;
    }

    const std::vector<double>& EclipseGrid::getCOORD() const {

        return m_coord;
    }

    std::size_t EclipseGrid::fixupZCORN() {

        ZcornMapper mapper( getNX(), getNY(), getNZ());

        return mapper.fixupZCORN( m_zcorn );
    }

    const std::vector<double>& EclipseGrid::getZCORN( ) const {

        return m_zcorn;
    }

    void EclipseGrid::save(const std::string& filename, bool formatted, const std::vector<Opm::NNCdata>& nnc, const Opm::UnitSystem& units) const {

        Opm::UnitSystem::UnitType unitSystemType = units.getType();
        constexpr auto length = ::Opm::UnitSystem::measure::length;

        const std::array<int, 3> dims = getNXYZ();

        // Preparing vectors to be saved

        // create coord vector of floats with input units, converted from SI
        std::vector<float> coord_f;
        coord_f.resize(m_coord.size());

        auto convert_length = [&units](const double x) { return static_cast<float>(units.from_si(length, x)); };

        if (m_input_coord.has_value()) {
            std::transform(m_input_coord.value().begin(), m_input_coord.value().end(), coord_f.begin(), convert_length);
        } else {
            std::transform(m_coord.begin(), m_coord.end(), coord_f.begin(), convert_length);
        }

        // create zcorn vector of floats with input units, converted from SI
        std::vector<float> zcorn_f;
        zcorn_f.resize(m_zcorn.size());

        if (m_input_zcorn.has_value()) {
            std::transform(m_input_zcorn.value().begin(), m_input_zcorn.value().end(), zcorn_f.begin(), convert_length);
        } else {
            std::transform(m_zcorn.begin(), m_zcorn.end(), zcorn_f.begin(), convert_length);
        }

        m_input_coord.reset();
        m_input_zcorn.reset();

        std::vector<int> filehead(100,0);
        filehead[0] = 3;                     // version number
        filehead[1] = 2007;                  // release year
        filehead[6] = 1;                     // corner point grid

        std::vector<int> gridhead(100,0);
        gridhead[0] = 1;                    // corner point grid
        gridhead[1] = dims[0];              // nI
        gridhead[2] = dims[1];              // nJ
        gridhead[3] = dims[2];              // nK
        gridhead[24] = 1;                   // NUMRES (number of reservoirs)
        //gridhead[25] = 1;                 // TODO: This value depends on LGRs?

        std::vector<int> nnchead(10, 0);
        std::vector<int> nnc1;
        std::vector<int> nnc2;

        for (const NNCdata& n : nnc ) {
            nnc1.push_back(n.cell1 + 1);
            nnc2.push_back(n.cell2 + 1);
        }

        nnchead[0] = nnc1.size();

        std::vector<std::string> gridunits;

        switch (unitSystemType) {
        case Opm::UnitSystem::UnitType::UNIT_TYPE_METRIC:
            gridunits.push_back("METRES");
            break;
        case Opm::UnitSystem::UnitType::UNIT_TYPE_FIELD:
            gridunits.push_back("FEET");
            break;
        case Opm::UnitSystem::UnitType::UNIT_TYPE_LAB:
            gridunits.push_back("CM");
            break;
        default:
            OPM_THROW(std::runtime_error, "Unit system not supported when writing to EGRID file");
            break;
        }
        gridunits.push_back("");

        std::vector<int> endgrid = {};

        // Writing vectors to egrid file

        Opm::EclIO::EclOutput egridfile(filename, formatted);
        egridfile.write("FILEHEAD", filehead);

        if (this->m_mapaxes.has_value()) {
            const auto& mapunits = this->m_mapaxes.value().mapunits();
            if (mapunits.has_value())
                egridfile.write("MAPUNITS", std::vector<std::string>{ mapunits.value() });

            egridfile.write("MAPAXES", this->m_mapaxes.value().input());
        }

        egridfile.write("GRIDUNIT", gridunits);
        egridfile.write("GRIDHEAD", gridhead);

        egridfile.write("COORD", coord_f);
        egridfile.write("ZCORN", zcorn_f);

        egridfile.write("ACTNUM", m_actnum);
        egridfile.write("ENDGRID", endgrid);
        // nnc for LGR is not currently supported.
        for (std::size_t index : m_print_order_lgr_cells) {
            lgr_children_cells[index].save(egridfile, nnc, units);
        }

        if (nnc1.size() > 0){
            egridfile.write("NNCHEAD", nnchead);
            egridfile.write("NNC1", nnc1);
            egridfile.write("NNC2", nnc2);
            // Placeholder the method that handles  LGR nnc
            // nnc in LGR is currently not supportd
            // for (const EclipseGridLGR& lgr_cell : lgr_children_cells) {
            //     lgr_cell.save_nnc(egridfile);
            // }

        }
    }

    EclipseGridLGR& EclipseGrid::getLGRCell(std::size_t index){
        return lgr_children_cells[index];
      };

    EclipseGridLGR& EclipseGrid::getLGRCell(const std::string& lgr_tag) const
    {
        std::optional<std::reference_wrapper<EclipseGridLGR>>lgr_found = std::nullopt;
        for (auto& lgr_cell : lgr_children_cells) 
        {
            auto lgr_temp = lgr_cell.get_child_LGR_cell(lgr_tag);          
            if (lgr_temp != std::nullopt)
                lgr_found  = lgr_temp;
        }
        if (lgr_found == std::nullopt)
        {
            throw std::runtime_error("No EclipseGridLGR found with tag: " + lgr_tag);
        }
        else
        {
            return lgr_found.value().get();
        }        
    }

    int EclipseGrid::getLGR_global_father(std::size_t global_index,  const std::string& lgr_tag) const
    {
        const EclipseGridLGR& lgr_cell = getLGRCell(lgr_tag);
        const std::string& father_label  = lgr_cell.get_father_label(); 
        if (father_label == "GLOBAL")
        {
            return lgr_cell.get_hostnum(global_index);
        }
        else
        {
            return getLGR_global_father(lgr_cell.get_hostnum(global_index), father_label);
        }        

    }

    /**
    * @brief Computes the dimensions of a local grid refinement (LGR) cell.
    * 
    * This function calculates the dimensions of an LGR cell by dividing 
    * the parent cell's dimensions by the subdivision ratio for the given 
    * LGR tag.
    * 
    * @param i Local index of the parent cell in the x-direction.
    * @param j Local index of the parent cell in the y-direction.
    * @param k Local index of the parent cell in the z-direction.
    * @param lgr_tag Identifier for the LGR region.
    * @return std::array<double,3> The computed dimensions of the LGR cell 
    *         in the x, y, and z directions.
    */
    std::array<double,3> EclipseGrid::getCellDimensionsLGR(size_t i, size_t j, size_t k, const std::string& lgr_tag) const
    {
        std::array<int,3> subdivision = getCellSubdivisionRatioLGR(lgr_tag);
        std::array<double,3> cell_dims = getCellDims(i, j, k);
        
        return {cell_dims[0] / static_cast<double>(subdivision[0]), 
                cell_dims[1] / static_cast<double>(subdivision[1]), 
                cell_dims[2] / static_cast<double>(subdivision[2])};
    }

    double EclipseGrid::getCellDepthLGR(size_t i, size_t j, size_t k, const std::string& lgr_tag) const
    {
        auto split_domain = [](double start_val, double end_val, int num_points, int index) {
            if (index < 0 || index >= num_points) 
            {
                throw std::out_of_range("Index must be between 1 and num_points.");
            }
            double step = (end_val - start_val) / (num_points - 1);
            return start_val + (index * step);
        };
        
        auto refine_cell = [&split_domain](std::array<double, 8>& Z, int nz, int k_pos) {
            std::array<double, 8> Zdiv;        
            for (int index = 0; index < 4; index++) 
            {
                Zdiv[index] = split_domain(Z[index], Z[index + 4], nz + 1, k_pos);
                Zdiv[index + 4] = split_domain(Z[index], Z[index + 4], nz + 1, k_pos + 1);
            }        
            Z = Zdiv;
        };

        std::vector<std::reference_wrapper<const std::string>> label_list;
        std::vector<std::size_t> recurrent_global_index;   
        std::size_t fater_global_index; 

        std::array<double,8> X;
        std::array<double,8> Y;
        std::array<double,8> Z;

        getLGRCell(lgr_tag).get_label_child_to_top_father(label_list);
        getLGRCell(lgr_tag).get_global_index_child_to_top_father(recurrent_global_index, i,j,k);
        fater_global_index = recurrent_global_index.back();
        recurrent_global_index.pop_back();        

        this->getCellCorners(fater_global_index, X, Y, Z );
        for (size_t index = 0; index < label_list.size(); index++)
        {
            int nz = getLGRCell(label_list[index].get()).getNZ();
            int k_pos = getLGRCell(label_list[index].get()).getIJK(recurrent_global_index[index])[2];
            refine_cell(Z, nz, k_pos);
        }
        double z1 = (Z[0]+Z[1]+Z[2]+Z[3])/4.0;
        double z2 = (Z[4]+Z[5]+Z[6]+Z[7])/4.0;
        return (z1 + z2) / 2.0;
    }


    std::array<int,3> EclipseGrid::getCellSubdivisionRatioLGR(const std::string& lgr_tag,
                                                                    std::array<int,3> acum ) const
    {
        const EclipseGridLGR& lgr_cell = getLGRCell(lgr_tag);
        const std::string& father_label  = lgr_cell.get_father_label(); 
        acum = {static_cast<int>(acum[0]*lgr_cell.getNX()), 
                static_cast<int>(acum[1]*lgr_cell.getNY()), 
                static_cast<int>(acum[2]*lgr_cell.getNZ())};
        if (father_label == "GLOBAL")
        {
            return acum;
        }
        else
        {            
            return getCellSubdivisionRatioLGR(father_label, acum);
        }    
    }


    std::vector<GridDims> EclipseGrid::get_lgr_children_gridim() const
    {
        std::vector<GridDims> lgr_children_gridim;
        for (std::string lgr_tag : get_all_lgr_labels()) 
        {
            EclipseGridLGR& lgr_cell =  getLGRCell(lgr_tag);
            lgr_children_gridim.emplace_back(lgr_cell.getNX(), lgr_cell.getNY(), lgr_cell.getNZ());            
        }
        return lgr_children_gridim;
    }

    void EclipseGrid::set_lgr_refinement(const std::string& lgr_tag, const std::vector<double>& coords, 
                                                                            const std::vector<double> & zcorn){              
        for (auto& lgr_cell : lgr_children_cells) {
            lgr_cell.set_lgr_refinement(lgr_tag, coords, zcorn);
        }                                                                             
    }


    void EclipseGrid::init_children_host_cells(bool logical){
        if (logical)
            init_children_host_cells_logical();
        else
            init_children_host_cells_geometrical();
    }

    void EclipseGrid::init_children_host_cells_geometrical(){
        auto getAllCellCorners = [this](const auto& father_list){
            std::vector<std::array<double, 8>> X(father_list.size());
            std::vector<std::array<double, 8>> Y(father_list.size());
            std::vector<std::array<double, 8>> Z(father_list.size());
            for (std::size_t index = 0; index < father_list.size(); index++) {
                getCellCorners(father_list[index],X[index],Y[index],Z[index]);
            }
            return std::make_tuple(X, Y,Z);
        };

        std::vector<double> element_centerX, element_centerY, element_centerZ;
        for (EclipseGridLGR& lgr_cell : lgr_children_cells) {
            std::tie(element_centerX, element_centerY,element_centerZ) =
                VectorUtil::callMethodForEachInputOnObjectXYZ<EclipseGridLGR, std::array<double,3>, int,std::array<double, 3> (EclipseGridLGR::*)(size_t) const>
                (lgr_cell, &EclipseGridLGR::getCellCenter, lgr_cell.getActiveMap());
            auto [host_cellX, host_cellY, host_cellZ]  =  getAllCellCorners(lgr_cell.get_father_global());
            auto inside_el = GeometryUtil::isInsideElement(element_centerX, element_centerY, element_centerZ, host_cellX, host_cellY, host_cellZ);
            std::vector<int> host_cells_global_ref = VectorUtil::filterArray<int>(lgr_cell.get_father_global(), inside_el);
            lgr_cell.set_hostnum(host_cells_global_ref);
            lgr_cell.init_children_host_cells();
        }
    }


    void EclipseGrid::init_children_host_cells_logical(){       
        auto  IJK_location = [](const std::size_t&  nx, const std::size_t& ny,const std::size_t& nz, 
                                          const std::size_t& host_nx, const std::size_t& host_ny, const std::size_t& host_nz,
                                          const std::size_t& base_host_nx, const std::size_t& base_host_ny, const std::size_t& base_host_nz){
            auto [i_list, j_list, k_list] = VectorUtil::generate_cartesian_product(0, nx-1, 0, ny-1, 0, nz-1);
            std::vector<std::size_t> resultI = VectorUtil::scalarVectorOperation(nx/host_nx, i_list,  std::divides<std::size_t>{});
            resultI = VectorUtil::vectorScalarOperation(resultI, base_host_nx, std::plus<std::size_t>{});
            std::vector<std::size_t> resultJ = VectorUtil::scalarVectorOperation(ny/host_ny, j_list,  std::divides<std::size_t>{});
            resultJ = VectorUtil::vectorScalarOperation(resultJ, base_host_ny , std::plus<std::size_t>{});
            std::vector<std::size_t> resultK = VectorUtil::scalarVectorOperation(nz/host_nz, k_list,  std::divides<std::size_t>{});
            resultK = VectorUtil::vectorScalarOperation(resultK, base_host_nz, std::plus<std::size_t>{});
            return std::make_tuple(resultI, resultJ, resultK);
        }; 
        auto getAllGlobalIndex  = [this](const std::vector<std::size_t>& i_list, 
                                                   const std::vector<std::size_t>& j_list, 
                                                   const std::vector<std::size_t>& k_list){
            std::vector<int> global_index(i_list.size());
            for (std::size_t idx = 0; idx < i_list.size(); ++idx) {
                global_index[idx] = getGlobalIndex(i_list[idx], j_list[idx], k_list[idx]);
            }
            return global_index;
        };
        for (EclipseGridLGR& lgr_cell : lgr_children_cells) {
            const std::array<int,3>& host_low_fatherIJK = lgr_cell.get_low_fatherIJK();            
            const std::array<int,3>& host_up_fatherIJK = lgr_cell.get_up_fatherIJK();    
            const std::array<int,3> host_IJK = {host_up_fatherIJK[0] - host_low_fatherIJK[0] + 1,
                                                host_up_fatherIJK[1] - host_low_fatherIJK[1] + 1,
                                                host_up_fatherIJK[2] - host_low_fatherIJK[2] + 1};
            auto [i_list, j_list, k_list] = IJK_location(lgr_cell.getNX(), lgr_cell.getNY(), lgr_cell.getNZ(), 
                          static_cast<std::size_t>(host_IJK[0]), static_cast<std::size_t>(host_IJK[1]), static_cast<std::size_t>(host_IJK[2]),
                            host_low_fatherIJK[0], host_low_fatherIJK[1], host_low_fatherIJK[2]);
            std::vector<int> host_cells_global_ref = getAllGlobalIndex(i_list, j_list, k_list);   
            lgr_cell.set_hostnum(host_cells_global_ref);
            lgr_cell.init_children_host_cells();
        }
    }

    const std::vector<int>& EclipseGrid::getActiveMap() const {
        return m_active_to_global;
    }

    void EclipseGrid::assertLabelLGR(const std::string& label) const{
        if (std::find(all_lgr_labels.begin(), all_lgr_labels.end(), label) == all_lgr_labels.end()){
            throw std::invalid_argument("LGR label not found");
        }
    }


    void EclipseGrid::assertIndexLGR(std::size_t localIndex) const {
     if (std::find(lgr_active_index.begin(), lgr_active_index.end(), localIndex) != lgr_active_index.end()) {
        throw std::invalid_argument("input provided is an LGR refined cell");
      }
    }

    std::size_t EclipseGrid::getActiveIndexLGR(const std::string& label, std::size_t localIndex) const {
        assertLabelLGR(label);
        return activeIndexLGR(label, localIndex);
    }

    std::size_t EclipseGrid::getActiveIndexLGR(const std::string& label, std::size_t i, std::size_t j, std::size_t k) const {
        assertLabelLGR(label);
        return activeIndexLGR(label,i,j,k);
    }

    std::size_t EclipseGrid::activeIndexLGR(const std::string& label, std::size_t localIndex) const {
        if (lgr_label == label){
            std::size_t local_global_ind = getActiveIndex(localIndex);
            this->assertIndexLGR(local_global_ind);
            return lgr_level_active_map[local_global_ind] + lgr_global_counter;
        }
        else if (lgr_children_cells.empty()) {
            return 0;
        }
        else {
            return std::accumulate(lgr_children_cells.begin(), lgr_children_cells.end(), std::size_t{0},
                                   [&label,localIndex](std::size_t val, const auto& lgr_cell) {
                                       return val + lgr_cell.activeIndexLGR(label, localIndex);
                                   });
        }
    }


    std::size_t EclipseGrid::activeIndexLGR(const std::string& label, std::size_t i, std::size_t j, std::size_t k) const {
        if (lgr_label == label) {
            this->assertIJK(i,j,k);
            std::size_t local_global_ind = getActiveIndex(i,j,k);
            this->assertIndexLGR(local_global_ind);
            return lgr_level_active_map[local_global_ind] + lgr_global_counter;
        }
        else if (lgr_children_cells.empty()) {
            return 0;
        }
        else {
            return std::accumulate(lgr_children_cells.begin(), lgr_children_cells.end(), std::size_t{0},
                                   [&label,i,j,k](std::size_t val, const auto& lgr_cell) {
                                       return val + lgr_cell.activeIndexLGR(label, i, j, k);
                                   });
        }
    }


    std::size_t EclipseGrid::getTotalActiveLGR() const
    {
        std::size_t num_coarse_level =  getNumActive();
        std::size_t num_fine_cells = 0;
        std::size_t num_spanned_cells = 0;
        for (const auto&  fine_cells : lgr_children_cells) {
            num_fine_cells += fine_cells.getTotalActiveLGR();
            num_spanned_cells += fine_cells.get_father_global().size();

        }
        return num_coarse_level + num_fine_cells - num_spanned_cells;
    }


    void EclipseGrid::init_lgr_cells(const LgrCollection& lgr_input) {
        save_all_lgr_labels(lgr_input);
        create_lgr_cells_tree(lgr_input);
        //initialize LGRObjbect Indices
        initializeLGRObjectIndices(0);
        //parse fatherLGRObjbect Indices to children
        propagateParentIndicesToLGRChildren(0);
        // initialize the LGR tree indices for each refined cell.
        initializeLGRTreeIndices();
        // parse the reference indices to object in the global level.
        parseGlobalReferenceToChildren();
        // initialize the host cells for each LGR cell.
        // because the standard algorithm is based on topological information
        // it does not need for the refinement information to be parsed.
        init_children_host_cells();
    }

    void EclipseGrid::propagateParentIndicesToLGRChildren(int index){
        lgr_level_father = index;
        for (auto& cell : lgr_children_cells) {
            cell.propagateParentIndicesToLGRChildren(lgr_level);
        }
    }

    int EclipseGrid::initializeLGRObjectIndices(int num){
        lgr_level = num;
        num++;
        std::accumulate(m_print_order_lgr_cells.begin(),
                        m_print_order_lgr_cells.end(), num,
                        [this](int n, const auto& it)
                        {
                            return lgr_children_cells[it].initializeLGRObjectIndices(n);
                        });
        return num;
    }

    void EclipseGrid::save_all_lgr_labels(const LgrCollection& lgr_input) {
        all_lgr_labels.reserve(lgr_input.size()+1);
        all_lgr_labels.push_back("GLOBAL");
        for (std::size_t index = 0; index < lgr_input.size(); index++) {
            all_lgr_labels.push_back(lgr_input.getLgr(index).NAME());
        }
    }

    void EclipseGrid::create_lgr_cells_tree(const LgrCollection& lgr_input) {
          auto IJK_global = [this](const auto& i_list, const auto& j_list, const auto& k_list){
            if (!(i_list.size() == j_list.size()) && (j_list.size() == k_list.size()) ){
                 throw std::invalid_argument("Sizes are not compatible.");
            }
            std::vector<std::size_t> global_ind_active(i_list.size());
            for (std::size_t index = 0; index < i_list.size(); index++) {
                global_ind_active[index] = this->getActiveIndex(i_list[index],j_list[index],k_list[index]);
            }
            return global_ind_active;
        };
         for (std::size_t index = 0; index < lgr_input.size(); index++) {
            const auto& lgr_cell = lgr_input.getLgr(index);
            if (this->lgr_label == lgr_cell.PARENT_NAME()){
                lgr_grid = true;
                // auto [i_list, j_list, k_list] = lgr_cell.parent_cellsIJK();
                auto [i_list, j_list, k_list] = VectorUtil::generate_cartesian_product(lgr_cell.I1(), lgr_cell.I2(),
                                                                                                                           lgr_cell.J1(), lgr_cell.J2(),
                                                                                                                           lgr_cell.K1(), lgr_cell.K2());

                auto father_lgr_index = IJK_global(i_list, j_list, k_list);

                std::array<int,3> lowIJK = {lgr_cell.I1(), lgr_cell.J1(),lgr_cell.K1()};
                std::array<int,3> upIJK  = {lgr_cell.I2(), lgr_cell.J2(),lgr_cell.K2()};

                lgr_children_cells.emplace_back(lgr_cell.NAME(), this->lgr_label,
                                                lgr_cell.NX(), lgr_cell.NY(), lgr_cell.NZ(), father_lgr_index,
                                                lowIJK,upIJK);

                lgr_children_cells.back().create_lgr_cells_tree(lgr_input);
            }
        }
        EclipseGridLGR::vec_size_t father_label_sorting(lgr_children_cells.size(),0);
        m_print_order_lgr_cells.resize(lgr_children_cells.size());
        std::iota(m_print_order_lgr_cells.begin(), m_print_order_lgr_cells.end(), 0); //
        std::transform(lgr_children_cells.begin(), lgr_children_cells.end(), father_label_sorting.begin(),
                       [](const auto& cell){return cell.get_father_global()[0];});

        std::sort(m_print_order_lgr_cells.begin(), m_print_order_lgr_cells.end(), [&](std::size_t i1, std::size_t i2) {
            return father_label_sorting[i1] < father_label_sorting[i2]; //
        });

        std::sort(lgr_children_cells.begin(), lgr_children_cells.end(),[](const EclipseGridLGR& a, const EclipseGridLGR& b) {
                return a.get_father_global()[0] < b.get_father_global()[0]; //
        });

        lgr_children_labels.reserve(lgr_children_cells.size());
        std::transform(lgr_children_cells.begin(),
                       lgr_children_cells.end(),
                       std::back_inserter(lgr_children_labels),
                       [](const auto& it) { return it.lgr_label; });
        lgr_active_index.resize(lgr_children_cells.size(),0);
    }

    void EclipseGrid::initializeLGRTreeIndices(){
        // initialize the LGR tree indices for each refined cell.
        auto set_map_scalar = [&](const auto& vec, const auto& value){
             num_lgr_children_cells[vec] = value;
        };
        auto set_vec_value = []( auto& vec, const auto& num_ref, const auto& data){
            for (std::size_t index = 0; index < num_ref.size(); index++) {
                vec[num_ref[index]] = data;
            }
        };
        std::vector<std::size_t> lgr_level_numbering_counting(getNumActive(),1);
        lgr_level_active_map.resize(getNumActive(),0);
        for (const auto& cell:lgr_children_cells) {
            set_map_scalar(cell.getFatherGlobalID(), cell.getTotalActiveLGR());
        }
        std::vector<std::size_t> bottom_lgr_cells;
        std::size_t index = 0;
        for (const auto& [key, value] : num_lgr_children_cells) {
            const std::size_t head_lgr_cell = key[0];
            lgr_level_numbering_counting[head_lgr_cell] = value;
            lgr_active_index[index] = head_lgr_cell;
            bottom_lgr_cells.insert(bottom_lgr_cells.end(), key.begin() + 1, key.end());
            set_vec_value(lgr_level_numbering_counting, bottom_lgr_cells, 0);
            index++;
        }
        std::partial_sum(lgr_level_numbering_counting.begin(), lgr_level_numbering_counting.end(),
                         lgr_level_active_map.begin());
        lgr_level_active_map.reserve(lgr_level_active_map.size()+1);
        lgr_level_active_map.insert(lgr_level_active_map.begin(),0);
        for (auto& lgr_cell : lgr_children_cells) {
            lgr_cell.initializeLGRTreeIndices();
        }
    }

    void EclipseGrid::parseGlobalReferenceToChildren(){
       for (std::size_t index = 0; index < lgr_children_cells.size(); index++)
        {
            lgr_children_cells[index].set_lgr_global_counter(lgr_level_active_map[lgr_active_index[index]] +
                                                             this->lgr_global_counter);
            lgr_children_cells[index].parseGlobalReferenceToChildren();
        }
    }


    void EclipseGrid::resetACTNUM() {
        std::size_t global_size = this->getCartesianSize();
        this->m_actnum.assign(global_size, 1);
        this->m_nactive = global_size;

        this->m_global_to_active.resize(global_size);
        std::iota(this->m_global_to_active.begin(), this->m_global_to_active.end(), 0);
        this->m_active_to_global = this->m_global_to_active;
        this->active_volume = std::nullopt;
    }

    void EclipseGrid::resetACTNUM(const int* actnum) {
        if (actnum == nullptr)
            this->resetACTNUM();
        else {
            auto global_size = this->getCartesianSize();
            this->m_global_to_active.clear();
            this->m_active_to_global.clear();
            this->m_actnum.resize(global_size);
            this->m_nactive = 0;

            for (std::size_t n = 0; n < global_size; n++) {
                this->m_actnum[n] = actnum[n];
                // numerical aquifer cells need to be active
                if (this->m_aquifer_cells.count(n) > 0) {
                    this->m_actnum[n] = 1;
                }
                if (this->m_actnum[n] > 0) {
                    this->m_global_to_active.push_back(this->m_nactive);
                    this->m_active_to_global.push_back(n);
                    this->m_nactive++;
                } else {
                    this->m_global_to_active.push_back(-1);

                }
            }
            this->active_volume = std::nullopt;
        }
    }


    void EclipseGrid::setMINPVV(const std::vector<double>& minpvv) {
        if (m_minpvMode == MinpvMode::Inactive || m_minpvMode == MinpvMode::EclSTD) {
            if (minpvv.size() != this->getCartesianSize()) {
                throw std::runtime_error("EclipseGrid::setMINPVV(): MINPVV vector size differs from logical cartesian size of grid.");
            }
            m_minpvVector = minpvv;
            m_minpvMode = MinpvMode::EclSTD;
        }
    }
    void EclipseGrid::resetACTNUM(const std::vector<int>& actnum) {
        if (actnum.size() != getCartesianSize())
            throw std::runtime_error("resetACTNUM(): actnum vector size differs from logical cartesian size of grid.");

        this->resetACTNUM(actnum.data());
    }

    ZcornMapper EclipseGrid::zcornMapper() const {
        return ZcornMapper( getNX() , getNY(), getNZ() );
    }

    void EclipseGrid::updateNumericalAquiferCells(const Deck& deck) {
        using AQUNUM = ParserKeywords::AQUNUM;
        if ( !deck.hasKeyword<AQUNUM>() ) {
            return;
        }
        const auto &aqunum_keywords = deck.getKeywordList<AQUNUM>();
        for (const auto &keyword : aqunum_keywords) {
            for (const auto &record : *keyword) {
                const std::size_t i = record.getItem<AQUNUM::I>().get<int>(0) - 1;
                const std::size_t j = record.getItem<AQUNUM::J>().get<int>(0) - 1;
                const std::size_t k = record.getItem<AQUNUM::K>().get<int>(0) - 1;
                const std::size_t global_index = this->getGlobalIndex(i, j, k);
                this->m_aquifer_cells.insert(global_index);

                if (! record.getItem<AQUNUM::DEPTH>().defaultApplied(0))
                    this->m_aquifer_cell_depths.insert_or_assign(global_index, record.getItem<AQUNUM::DEPTH>().getSIDouble(0));

                // Create map global_index -> (PVTNUM, SATNUM) to allow QC during FieldProps creation
                const int pvtnum = record.getItem<AQUNUM::PVT_TABLE_NUM>().defaultApplied(0) ? 0 : record.getItem<AQUNUM::PVT_TABLE_NUM>().get<int>(0);
                const int satnum = record.getItem<AQUNUM::SAT_TABLE_NUM>().defaultApplied(0) ? 0 : record.getItem<AQUNUM::SAT_TABLE_NUM>().get<int>(0);
                this->m_aquifer_cell_tabnums.insert_or_assign(global_index, std::array<int,2>{pvtnum, satnum});
            }
        }
    }

    const std::optional<MapAxes>& EclipseGrid::getMapAxes() const
    {
        return this->m_mapaxes;
    }


    ZcornMapper::ZcornMapper(std::size_t nx , std::size_t ny, std::size_t nz)
        : dims( {{nx,ny,nz}} ),
          stride( {{2 , 4*nx, 8*nx*ny}} ),
          cell_shift( {{0 , 1 , 2*nx , 2*nx + 1 , 4*nx*ny , 4*nx*ny + 1, 4*nx*ny + 2*nx , 4*nx*ny + 2*nx + 1 }})
    {
    }


    /* lower layer:   upper layer  (higher value of z - i.e. lower down in resrvoir).

         2---3           6---7
         |   |           |   |
         0---1           4---5
    */

    std::size_t ZcornMapper::index(std::size_t i, std::size_t j, std::size_t k, int c) const {
        if ((i >= dims[0]) || (j >= dims[1]) || (k >= dims[2]) || (c < 0) || (c >= 8))
            throw std::invalid_argument("Invalid cell argument");

        return i*stride[0] + j*stride[1] + k*stride[2] + cell_shift[c];
    }

    std::size_t ZcornMapper::size() const {
        return dims[0] * dims[1] * dims[2] * 8;
    }

    std::size_t ZcornMapper::index(std::size_t g, int c) const {
        int k = g / (dims[0] * dims[1]);
        g -= k * dims[0] * dims[1];

        int j = g / dims[0];
        g -= j * dims[0];

        int i = g;

        return index(i,j,k,c);
    }

    bool ZcornMapper::validZCORN( const std::vector<double>& zcorn) const {
        int sign = zcorn[ this->index(0,0,0,0) ] <= zcorn[this->index(0,0, this->dims[2] - 1,4)] ? 1 : -1;
        for (std::size_t j=0; j < this->dims[1]; j++)
            for (std::size_t i=0; i < this->dims[0]; i++)
                for (std::size_t c=0; c < 4; c++)
                    for (std::size_t k=0; k < this->dims[2]; k++) {
                        /* Between cells */
                        if (k > 0) {
                            std::size_t index1 = this->index(i,j,k-1,c+4);
                            std::size_t index2 = this->index(i,j,k,c);
                            if ((zcorn[index2] - zcorn[index1]) * sign < 0)
                                return false;
                        }

                        /* In cell */
                        {
                            std::size_t index1 = this->index(i,j,k,c);
                            std::size_t index2 = this->index(i,j,k,c+4);
                            if ((zcorn[index2] - zcorn[index1]) * sign < 0)
                                return false;
                        }
                    }

        return true;
    }


    std::size_t ZcornMapper::fixupZCORN( std::vector<double>& zcorn) {
        int sign = zcorn[ this->index(0,0,0,0) ] <= zcorn[this->index(0,0, this->dims[2] - 1,4)] ? 1 : -1;
        std::size_t cells_adjusted = 0;

        for (std::size_t k=0; k < this->dims[2]; k++)
            for (std::size_t j=0; j < this->dims[1]; j++)
                for (std::size_t i=0; i < this->dims[0]; i++)
                    for (std::size_t c=0; c < 4; c++) {
                        /* Cell to cell */
                        if (k > 0) {
                            std::size_t index1 = this->index(i,j,k-1,c+4);
                            std::size_t index2 = this->index(i,j,k,c);

                            if ((zcorn[index2] - zcorn[index1]) * sign < 0 ) {
                                zcorn[index2] = zcorn[index1];
                                cells_adjusted++;
                            }
                        }

                        /* Cell internal */
                        {
                            std::size_t index1 = this->index(i,j,k,c);
                            std::size_t index2 = this->index(i,j,k,c+4);

                            if ((zcorn[index2] - zcorn[index1]) * sign < 0 ) {
                                zcorn[index2] = zcorn[index1];
                                cells_adjusted++;
                            }
                        }
                    }
        return cells_adjusted;
    }

    CoordMapper::CoordMapper(std::size_t nx_, std::size_t ny_) :
        nx(nx_),
        ny(ny_)
    {
    }

    std::size_t CoordMapper::size() const {
        return (this->nx + 1) * (this->ny + 1) * 6;
    }


    std::size_t CoordMapper::index(std::size_t i, std::size_t j, std::size_t dim, std::size_t layer) const {
        if (i > this->nx)
            throw std::invalid_argument("Out of range");

        if (j > this->ny)
            throw std::invalid_argument("Out of range");

        if (dim > 2)
            throw std::invalid_argument("Out of range");

        if (layer > 1)
            throw std::invalid_argument("Out of range");

        return 6*( i + j*(this->nx + 1) ) +  layer * 3 + dim;
    }
}



namespace Opm {
    EclipseGridLGR::EclipseGridLGR(const std::string& self_label, const std::string& father_label_,
                                   std::size_t nx, std::size_t ny, std::size_t nz,
                                   const vec_size_t& father_lgr_index, const std::array<int,3>& low_fatherIJK_,
                                   const std::array<int,3>& up_fatherIJK_)
    : EclipseGrid(nx,ny,nz), father_label(father_label_), father_global(father_lgr_index), 
                             low_fatherIJK(low_fatherIJK_), up_fatherIJK(up_fatherIJK_)
    {
        init_father_global();
        lgr_label= self_label;
    }

    std::vector<int> EclipseGridLGR::save_hostnum(void) const
    {
        std::vector<int> hostnum{this->m_hostnum};
        std::transform(hostnum.begin(),hostnum.end(), hostnum.begin(), [](int a){return a+1;});
        return hostnum;
    }

    void EclipseGridLGR::get_label_child_to_top_father(std::vector<std::reference_wrapper<const std::string>>& list) const
    {
        list.push_back(std::ref(lgr_label));
        if (father_label != "GLOBAL")
        {
            getLGRCell(father_label).get_label_child_to_top_father(list);
        }
    }

    void EclipseGridLGR::get_global_index_child_to_top_father(std::vector<std::size_t> & list,
                                                              std::size_t i, std::size_t j, std::size_t k) const
    {        
        get_global_index_child_to_top_father(list, activeIndex(i,j,k));
    }

    void EclipseGridLGR::get_global_index_child_to_top_father(std::vector<std::size_t> & list, std::size_t global_ind) const
    {    
        if (list.empty())
        {
            list.push_back(global_ind);
        }
        std::size_t father_id = get_hostnum(global_ind);
        list.push_back(father_id);        
        if (father_label != "GLOBAL")
        {
            getLGRCell(father_label).get_global_index_child_to_top_father(list, father_id);
        }
    }


    std::optional<std::reference_wrapper<EclipseGridLGR>> 
    EclipseGridLGR::get_child_LGR_cell(const std::string& lgr_tag) const
    {
        std::optional<std::reference_wrapper<EclipseGridLGR>> lgr_found = std::nullopt;
        // Otherwise, search recursively within children.
        if (lgr_tag == lgr_label) 
        {
            lgr_found = std::ref(const_cast<EclipseGridLGR&>(*this));
        }
        else
        {
            for (auto &child : this->lgr_children_cells) {
                if ((child.get_lgr_tag() == lgr_tag) && (lgr_found == std::nullopt)) 
                {
                    lgr_found = std::ref(const_cast<EclipseGridLGR&>(child));
                }
                else
                {
                    return child.get_child_LGR_cell(lgr_tag);
                }
            }
        }
        return lgr_found;
    }

    void EclipseGridLGR::set_hostnum(std::vector<int>& hostnum)
    {
        m_hostnum = hostnum;
    }
    void EclipseGridLGR::set_lgr_refinement(const std::string& lgr_tag, const std::vector<double>& coord, const std::vector<double>& zcorn)
    {
        if (lgr_tag == lgr_label)
        {
            this->set_lgr_refinement(coord, zcorn);
        } else
        {
            for (auto& lgr_cell : lgr_children_cells) {
                lgr_cell.set_lgr_refinement(lgr_tag, coord, zcorn);
            }   
        }
    }
    void EclipseGridLGR::set_lgr_refinement(const std::vector<double>& coord, const std::vector<double>& zcorn)
    {
        m_coord = coord;
        m_zcorn = zcorn;
    }

    void EclipseGridLGR::init_father_global()
    {
        std::sort(father_global.begin(),father_global.end());
    }

    const EclipseGridLGR::vec_size_t& EclipseGridLGR::getFatherGlobalID() const
    {
        return father_global;
    }

    void EclipseGridLGR::save(Opm::EclIO::EclOutput& egridfile, const std::vector<Opm::NNCdata>& nnc, const Opm::UnitSystem& units) const {
        const auto lgr_name_label = std::vector{ Opm::EclIO::PaddedOutputString<8>{ lgr_label }};
        egridfile.write("LGR",lgr_name_label);
        //std::vector<Opm::EclIO::PaddedOutputString<8>> lgr_father_name_label;
        std::vector<std::string> lgr_father_name_label;
        if (lgr_level_father == 0){
            lgr_father_name_label.push_back("");
        }
        else {
            lgr_father_name_label.push_back(father_label);
        }
        egridfile.write("LGRPARNT", lgr_father_name_label);

        constexpr auto length = ::Opm::UnitSystem::measure::length;

        const std::array<int, 3> dims = getNXYZ();

        // Preparing vectors to be saved

        // create coord vector of floats with input units, converted from SI
        std::vector<float> coord_f;
        coord_f.resize(m_coord.size());

        auto convert_length = [&units](const double x) { return static_cast<float>(units.from_si(length, x)); };

        if (m_input_coord.has_value()) {
            std::transform(m_input_coord.value().begin(), m_input_coord.value().end(), coord_f.begin(), convert_length);
        } else {
            std::transform(m_coord.begin(), m_coord.end(), coord_f.begin(), convert_length);
        }

        // create zcorn vector of floats with input units, converted from SI
        std::vector<float> zcorn_f;
        zcorn_f.resize(m_zcorn.size());

        if (m_input_zcorn.has_value()) {
            std::transform(m_input_zcorn.value().begin(), m_input_zcorn.value().end(), zcorn_f.begin(), convert_length);
        } else {
            std::transform(m_zcorn.begin(), m_zcorn.end(), zcorn_f.begin(), convert_length);
        }

        m_input_coord.reset();
        m_input_zcorn.reset();

        // corner point grid

        std::vector<int> gridhead(100,0);
        // GLOBAL and LGR Gridhead
        gridhead[0] = 1;                    // corner point grid
        gridhead[1] = dims[0];              // nI
        gridhead[2] = dims[1];              // nJ
        gridhead[3] = dims[2];              // nK
        gridhead[4] = lgr_level;            // LGR index

        // LGR Exclusive Gridhead Flags
        gridhead[24] = 1;                   // number of reservoirs
        gridhead[25] = 1;                   // number of coordinate line seg
        gridhead[26] = 0;                   // NTHETA =0 non-radial
        gridhead[27] = low_fatherIJK[0] + 1;// Lower I-index-host
        gridhead[28] = low_fatherIJK[1] + 1;// Lower J-index-host
        gridhead[29] = low_fatherIJK[2] + 1;// Lower K-index-host
        gridhead[30] = up_fatherIJK[0] + 1; // Upper I-index-host
        gridhead[31] = up_fatherIJK[1] + 1; // Upper J-index-host
        gridhead[32] = up_fatherIJK[2] + 1; // Upper K-index-host

        [[maybe_unused]] std::vector<int> nnchead(10, 0);
        [[maybe_unused]] std::vector<int> nnc1;
        [[maybe_unused]] std::vector<int> nnc2;

        for (const NNCdata& n : nnc ) {
            nnc1.push_back(n.cell1 + 1);
            nnc2.push_back(n.cell2 + 1);
        }

        nnchead[0] = nnc1.size();

        std::vector<int> endgrid = {};

        // Writing vectors to egrid file

        egridfile.write("GRIDHEAD", gridhead);

        egridfile.write("COORD", coord_f);
        egridfile.write("ZCORN", zcorn_f);

        egridfile.write("ACTNUM", m_actnum);
        egridfile.write("HOSTNUM", save_hostnum());
        egridfile.write("ENDGRID", endgrid);
        egridfile.write("ENDLGR", endgrid);
        for (std::size_t index: m_print_order_lgr_cells ){
            lgr_children_cells[index].save(egridfile, nnc, units);
        }
    }


      void EclipseGridLGR::save_nnc(Opm::EclIO::EclOutput& egridfile) const{

        //std::vector<int> nnchead(10, 0);
        std::vector<int> nnc1;
        std::vector<int> nnc2;

        /* to be implemented based on how the NNC is input */
        // for (const NNCdata& n : nnc ) {
        //     nnc1.push_back(n.cell1 + 1);
        //     nnc2.push_back(n.cell2 + 1);
        // }
        egridfile.write("NNCL", nnc1);
        egridfile.write("NNCG", nnc2);
        //nnchead[0] = nnc1.size();
      }

}
