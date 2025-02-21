/*
  Copyright 2013 Statoil ASA.

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

#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>

#include <opm/io/eclipse/rst/connection.hpp>

#include <opm/common/OpmLog/KeywordLocation.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/utility/ActiveGridCells.hpp>
#include <opm/common/utility/numeric/linearInterpolation.hpp>

#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FieldPropsManager.hpp>
#include <opm/input/eclipse/EclipseState/Grid/GridDims.hpp>

#include <opm/input/eclipse/Schedule/ScheduleGrid.hpp>
#include <opm/input/eclipse/Schedule/Well/Connection.hpp>
#include <opm/input/eclipse/Schedule/Well/WDFAC.hpp>
#include <opm/input/eclipse/Schedule/WellTraj/RigEclipseWellLogExtractor.hpp>

#include <opm/input/eclipse/Units/Units.hpp>

#include <opm/input/eclipse/Deck/DeckRecord.hpp>

#include <external/resinsight/LibCore/cvfVector3.h>
#include <external/resinsight/ReservoirDataModel/RigHexIntersectionTools.h>
#include <external/resinsight/ReservoirDataModel/RigWellLogExtractionTools.h>
#include <external/resinsight/ReservoirDataModel/RigWellLogExtractor.h>
#include <external/resinsight/ReservoirDataModel/RigWellPath.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <fmt/format.h>

#include <stddef.h>

namespace Opm {

    // Compute direction permutation corresponding to completion's
    // direction.  First two elements of return value are directions
    // perpendicular to completion while last element is direction
    // along completion.
    std::array<std::size_t, 3>
    directionIndices(const Opm::Connection::Direction direction)
    {
        switch (direction) {
        case Opm::Connection::Direction::X:
            return {{ 1, 2, 0 }};

        case Opm::Connection::Direction::Y:
            return {{ 2, 0, 1}};

        case Opm::Connection::Direction::Z:
            return {{ 0, 1, 2 }};
        }

        // All enum values should be handled above. Therefore we should
        // never reach this one. Anyway for the sake of reduced warnings we
        // throw an exception.
        throw std::invalid_argument {
            fmt::format("Unhandled direction enumeration value '{}'",
                        Opm::Connection::Direction2String(direction))
        };
    }

    // Permute (diagonal) permeability components according to
    // completion's direction.
    std::array<double, 3>
    permComponents(const Opm::Connection::Direction direction,
                   const std::array<double,3>&      perm)
    {
        const auto p = directionIndices(direction);

        return {{ perm[ p[0] ] ,
                  perm[ p[1] ] ,
                  perm[ p[2] ] }};
    }

    // Permute cell's geometric extent according to completion's
    // direction.  Honour net-to-gross ratio.
    //
    // Note: 'extent' is intentionally accepted by modifiable value
    // rather than reference-to-const to support NTG manipulation.
    std::array<double, 3>
    effectiveExtent(const Opm::Connection::Direction direction,
                    const double                     ntg,
                    std::array<double,3>             extent)
    {
        // Vertical extent affected by net-to-gross ratio.
        extent[2] *= ntg;

        const auto p = directionIndices(direction);

        return {{ extent[ p[0] ] ,
                  extent[ p[1] ] ,
                  extent[ p[2] ] }};
    }

    // Compute Peaceman's effective radius of single completion.
    double effectiveRadius(const std::array<double,3>& K,
                           const std::array<double,3>& D)
    {
        const double K01   = K[0] / K[1];
        const double K10   = K[1] / K[0];

        const double D0_sq = D[0] * D[0];
        const double D1_sq = D[1] * D[1];

        const double num = std::sqrt((std::sqrt(K10) * D0_sq) +
                                     (std::sqrt(K01) * D1_sq));
        const double den = std::pow(K01, 0.25) + std::pow(K10, 0.25);

        // Note: Analytic constant 0.28 derived for infintely sized
        // formation with repeating well placement.
        return 0.28 * (num / den);
    }

    double peacemanDenominator(const double r0,
                               const double rw,
                               const double skin_factor)
    {
        return std::log(r0 / std::min(rw, r0)) + skin_factor;
    }

    double peacemanDenominator(const Opm::Connection::CTFProperties& ctf_props)
    {
        return peacemanDenominator(ctf_props.r0, ctf_props.rw, ctf_props.skin_factor);
    }

    double staticForchheimerCoefficient(const Opm::Connection::CTFProperties& ctf_props,
                                        const double                          porosity,
                                        const Opm::WDFAC&                     wdfac)
    {
        // Reference/background permeability against which to scale cell
        // permeability for model evaluation.
        constexpr auto Kref = 1.0*Opm::prefix::milli*Opm::unit::darcy;

        const auto& corr_coeff = wdfac.getDFactorCorrelationCoefficients();

        return corr_coeff.coeff_a
            * std::pow(ctf_props.Ke / Kref, corr_coeff.exponent_b)
            * std::pow(porosity           , corr_coeff.exponent_c)
            * ctf_props.Ke / (ctf_props.connection_length * ctf_props.rw);
    }

    // Calculate permeability thickness Kh for line segment in a cell for x,y,z directions
    std::array<double, 3>
    permThickness(const external::cvf::Vec3d& effective_connection,
                  const std::array<double,3>& cell_perm,
                  const double ntg)
    {
        auto perm_thickness = std::array {
            effective_connection[0],
            effective_connection[1],
            effective_connection[2] * ntg,
        };

        const auto direction = std::array {
            Opm::Connection::Direction::X,
            Opm::Connection::Direction::Y,
            Opm::Connection::Direction::Z,
        };

        for (size_t i = 0; i < 3; ++i) {
            const auto K = permComponents(direction[i], cell_perm);
            perm_thickness[i] *= std::sqrt(K[0] * K[1]);
        }

        return perm_thickness;
    }

    // Calculate directional (x,y,z) peaceman connection factors CFx, CFy, CFz
    std::array<double, 3>
    connectionFactor(const std::array<double,3>& cell_perm,
                     const std::array<double,3>& cell_size,
                     const double ntg,
                     const std::array<double,3>& Kh,
                     const double rw,
                     const double skin_factor)
    {
        auto connection_factor = std::array<double, 3>{};

        const auto direction = std::array {
            Opm::Connection::Direction::X,
            Opm::Connection::Direction::Y,
            Opm::Connection::Direction::Z,
        };

        const double angle = 6.2831853071795864769252867665590057683943387987502116419498;

        for (size_t i = 0; i < 3; ++i) {
            const auto K = permComponents(direction[i], cell_perm);
            const auto D = effectiveExtent(direction[i], ntg, cell_size);
            const auto r0 = effectiveRadius(K, D);

            connection_factor[i] = angle * Kh[i]
                / peacemanDenominator(r0, rw, skin_factor);
        }

        return connection_factor;
    }

} // anonymous namespace

namespace Opm {

    WellConnections::WellConnections(const Connection::Order order,
                                     const int               headIArg,
                                     const int               headJArg)
        : m_ordering(order)
        , headI     (headIArg)
        , headJ     (headJArg)
    {}

    WellConnections::WellConnections(const Connection::Order        order,
                                     const int                      headIArg,
                                     const int                      headJArg,
                                     const std::vector<Connection>& connections)
        : m_ordering   (order)
        , headI        (headIArg)
        , headJ        (headJArg)
        , m_connections(connections)
    {}

    WellConnections WellConnections::serializationTestObject()
    {
        WellConnections result;
        result.m_ordering = Connection::Order::DEPTH;
        result.headI = 1;
        result.headJ = 2;
        result.m_connections = {Connection::serializationTestObject()};

        return result;
    }

    std::vector<const Connection*>
    WellConnections::output(const EclipseGrid& grid) const
    {
        auto out = std::vector<const Connection*>{};
        out.reserve(this->m_connections.size());

        if (this->m_connections.empty()) {
            return out;
        }

        for (const auto& conn : this->m_connections) {
            if (grid.isCellActive(conn.getI(), conn.getJ(), conn.getK())) {
                out.push_back(&conn);
            }
        }

        if (! this->m_connections[0].attachedToSegment() &&
            (this->m_ordering != Connection::Order::INPUT))
        {
            std::sort(out.begin(), out.end(), [](const Opm::Connection* conn1,
                                                 const Opm::Connection* conn2)
            {
                return conn1->sort_value() < conn2->sort_value();
            });
        }

        return out;
    }

    bool WellConnections::prepareWellPIScaling()
    {
        auto update = false;
        for (auto& conn : this->m_connections) {
            update = conn.prepareWellPIScaling() || update;
        }

        return update;
    }

    void WellConnections::applyWellPIScaling(const double       scaleFactor,
                                             std::vector<bool>& scalingApplicable)
    {
        scalingApplicable.resize(std::max(scalingApplicable.size(), this->m_connections.size()), true);

        auto i = std::size_t{0};
        for (auto& conn : this->m_connections) {
            if (scalingApplicable[i]) {
                scalingApplicable[i] = conn.applyWellPIScaling(scaleFactor);
            }

            ++i;
        }
    }

    void WellConnections::addConnection(const int i, const int j, const int k,
                                        const std::size_t                global_index,
                                        const int                        complnum,
                                        const Connection::State          state,
                                        const double                     depth,
                                        const Connection::CTFProperties& ctf_props,
                                        const int                        satTableId,
                                        const Connection::Direction      direction,
                                        const Connection::CTFKind        ctf_kind,
                                        const std::size_t                seqIndex,
                                        const bool                       defaultSatTabId)
    {
        const int conn_i = (i < 0) ? this->headI : i;
        const int conn_j = (j < 0) ? this->headJ : j;

        this->m_connections.emplace_back(conn_i, conn_j, k, global_index, complnum,
                                         state, direction, ctf_kind, satTableId,
                                         depth, ctf_props, seqIndex, defaultSatTabId);
    }

    void WellConnections::addConnection(const int i, const int j, const int k,
                                        const std::size_t                global_index,
                                        const Connection::State          state,
                                        const double                     depth,
                                        const Connection::CTFProperties& ctf_props,
                                        const int                        satTableId,
                                        const Connection::Direction      direction,
                                        const Connection::CTFKind        ctf_kind,
                                        const std::size_t                seqIndex,
                                        const bool                       defaultSatTabId)
    {
        const auto complnum = static_cast<int>(this->m_connections.size()) + 1;

        this->addConnection(i, j, k, global_index, complnum, state,
                            depth, ctf_props,
                            satTableId, direction, ctf_kind,
                            seqIndex, defaultSatTabId);
    }

    void WellConnections::loadCOMPDAT(const DeckRecord&      record,
                                      const ScheduleGrid&    grid,
                                      const std::string&     wname,
                                      const WDFAC&           wdfac,
                                      const KeywordLocation& location)
    {
        const auto& itemI = record.getItem("I");
        const auto defaulted_I = itemI.defaultApplied(0) || (itemI.get<int>(0) == 0);
        const auto I = !defaulted_I ? itemI.get<int>(0) - 1 : this->headI;

        const auto& itemJ = record.getItem("J");
        const auto defaulted_J = itemJ.defaultApplied(0) || (itemJ.get<int>(0) == 0);
        const auto J = !defaulted_J ? itemJ.get<int>(0) - 1 : this->headJ;

        const auto K1 = record.getItem("K1").get<int>(0) - 1;
        const auto K2 = record.getItem("K2").get<int>(0) - 1;
        const auto state = Connection::StateFromString(record.getItem("STATE").getTrimmedString(0));

        const auto& r0Item = record.getItem("PR");
        const auto& CFItem = record.getItem("CONNECTION_TRANSMISSIBILITY_FACTOR");
        const auto& diameterItem = record.getItem("DIAMETER");
        const auto& KhItem = record.getItem("Kh");
        const auto& satTableIdItem = record.getItem("SAT_TABLE");
        const auto direction = Connection::DirectionFromString(record.getItem("DIR").getTrimmedString(0));

        const auto skin_factor = record.getItem("SKIN").getSIDouble(0);
        const auto d_factor = record.getItem("D_FACTOR").getSIDouble(0);

        int satTableId = -1;
        bool defaultSatTable = true;
        if (satTableIdItem.hasValue(0) && (satTableIdItem.get<int>(0) > 0)) {
            satTableId = satTableIdItem.get<int>(0);
            defaultSatTable = false;
        }

        double rw{};
        if (diameterItem.hasValue(0)) {
            rw = diameterItem.getSIDouble(0) / 2;
        }
        else {
            // The Eclipse100 manual does not specify a default value for the wellbore
            // diameter, but the Opm codebase has traditionally implemented a default
            // value of one foot. The same default value is used by Eclipse300.
            rw = 0.5*unit::feet;
        }

        // Angle of completion exposed to flow.  We assume centre
        // placement so there's complete exposure (= 2\pi).
        const auto angle = 6.2831853071795864769252867665590057683943387987502116419498;

        for (int k = K1; k <= K2; ++k) {
            const auto& cell = grid.get_cell(I, J, k);
            if (!cell.is_active()) {
                auto msg = fmt::format(R"(Problem with COMPDAT keyword
In {} line {}
The cell ({},{},{}) in well {} is not active and the connection will be ignored)",
                                       location.filename, location.lineno,
                                       I + 1, J + 1, k + 1, wname);

                OpmLog::warning(msg);
                continue;
            }

            const auto& props = cell.props;

            auto ctf_props = ::Opm::Connection::CTFProperties{};
            ctf_props.rw = rw;
            ctf_props.skin_factor = skin_factor;
            ctf_props.d_factor = d_factor;

            if (defaultSatTable) {
                satTableId = props->satnum;
            }

            ctf_props.r0 = -1.0;
            if (r0Item.hasValue(0)) {
                ctf_props.r0 = r0Item.getSIDouble(0);
            }

            ctf_props.Kh = -1.0;
            if (KhItem.hasValue(0) && (KhItem.getSIDouble(0) > 0.0)) {
                ctf_props.Kh = KhItem.getSIDouble(0);
            }

            ctf_props.CF = -1.0;
            if (CFItem.hasValue(0) && (CFItem.getSIDouble(0) > 0.0)) {
                ctf_props.CF = CFItem.getSIDouble(0);
            }

            const auto D = effectiveExtent(direction, props->ntg, cell.dimensions);
            const auto K = permComponents(direction, { props->permx, props->permy, props->permz });
            ctf_props.Ke = std::sqrt(K[0] * K[1]);

            const auto ctf_kind = (ctf_props.CF < 0.0)
                ? ::Opm::Connection::CTFKind::Defaulted
                : ::Opm::Connection::CTFKind::DeckValue;

            // We start with the absolute happy path; both CF and Kh are
            // explicitly given in the deck.
            if ((ctf_props.CF > 0.0) && (ctf_props.Kh > 0.0)) {
                ctf_props.peaceman_denom = angle * ctf_props.Kh / ctf_props.CF;
                goto CF_done;
            }

            // We must calculate CF and Kh from the items in the COMPDAT
            // record and cell properties.
            if (ctf_props.r0 < 0.0) {
                ctf_props.r0 = effectiveRadius(K, D);
            }

            if (const auto peaceman_denom = peacemanDenominator(ctf_props);
                ctf_props.Kh > 0.0)
            {
                // CF < 0
                ctf_props.CF = angle * ctf_props.Kh / peaceman_denom;
                ctf_props.peaceman_denom = peaceman_denom;
            }
            else if (ctf_props.CF > 0.0) { // Kh < 0
                if (KhItem.defaultApplied(0) || (KhItem.get<double>(0) < 0.0)) {
                    // Kh explicitly defaulted.  Derive compatible Kh value
                    // from specified CTF, r0, rw, and skin factor.
                    ctf_props.Kh = ctf_props.CF * peaceman_denom / angle;
                }
                else {
                    // Kh = 0 entered in item 10 of COMPDAT.  Compute Kh
                    // from permeability and length of perforation interval
                    // and request a compatible pressure equivalent radius
                    // (r0) be calculated.
                    ctf_props.Kh = ctf_props.Ke * D[2];

                    // Defer r0 calculation to CF_done: r0 < 0 case.
                    ctf_props.r0 = -1.0;
                }

                ctf_props.peaceman_denom = angle * ctf_props.Kh / ctf_props.CF;
            }
            else {              // (CF < 0) && (Kh < 0)
                ctf_props.Kh = ctf_props.Ke * D[2];
                ctf_props.CF = angle * ctf_props.Kh / peaceman_denom;
                ctf_props.peaceman_denom = peaceman_denom;
            }

        CF_done:
            if (ctf_props.r0 < 0.0) {
                ctf_props.r0 = RestartIO::RstConnection::
                    inverse_peaceman(ctf_props.CF, ctf_props.Kh,
                                     ctf_props.rw,
                                     ctf_props.skin_factor);
            }

            // Length of the well perforation interval.
            ctf_props.connection_length = ctf_props.Kh / ctf_props.Ke;

            // Area equivalent radius of the grid block.  Used by the
            // PolymerMW module.
            ctf_props.re = std::sqrt(D[0] * D[1] / angle * 2);

            ctf_props.static_dfac_corr_coeff =
                staticForchheimerCoefficient(ctf_props, props->poro, wdfac);

            auto prev = std::find_if(this->m_connections.begin(),
                                     this->m_connections.end(),
                                     [I, J, k](const Connection& c)
                                     {
                                         return c.sameCoordinate(I, J, k);
                                     });

            if (prev == this->m_connections.end()) {
                const std::size_t noConn = this->m_connections.size();
                this->addConnection(I, J, k, cell.global_index, state,
                                    cell.depth, ctf_props, satTableId,
                                    direction, ctf_kind,
                                    noConn, defaultSatTable);
            }
            else {
                const auto compl_num = prev->complnum();
                const auto css_ind = prev->sort_value();
                const auto conSegNo = prev->segment();
                const auto perf_range = prev->perf_range();

                *prev = Connection {
                    I, J, k, cell.global_index, compl_num,
                    state, direction, ctf_kind, satTableId,
                    cell.depth, ctf_props,
                    css_ind, defaultSatTable
                };

                prev->updateSegment(conSegNo, cell.depth, css_ind, perf_range);
            }
        }
    }

    void WellConnections::loadCOMPTRAJ(const DeckRecord&      record,
                                       const ScheduleGrid&    grid,
                                       const std::string&     wname,
                                       const KeywordLocation& location,
                                       external::cvf::ref<external::cvf::BoundingBoxTree>& cellSearchTree)
    {
        const auto& perf_top = record.getItem("PERF_TOP");
        const auto& perf_bot = record.getItem("PERF_BOT");

        const auto& CFItem = record.getItem("CONNECTION_TRANSMISSIBILITY_FACTOR");
        const auto& diameterItem = record.getItem("DIAMETER");
        const auto& KhItem = record.getItem("Kh");
        const auto skin_factor = record.getItem("SKIN").getSIDouble(0);
        const auto d_factor = record.getItem("D_FACTOR").getSIDouble(0);
        const auto& satTableIdItem = record.getItem("SAT_TABLE");
        const auto state = Connection::StateFromString(record.getItem("STATE").getTrimmedString(0));
	//
	m_top.push_back(perf_top.getSIDouble(0));
	m_bot.push_back(perf_bot.getSIDouble(0));
	m_skinn.push_back(skin_factor);
	m_d_factor.push_back(d_factor);

	int copml = -1;
        m_copml.push_back(copml);
	m_state.push_back(state);


	
        int satTableId = -1;
        bool defaultSatTable = true;
        if (satTableIdItem.hasValue(0) && (satTableIdItem.get<int>(0) > 0)) {
            satTableId = satTableIdItem.get<int>(0);
            defaultSatTable = false;
        }
	m_tab.push_back(satTableId);
        double rw{};
        if (diameterItem.hasValue(0)) {
            rw = diameterItem.getSIDouble(0) / 2;
        }
	
        else {
            // The Eclipse100 manual does not specify a default value for the wellbore
            // diameter, but the Opm codebase has traditionally implemented a default
            // value of one foot. The same default value is used by Eclipse300.
	  rw = 0.5*unit::feet; // NB probably should be handled in input
        }
	m_rw.push_back(rw);
        // Get the grid
        const auto& ecl_grid = grid.get_grid();

        std::vector<external::cvf::Vec3d> points;
        std::vector<double> measured_depths;

        // Calulate the x,y,z coordinates of the begin and end of a perforation
        external::cvf::Vec3d p_top, p_bot;
        double top = perf_top.getSIDouble(0), bot = perf_bot.getSIDouble(0);
        for (size_t i = 0; i < 3 ; ++i) {
            p_top[i] = linearInterpolation(this->md, this->coord[i], top);
            p_bot[i] = linearInterpolation(this->md, this->coord[i], bot);
        }
        points.push_back(p_top);
        measured_depths.push_back(top);

        external::cvf::ref<external::RigWellPath> wellPathGeometry { new external::RigWellPath };
        points.reserve(this->coord[0].size());
        measured_depths.reserve(this->coord[0].size());
        for (size_t i = 0; i < coord[0].size(); ++i) {
            if (this->md[i] > top and this->md[i] < bot) {
                points.push_back(external::cvf::Vec3d(coord[0][i], coord[1][i], coord[2][i]));
                measured_depths.push_back(this->md[i]);
            }
        }

        points.push_back(p_bot);
        measured_depths.push_back(bot);

        wellPathGeometry->setWellPathPoints(points);
        wellPathGeometry->setMeasuredDepths(measured_depths);

        external::cvf::ref<external::RigEclipseWellLogExtractor> e {
            new external::RigEclipseWellLogExtractor {
                wellPathGeometry.p(), *ecl_grid, cellSearchTree
            }
        };

        // Keep the AABB search tree of the grid to avoid redoing an
        // expensive calulation.
        cellSearchTree = e->getCellSearchTree();

        // This gives the intersected grid cells IJK, cell face entrance &
        // exit cell face point and connection length.
        auto intersections = e->cellIntersectionInfosAlongWellPath();

	
	
        for (size_t is = 0; is < intersections.size(); ++is) {
            const auto ijk = ecl_grid->getIJK(intersections[is].globCellIndex);

            // When using WELTRAJ & COMPTRAJ one may use default settings in
            // WELSPECS for headI/J and let the headI/J be calculated by the
            // trajectory data.
            //
            // If these defaults are used the headI/J are set to the first
            // intersection.
            if (is == 0) {
                if (this->headI < 0) { this->headI = ijk[0]; }
                if (this->headJ < 0) { this->headJ = ijk[1]; }
            }

            const auto& cell = grid.get_cell(ijk[0], ijk[1], ijk[2]);

            if (!cell.is_active()) {
                const auto msg = fmt::format(R"(Problem with COMPTRAJ keyword
In {} line {}
The cell ({},{},{}) in well {} is not active and the connection will be ignored)",
                                             location.filename,
                                             location.lineno,
                                             ijk[0] + 1,
                                             ijk[1] + 1,
                                             ijk[2] + 1, wname);
                OpmLog::warning(msg);

                continue;
            }

            const auto& props = cell.props;

            auto ctf_props = Connection::CTFProperties{};
            ctf_props.rw = rw;
            ctf_props.skin_factor = skin_factor;
            ctf_props.d_factor = d_factor;

            if (defaultSatTable) {
                satTableId = props->satnum;
            }

            ctf_props.r0 = -1.0;//NB never set later
            ctf_props.Kh = -1.0;
            if (KhItem.hasValue(0) && (KhItem.getSIDouble(0) > 0.0)) {
                ctf_props.Kh = KhItem.getSIDouble(0);
            }
	    if(is==0){
	      m_Kh.push_back(ctf_props.Kh);
	    }
            ctf_props.CF = -1.0;
            if (CFItem.hasValue(0) && (CFItem.getSIDouble(0) > 0.0)) {
                ctf_props.CF = CFItem.getSIDouble(0);
            }
	    if(is == 0){
	      m_conn.push_back(ctf_props.CF);
	    }
            const auto cell_perm = std::array {
                props->permx, props->permy, props->permz
            };

            auto ctf_kind = ::Opm::Connection::CTFKind::DeckValue;
            if ((ctf_props.CF < 0.0) && (ctf_props.Kh < 0.0)) {
                // We must calculate CF and Kh from the items in the
                // COMPTRAJ record and cell properties.
                ctf_kind = ::Opm::Connection::CTFKind::Defaulted;

                const auto& connection_vector =
                    intersections[is].intersectionLengthsInCellCS;

                const auto perm_thickness =
                    permThickness(connection_vector, cell_perm, props->ntg);

                const auto connection_factor =
                    connectionFactor(cell_perm, cell.dimensions, props->ntg,
                                     perm_thickness, rw, skin_factor);

                ctf_props.connection_length = connection_vector.length();

                ctf_props.CF = std::hypot(connection_factor[0],
                                          connection_factor[1],
                                          connection_factor[2]);

                ctf_props.Kh = std::hypot(perm_thickness[0],
                                          perm_thickness[1],
                                          perm_thickness[2]);
            }
            else if (! ((ctf_props.CF > 0.0) && (ctf_props.Kh > 0.0))) {
                auto msg = fmt::format(R"(Problem with COMPTRAJ keyword
In {} line {}
CF and Kh items for well {} must both be specified or both defaulted/negative)",
                                       location.filename, location.lineno, wname);

                throw std::logic_error(msg);
            }

            // Todo: check what needs to be done for polymerMW module, see
            // loadCOMPDAT used by the PolymerMW module

            const auto direction = ::Opm::Connection::Direction::Z;

            ctf_props.re = -1;

            {
                const auto K = permComponents(direction, cell_perm);
                ctf_props.Ke = std::sqrt(K[0] * K[1]);
            }

            auto prev = std::find_if(this->m_connections.begin(),
                                     this->m_connections.end(),
                                     [&ijk](const Connection& c)
                                     { return c.sameCoordinate(ijk[0], ijk[1], ijk[2]); });

            if (prev == this->m_connections.end()) {
                const std::size_t noConn = this->m_connections.size();
                this->addConnection(ijk[0], ijk[1], ijk[2],
                                    cell.global_index, state,
                                    cell.depth, ctf_props, satTableId,
                                    direction, ctf_kind,
                                    noConn, defaultSatTable);
		// this->addConnection(-ijk[0], -ijk[1], -ijk[2],
                //                     -cell.global_index, state,
                //                     cell.depth, ctf_props, satTableId,
                //                     direction, ctf_kind,
                //                     noConn, defaultSatTable);
            }
            else {
                const auto compl_num = prev->complnum();
                const auto css_ind = prev->sort_value();
                const auto conSegNo = prev->segment();
                const auto perf_range = prev->perf_range();

                *prev = Connection {
                    ijk[0], ijk[1], ijk[2],
                    cell.global_index, compl_num,
                    state, direction, ctf_kind, satTableId,
                    cell.depth, ctf_props,
                    css_ind, defaultSatTable
                };

                prev->updateSegment(conSegNo, cell.depth, css_ind, *perf_range);
            }
        }
    }


  
    void WellConnections::loadWELTRAJ(const DeckRecord& record,
                                      [[maybe_unused]] const ScheduleGrid&    grid,
                                      [[maybe_unused]] const std::string&     wname,
                                      [[maybe_unused]] const KeywordLocation& location)
    {
        double x = record.getItem("X").getSIDouble(0);
        double y = record.getItem("Y").getSIDouble(0);
        const auto& mapaxes = grid.get_grid()->getMapAxes();
        if (mapaxes.has_value())
            mapaxes.value().inv_transform(x, y);

        this->coord[0].push_back(x);
        this->coord[1].push_back(y);
        this->coord[2].push_back(record.getItem("TVD").getSIDouble(0));

        this->md.push_back(record.getItem("MD").getSIDouble(0));
    }

    void WellConnections::applyDFactorCorrelation(const ScheduleGrid& grid,
                                                  const WDFAC&        wdfac)
    {
        for (auto& conn : this->m_connections) {
            const auto& complCell = grid.get_cell(conn.getI(), conn.getJ(), conn.getK());
            if (! complCell.is_active()) {
                continue;
            }

            conn.setStaticDFacCorrCoeff
                (staticForchheimerCoefficient(conn.ctfProperties(),
                                              complCell.props->poro, wdfac));
        }
    }

    std::size_t WellConnections::size() const
    {
        return m_connections.size();
    }

    std::size_t WellConnections::num_open() const
    {
        return std::count_if(this->m_connections.begin(),
                             this->m_connections.end(),
                             [] (const Connection& c) { return c.state() == Connection::State::OPEN; });
    }

    bool WellConnections::empty() const
    {
        return this->size() == size_t{0};
    }

    const Connection& WellConnections::get(size_t index) const
    {
        return (*this)[index];
    }

    const Connection& WellConnections::operator[](size_t index) const
    {
        return this->m_connections.at(index);
    }

    const Connection& WellConnections::lowest() const
    {
        if (this->m_connections.empty()) {
            throw std::logic_error("Tried to get lowest connection from empty set");
        }

        const auto max_iter = std::max_element(this->m_connections.begin(),
                                               this->m_connections.end(),
                                               [](const Connection& c1, const Connection& c2)
                                               {
                                                   return c1.depth() < c2.depth();
                                               });

        return *max_iter;
    }

    bool WellConnections::hasGlobalIndex(std::size_t global_index) const
    {
        return std::any_of(this->begin(), this->end(),
                           [global_index](const Connection& conn)
                           {
                               return conn.global_index() == global_index;
                           });
    }

    const Connection&
    WellConnections::getFromIJK(const int i, const int j, const int k) const
    {
        for (size_t ic = 0; ic < size(); ++ic) {
            if (get(ic).sameCoordinate(i, j, k)) {
                return get(ic);
            }
        }

        throw std::runtime_error(" the connection is not found! \n ");
    }

    const Connection& WellConnections::getFromGlobalIndex(std::size_t global_index) const
    {
        auto conn_iter =
            std::find_if(this->begin(), this->end(),
                         [global_index] (const Connection& conn)
                         { return conn.global_index() == global_index; });

        if (conn_iter == this->end()) {
            throw std::logic_error(fmt::format("No connection with global index {}", global_index));
        }

        return *conn_iter;
    }

    Connection& WellConnections::getFromIJK(const int i, const int j, const int k)
    {
        for (size_t ic = 0; ic < size(); ++ic) {
            if (get(ic).sameCoordinate(i, j, k)) {
                return this->m_connections[ic];
            }
        }

        throw std::runtime_error(" the connection is not found! \n ");
    }

    Connection* WellConnections::maybeGetFromGlobalIndex(const std::size_t global_index)
    {
        auto conn_iter =
            std::find_if(this->begin(), this->end(),
                         [global_index] (const Connection& conn)
                         { return conn.global_index() == global_index; });

        if (conn_iter == this->end()) {
            return nullptr;
        }

        return &*conn_iter;
    }

    bool WellConnections::allConnectionsShut() const
    {
        if (this->empty()) {
            return false;
        }

        return std::all_of(this->m_connections.begin(),
                           this->m_connections.end(),
                           [](const Connection& c)
                           { return c.state() == Connection::State::SHUT; });
    }

    void WellConnections::order()
    {
        if (m_connections.empty()) {
            return;
        }

        if (this->m_connections[0].attachedToSegment()) {
            this->orderMSW();
        }
        else if (this->m_ordering == Connection::Order::TRACK) {
            this->orderTRACK();
        }
        else if (this->m_ordering == Connection::Order::DEPTH) {
            this->orderDEPTH();
        }
    }

    void WellConnections::orderMSW()
    {
        std::sort(this->m_connections.begin(), this->m_connections.end(),
                  [](const Opm::Connection& conn1, const Opm::Connection& conn2)
                  {
                      return conn1.sort_value() < conn2.sort_value();
                  });
    }

    void WellConnections::orderTRACK()
    {
        // Find the first connection and swap it into the 0-position.
        const double surface_z = 0.0;
        size_t first_index = findClosestConnection(this->headI, this->headJ, surface_z, 0);
        std::swap(m_connections[first_index], m_connections[0]);

        // Repeat for remaining connections.
        //
        // Note that since findClosestConnection() is O(n), this is an
        // O(n^2) algorithm.  However, it should be acceptable since the
        // expected number of connections is fairly low (< 100).

        if (this->m_connections.empty()) {
            return;
        }

        for (size_t pos = 1; pos < m_connections.size() - 1; ++pos) {
            const auto& prev = m_connections[pos - 1];
            const double prevz = prev.depth();
            size_t next_index = findClosestConnection(prev.getI(), prev.getJ(), prevz, pos);
            std::swap(m_connections[next_index], m_connections[pos]);
        }
    }

    size_t WellConnections::findClosestConnection(int oi, int oj, double oz, size_t start_pos)
    {
        size_t closest = std::numeric_limits<size_t>::max();
        int min_ijdist2 = std::numeric_limits<int>::max();
        double min_zdiff = std::numeric_limits<double>::max();
        for (size_t pos = start_pos; pos < m_connections.size(); ++pos) {
            const auto& connection = m_connections[ pos ];

            const double depth = connection.depth();
            const int ci = connection.getI();
            const int cj = connection.getJ();
            // Using square of distance to avoid non-integer arithmetics.
            const int ijdist2 = (ci - oi) * (ci - oi) + (cj - oj) * (cj - oj);
            if (ijdist2 < min_ijdist2) {
                min_ijdist2 = ijdist2;
                min_zdiff = std::abs(depth - oz);
                closest = pos;
            } else if (ijdist2 == min_ijdist2) {
                const double zdiff = std::abs(depth - oz);
                if (zdiff < min_zdiff) {
                    min_zdiff = zdiff;
                    closest = pos;
                }
            }
        }
        assert(closest != std::numeric_limits<size_t>::max());
        return closest;
    }

    void WellConnections::orderDEPTH()
    {
        std::sort(this->m_connections.begin(), this->m_connections.end(),
                  [](const Opm::Connection& conn1, const Opm::Connection& conn2)
                  {
                      return conn1.depth() < conn2.depth();
                  });
    }

    bool WellConnections::operator==(const WellConnections& rhs) const
    {
        return (this->size() == rhs.size())
            && (this->m_ordering == rhs.m_ordering)
            && (this->coord == rhs.coord)
            && (this->md == rhs.md)
            && std::equal( this->begin(), this->end(), rhs.begin());
    }

    bool WellConnections::operator!=(const WellConnections& rhs) const
    {
        return ! (*this == rhs);
    }

    void WellConnections::filter(const ActiveGridCells& grid)
    {
        auto isInactive = [&grid](const Connection& c) {
            return !grid.cellActive(c.getI(), c.getJ(), c.getK());
        };

        auto new_end = std::remove_if(m_connections.begin(), m_connections.end(), isInactive);
        m_connections.erase(new_end, m_connections.end());
    }

    double WellConnections::segment_perf_length(int segment) const
    {
        double perf_length = 0;
        for (const auto& conn : this->m_connections) {
            if (conn.segment() == segment) {
                const auto& [start_seg, end_seg] = *conn.perf_range();
                perf_length += end_seg - start_seg;
            }
        }

        return perf_length;
    }

    int WellConnections::getHeadI() const
    {
        return this->headI;
    }

    int WellConnections::getHeadJ() const
    {
        return this->headJ;
    }

    const std::vector<double>& WellConnections::getMD() const
    {
        return this->md;
    }

    std::optional<int>
    getCompletionNumberFromGlobalConnectionIndex(const WellConnections& connections,
                                                 const std::size_t      global_index)
    {
        auto connPos = std::find_if(connections.begin(), connections.end(),
            [global_index](const Connection& conn)
        {
            return conn.global_index() == global_index;
        });

        if (connPos == connections.end()) {
            // No connection exists with the requisite 'global_index'
            return {};
        }

        return { connPos->complnum() };
    }
}
