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

#ifndef CONNECTIONSET_HPP_
#define CONNECTIONSET_HPP_

#include <opm/input/eclipse/Schedule/Well/Connection.hpp>
#include <external/resinsight/LibGeometry/cvfBoundingBoxTree.h>
#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <stddef.h>

//only for the template function
#include <opm/input/eclipse/Schedule/WellTraj/RigEclipseWellLogExtractorGrid.hpp>
#include <opm/common/utility/numeric/linearInterpolation.hpp>
#include <external/resinsight/LibCore/cvfVector3.h>
#include <external/resinsight/ReservoirDataModel/RigHexIntersectionTools.h>
#include <external/resinsight/ReservoirDataModel/RigWellLogExtractionTools.h>
#include <external/resinsight/ReservoirDataModel/RigWellLogExtractor.h>
#include <external/resinsight/ReservoirDataModel/RigWellPath.h>

#include <fmt/format.h>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <cassert>
namespace Opm {
    class ActiveGridCells;
    class DeckRecord;
    class EclipseGrid;
    class FieldPropsManager;
    class KeywordLocation;
    class ScheduleGrid;
    class WDFAC;
} // namespace Opm
namespace Opm {
  std::array<double, 3>
    permComponents(const Opm::Connection::Direction direction,
                   const std::array<double,3>&      perm);
  std::array<double, 3>
    permThickness(const external::cvf::Vec3d& effective_connection,
                  const std::array<double,3>& cell_perm,
                  const double ntg);
  std::array<double, 3>
    connectionFactor(const std::array<double,3>& cell_perm,
                     const std::array<double,3>& cell_size,
                     const double ntg,
                     const std::array<double,3>& Kh,
                     const double rw,
                     const double skin_factor);
}



namespace Opm {

    class WellConnections
    {
    public:
        using const_iterator = std::vector<Connection>::const_iterator;

        WellConnections() = default;
        WellConnections(const Connection::Order ordering, const int headI, const int headJ);
        WellConnections(const Connection::Order ordering, const int headI, const int headJ,
                        const std::vector<Connection>& connections);

        static WellConnections serializationTestObject();

        // cppcheck-suppress noExplicitConstructor
        template <class Grid>
        WellConnections(const WellConnections& src, const Grid& grid)
            : m_ordering(src.ordering())
            , headI     (src.headI)
            , headJ     (src.headJ)
        {
            for (const auto& c : src) {
                if (grid.isCellActive(c.getI(), c.getJ(), c.getK())) {
                    this->add(c);
                }
            }
        }

        void add(const Connection& conn)
        {
            this->m_connections.push_back(conn);
        }

        void addConnection(const int i, const int j, const int k,
                           const std::size_t global_index,
                           const Connection::State state,
                           const double depth,
                           const Connection::CTFProperties& ctf_props,
                           const int satTableId,
                           const Connection::Direction direction = Connection::Direction::Z,
                           const Connection::CTFKind ctf_kind = Connection::CTFKind::DeckValue,
                           const std::size_t seqIndex = 0,
                           const bool defaultSatTabId = true);

        void loadCOMPDAT(const DeckRecord&      record,
                         const ScheduleGrid&    grid,
                         const std::string&     wname,
                         const WDFAC&           wdfac,
                         const KeywordLocation& location);

        void loadCOMPTRAJ(const DeckRecord&      record,
                          const ScheduleGrid&    grid,
                          const std::string&     wname,
                          const KeywordLocation& location,
                          external::cvf::ref<external::cvf::BoundingBoxTree>& cellSearchTree);

        void loadWELTRAJ(const DeckRecord&      record,
                         const ScheduleGrid&    grid,
                         const std::string&     wname,
                         const KeywordLocation& location);

        void applyDFactorCorrelation(const ScheduleGrid& grid,
                                     const WDFAC&        wdfac);

        int getHeadI() const;
        int getHeadJ() const;
        const std::vector<double>& getMD() const;
        std::size_t size() const;
        bool empty() const;
        std::size_t num_open() const;
        const Connection& operator[](size_t index) const;
        const Connection& get(size_t index) const;
        const Connection& getFromIJK(const int i, const int j, const int k) const;
        const Connection& getFromGlobalIndex(std::size_t global_index) const;
        const Connection& lowest() const;
        Connection& getFromIJK(const int i, const int j, const int k);
        Connection* maybeGetFromGlobalIndex(const std::size_t global_index);
        bool hasGlobalIndex(std::size_t global_index) const;
        double segment_perf_length(int segment) const;

        const_iterator begin() const { return this->m_connections.begin(); }
        const_iterator end() const { return this->m_connections.end(); }
        auto begin() { return this->m_connections.begin(); }
        auto end() { return this->m_connections.end(); }
        void filter(const ActiveGridCells& grid);
        bool allConnectionsShut() const;
        /// Order connections irrespective of input order.
        /// The algorithm used is the following:
        ///     1. The connection nearest to the given (well_i, well_j)
        ///        coordinates in terms of the connection's (i, j) is chosen
        ///        to be the first connection. If non-unique, choose one with
        ///        lowest z-depth (shallowest).
        ///     2. Choose next connection to be nearest to current in (i, j) sense.
        ///        If non-unique choose closest in z-depth (not logical cartesian k).
        ///
        /// \param[in] well_i  logical cartesian i-coordinate of well head
        /// \param[in] well_j  logical cartesian j-coordinate of well head
        /// \param[in] grid    EclipseGrid object, used for cell depths
        void order();

        bool operator==( const WellConnections& ) const;
        bool operator!=( const WellConnections& ) const;

        Connection::Order ordering() const { return this->m_ordering; }
        std::vector<const Connection *> output(const EclipseGrid& grid) const;

        /// Activate or reactivate WELPI scaling for this connection set.
        ///
        /// Following this call, any WELPI-based scaling will apply to all
        /// connections whose properties are not reset in COMPDAT.
        ///
        /// Returns whether or not this call to prepareWellPIScaling() is
        /// a state change (e.g., no WELPI to active WELPI or WELPI for
        /// some connections to WELPI for all connections).
        bool prepareWellPIScaling();

        /// Scale pertinent connections' CF value by supplied value.  Scaling
        /// factor typically derived from 'WELPI' input keyword and a dynamic
        /// productivity index calculation.  Applicability array specifies
        /// whether or not a particular connection is exempt from scaling.
        /// Empty array means "apply scaling to all eligible connections".
        /// This array is updated on return (entries set to 'false' if
        /// corresponding connection is not eligible).
        void applyWellPIScaling(const double       scaleFactor,
                                std::vector<bool>& scalingApplicable);

        template <class Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(this->m_ordering);
            serializer(this->headI);
            serializer(this->headJ);
            serializer(this->m_connections);
            serializer(this->coord);
            serializer(this->md);
        }

      
      template<class Grid>
      void recomputeConnections(const Grid& grid,
				  external::cvf::ref<external::cvf::BoundingBoxTree>& cellSearchTree)
      {
        
        
        // recalculate well trac with new grid
      for(int seg=0; seg<m_top.size(); ++seg){
	int satTableId = m_tab[seg];
	bool defaultSatTable = false;
	if(satTableId == -1){
	  defaultSatTable = true;
	}
	double rw = m_rw[seg];

	double skin_factor = m_skinn[seg];
	double d_factor = m_d_factor[seg];
	const auto state = m_state[seg];
	
	
	std::vector<external::cvf::Vec3d> points;
        std::vector<double> measured_depths;

        // Calulate the x,y,z coordinates of the begin and end of a perforation
        external::cvf::Vec3d p_top, p_bot;
        double top = m_top[seg], bot = m_bot[seg];
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

        external::cvf::ref<external::RigEclipseWellLogExtractorGrid> e {
            new external::RigEclipseWellLogExtractorGrid {
                wellPathGeometry.p(), grid, cellSearchTree
            }
        };

        // Keep the AABB search tree of the grid to avoid redoing an
        // expensive calulation.
        cellSearchTree = e->getCellSearchTree(); 

        // This gives the intersected grid cells IJK, cell face entrance &
        // exit cell face point and connection length.
        auto intersections = e->cellIntersectionInfosAlongWellPath();

        for (size_t is = 0; is < intersections.size(); ++is) {
	  const auto& cellIdx = intersections[is].globCellIndex; //NB cell is uncompressed
	  
	    std::array<int,3> ijk = {cellIdx, -1, -1}; // just fak it
            
            //const auto& cell = grid.get_cell(cellIdx);//ned to be implemented
	    Opm::CompletedCells::Cell cell;
	    {
	    cell.depth = 1000;//this->grid->getCellDepth(i, j, k);
	    cell.dimensions = {10,10, 1};//this->grid->getCellDimensions(i, j, k);
	    cell.global_index = cellIdx;
	    auto& props = cell.props.emplace(CompletedCells::Cell::Props{});
	    props.active_index = cellIdx; //??
            props.permx = 1e-13;//try_get_value(*this->fp, "PERMX", props.active_index);
            props.permy = 1e-13;//try_get_value(*this->fp, "PERMY", props.active_index);
            props.permz = 1e-13;//try_get_value(*this->fp, "PERMZ", props.active_index);
            props.poro = 0.1;//try_get_value(*this->fp, "PORO", props.active_index);
            props.satnum = 1;//this->fp->get_int("SATNUM").at(props.active_index);
            props.pvtnum = 1;//this->fp->get_int("PVTNUM").at(props.active_index);
            props.ntg = 1;//try_get_ntg_value(*this->fp, "NTG", props.active_index);
	    }
            if (cellIdx < 0) {
                const auto msg = fmt::format(R"(Problem with recalculating COMPTRAJ keyword
In {} line {}
The cell ({},{},{}) in well {} is not active/ and the connection will be ignored)",
                                             cellIdx);
                OpmLog::warning(msg);
		assert(false);
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

            ctf_props.r0 = -1.0; //NB not set later
            ctf_props.Kh = m_Kh[seg];
            

            ctf_props.CF = m_conn[seg];
            
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
CF and Kh items for well {} must both be specified or both defaulted/negative)");

                throw std::logic_error(msg);
            }

            // Todo: check what needs to be done for polymerMW module, see
            // loadCOMPDAT used by the PolymerMW module

            const auto direction = ::Opm::Connection::Direction::Z;

            ctf_props.re = -1;//NB is it used

            {
                const auto K = permComponents(direction, cell_perm);
                ctf_props.Ke = std::sqrt(K[0] * K[1]);
            }

	    //NB when can this happen branched well?
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
      }

      bool hasTraj(){
	size_t np = m_top.size();
	
	assert( np == m_tab.size());
	assert( np == m_Kh.size());
	assert( np == m_conn.size());
	return np > 0;
      }
    private:
        Connection::Order m_ordering { Connection::Order::TRACK };
        int headI{0};
        int headJ{0};
        std::vector<Connection> m_connections{};

      // WELLTRAJ
        std::array<std::vector<double>, 3> coord{};
        std::vector<double> md{};
      // COMPTRAJ
      std::vector<double> m_top;
      std::vector<double> m_bot;
      std::vector<double> m_tab;
      std::vector<double> m_conn;
      std::vector<double> m_rw;
      std::vector<double> m_Kh;
      std::vector<double> m_skinn;
      std::vector<double> m_d_factor;
      std::vector<int> m_copml;
      std::vector<Connection::State> m_state;
      
        void addConnection(const int i, const int j, const int k,
                           const std::size_t global_index,
                           const int complnum,
                           const Connection::State state,
                           const double depth,
                           const Connection::CTFProperties& ctf_props,
                           const int satTableId,
                           const Connection::Direction direction,
                           const Connection::CTFKind ctf_kind,
                           const std::size_t seqIndex,
                           const bool defaultSatTabId);

        size_t findClosestConnection(int oi, int oj, double oz, size_t start_pos);
        void orderTRACK();
        void orderMSW();
        void orderDEPTH();
    };

    std::optional<int>
    getCompletionNumberFromGlobalConnectionIndex(const WellConnections& connections,
                                                 const std::size_t      global_index);
} // namespace Opm

#endif // CONNECTIONSET_HPP_
