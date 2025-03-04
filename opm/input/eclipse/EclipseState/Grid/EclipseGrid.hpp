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



#ifndef OPM_PARSER_ECLIPSE_GRID_HPP
#define OPM_PARSER_ECLIPSE_GRID_HPP
#include <opm/input/eclipse/EclipseState/Grid/GridDims.hpp>
#include <opm/input/eclipse/EclipseState/Grid/MapAxes.hpp>
#include <opm/input/eclipse/EclipseState/Grid/MinpvMode.hpp>
#include <opm/input/eclipse/EclipseState/Grid/PinchMode.hpp>

#include <array>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>
#include <map>


namespace Opm {

    class Deck;
    namespace EclIO {
        class EclFile;
        class EclOutput;
    }
    struct NNCdata;
    class UnitSystem;
    class ZcornMapper;
    class EclipseGridLGR;
    class LgrCollection;
    /**
       About cell information and dimension: The actual grid
       information is held in a pointer to an ERT ecl_grid_type
       instance. This pointer must be used for access to all cell
       related properties, including:

         - Size of cells
         - Real world position of cells
         - Active/inactive status of cells
    */

    class EclipseGrid : public GridDims {
    public:
        EclipseGrid();
        explicit EclipseGrid(const std::string& filename);

        /*
          These constructors will make a copy of the src grid, with
          zcorn and or actnum have been adjustments.
        */
        EclipseGrid(const EclipseGrid& src, const std::vector<int>& actnum);
        EclipseGrid(const EclipseGrid& src, const double* zcorn, const std::vector<int>& actnum);

        EclipseGrid(size_t nx, size_t ny, size_t nz,
                    double dx = 1.0, double dy = 1.0, double dz = 1.0,
                    double top = 0.0);
        explicit EclipseGrid(const GridDims& gd);

        EclipseGrid(const std::array<int, 3>& dims ,
                    const std::vector<double>& coord ,
                    const std::vector<double>& zcorn ,
                    const int * actnum = nullptr);

        virtual ~EclipseGrid() = default;

        /// EclipseGrid ignores ACTNUM in Deck, and therefore needs ACTNUM
        /// explicitly.  If a null pointer is passed, every cell is active.
        explicit EclipseGrid(const Deck& deck, const int * actnum = nullptr);

        static bool hasGDFILE(const Deck& deck);
        static bool hasRadialKeywords(const Deck& deck);
        static bool hasSpiderKeywords(const Deck& deck);
        static bool hasCylindricalKeywords(const Deck& deck);
        static bool hasCornerPointKeywords(const Deck&);
        static bool hasCartesianKeywords(const Deck&);
        size_t  getNumActive( ) const;
        bool allActive() const;

        size_t activeIndex(size_t i, size_t j, size_t k) const;
        size_t activeIndex(size_t globalIndex) const;

        size_t getTotalActiveLGR() const;
        size_t getActiveIndexLGR(const std::string& label, size_t i, size_t j, size_t k) const;
        size_t getActiveIndexLGR(const std::string& label, size_t localIndex) const;

        size_t activeIndexLGR(const std::string& label, size_t i, size_t j, size_t k) const;
        size_t activeIndexLGR(const std::string& label, size_t localIndex) const;

        size_t getActiveIndex(size_t i, size_t j, size_t k) const {
            return activeIndex(i, j, k);
        }

        size_t getActiveIndex(size_t globalIndex) const {
            return activeIndex(globalIndex);
        }

        const std::vector<std::string> get_all_lgr_labels() const {
            std::vector<std::string> sliced_labels(all_lgr_labels.begin() + 1, all_lgr_labels.end());
            return sliced_labels;
        }

        const std::vector<std::string> get_all_labels() const {
            return all_lgr_labels;
        }
        
        const std::string get_lgr_tag() const {
            return this->lgr_label;
        }

        std::vector<GridDims> get_lgr_children_gridim() const;
      
        
        void assertIndexLGR(size_t localIndex) const;

        void assertLabelLGR(const std::string& label) const;

        void save(const std::string& filename, bool formatted, const std::vector<Opm::NNCdata>& nnc, const Opm::UnitSystem& units) const;
        /*
          Observe that the there is a getGlobalIndex(i,j,k)
          implementation in the base class. This method - translating
          from an active index to a global index must be implemented
          in the current class.
        */
        void init_children_host_cells(bool logical = true);
        void init_children_host_cells_logical(void);
        void init_children_host_cells_geometrical(void);

        using GridDims::getGlobalIndex;
        size_t getGlobalIndex(size_t active_index) const;

        /*
          For RADIAL grids you can *optionally* use the keyword
          'CIRCLE' to denote that period boundary conditions should be
          applied in the 'THETA' direction; this will only apply if
          the theta keywords entered sum up to exactly 360 degrees!
        */
        bool circle() const;
        bool isPinchActive() const;
        double getPinchThresholdThickness() const;
        PinchMode getPinchOption() const;
        PinchMode getMultzOption() const;
        PinchMode getPinchGapMode() const;
        double getPinchMaxEmptyGap() const;
        MinpvMode getMinpvMode() const;
        const std::vector<double>& getMinpvVector( ) const;

        /*
          Will return a vector of nactive elements. The method will
          behave differently depending on the lenght of the
          input_vector:

             nx*ny*nz: only the values corresponding to active cells
               are copied out.

             nactive: The input vector is copied straight out again.

             ??? : Exception.
        */

        template<typename T>
        std::vector<T> compressedVector(const std::vector<T>& input_vector) const {
            if( input_vector.size() == this->getNumActive() ) {
                return input_vector;
            }

            if (input_vector.size() != getCartesianSize())
                throw std::invalid_argument("Input vector must have full size");

            {
                std::vector<T> compressed_vector( this->getNumActive() );
                const auto& active_map = this->getActiveMap( );

                for (size_t i = 0; i < this->getNumActive(); ++i)
                    compressed_vector[i] = input_vector[ active_map[i] ];

                return compressed_vector;
            }
        }


        /// Will return a vector a length num_active; where the value
        /// of each element is the corresponding global index.
        const std::vector<int>& getActiveMap() const;

        void init_lgr_cells(const LgrCollection& lgr_input);
        void create_lgr_cells_tree(const LgrCollection& );
        /// \brief get cell center, and center and normal of bottom face
        std::tuple<std::array<double, 3>,std::array<double, 3>,std::array<double, 3>>
        getCellAndBottomCenterNormal(size_t globalIndex) const;
        std::array<double, 3> getCellCenter(size_t i,size_t j, size_t k) const;
        std::array<double, 3> getCellCenter(size_t globalIndex) const;
        std::array<double, 3> getCornerPos(size_t i,size_t j, size_t k, size_t corner_index) const;
        const std::vector<double>& activeVolume() const;
        double getCellVolume(size_t globalIndex) const;
        double getCellVolume(size_t i , size_t j , size_t k) const;
        double getCellThickness(size_t globalIndex) const;
        double getCellThickness(size_t i , size_t j , size_t k) const;
        std::array<double, 3> getCellDims(size_t i,size_t j, size_t k) const;
        std::array<double, 3> getCellDims(size_t globalIndex) const;
        bool cellActive( size_t globalIndex ) const;
        bool cellActive( size_t i , size_t j, size_t k ) const;
        bool cellActiveAfterMINPV( size_t i , size_t j , size_t k, double cell_porv ) const;
        bool is_lgr() const {return lgr_grid;};
        std::array<double, 3> getCellDimensions(size_t i, size_t j, size_t k) const {
            return getCellDims(i, j, k);
        }
        std::array<double,3> getCellDimensionsLGR(size_t i, size_t j, size_t k, const std::string& lgr_tag) const;
        double getCellDepthLGR(size_t i, size_t j, size_t k, const std::string& lgr_tag) const;


        bool isCellActive(size_t i, size_t j, size_t k) const {
            return cellActive(i, j, k);
        }

        /// Whether or not given cell has a valid cell geometry
        ///
        /// Valid geometry is defined as all vertices have finite
        /// coordinates and at least one pair of coordinates are separated
        /// by a physical distance along a pillar.
        bool isValidCellGeomtry(const std::size_t globalIndex,
                                const UnitSystem& usys) const;

        double getCellDepth(size_t i,size_t j, size_t k) const;
        double getCellDepth(size_t globalIndex) const;
        ZcornMapper zcornMapper() const;

        const std::vector<double>& getCOORD() const;
        const std::vector<double>& getZCORN() const;
        const std::vector<int>& getACTNUM( ) const;

        const std::optional<MapAxes>& getMapAxes() const;

        const std::map<size_t, std::array<int,2>>& getAquiferCellTabnums() const;

        /*
          The fixupZCORN method is run as part of constructiong the grid. This will adjust the
          z-coordinates to ensure that cells do not overlap. The return value is the number of
          points which have been adjusted. The number of zcorn nodes that has ben fixed is
          stored in private member zcorn_fixed.
        */

        size_t fixupZCORN();
        size_t getZcornFixed() { return zcorn_fixed; };

        // resetACTNUM with no arguments will make all cells in the grid active.

        void resetACTNUM();
        void resetACTNUM( const std::vector<int>& actnum);

        /// \brief Sets MINPVV if MINPV and MINPORV are not used
        void setMINPVV(const std::vector<double>& minpvv);

        bool equal(const EclipseGrid& other) const;
        static bool hasDVDEPTHZKeywords(const Deck&);

        /*
          For ALugrid we can *only* use the keyword <DXV, DXYV, DZV, DEPTHZ> so to
          initialize a Regular Cartesian Grid; further we need equidistant mesh
          spacing in each direction to initialize ALuGrid (mandatory for
          mesh refinement!).
        */

        static bool hasEqualDVDEPTHZ(const Deck&);
        static bool allEqual(const std::vector<double> &v);
        EclipseGridLGR& getLGRCell(std::size_t index);
        EclipseGridLGR& getLGRCell(const std::string& lgr_tag) const;
        int getLGR_global_father(std::size_t global_index,  const std::string& lgr_tag) const;

        std::vector<EclipseGridLGR> lgr_children_cells;
        /**
        * @brief Sets Local Grid Refinement for the EclipseGrid.
        *
        * @param lgr_tag The string that contains the name of a given LGR cell.
        * @param coords The coordinates of a given LGR cell in  CPG COORDSformat.
        * @param zcorn The z-coordinates values of a given LGR cell in CPG ZCORN format.
        */
        virtual void set_lgr_refinement(const std::string& lgr_tag, const std::vector<double>& coords, const std::vector<double>& zcorn);

    protected:
        std::size_t lgr_global_counter = 0;
        std::string lgr_label = "GLOBAL";
        int lgr_level = 0;
        int lgr_level_father = 0;
        std::vector<std::string> lgr_children_labels;
        std::vector<std::size_t> lgr_active_index;
        std::vector<std::size_t> lgr_level_active_map;
        std::vector<std::string> all_lgr_labels;
        std::map<std::vector<std::size_t>, std::size_t> num_lgr_children_cells;
        std::vector<double> m_zcorn;
        std::vector<double> m_coord;
        std::vector<int> m_actnum;
        std::vector<std::size_t> m_print_order_lgr_cells;
        // Input grid data.
        mutable std::optional<std::vector<double>> m_input_zcorn;
        mutable std::optional<std::vector<double>> m_input_coord;

    private:
        std::vector<double> m_minpvVector;
        MinpvMode m_minpvMode;
        std::optional<double> m_pinch;

        // Option 4 of PINCH (TOPBOT/ALL), how to calculate TRANS
        PinchMode m_pinchoutMode;
        // Option 5 of PINCH (TOP/ALL), how to apply MULTZ
        PinchMode m_multzMode;
        // Option 2 of PINCH (GAP/NOGAP)
        PinchMode m_pinchGapMode;
        double    m_pinchMaxEmptyGap;
        bool lgr_grid = false;
        mutable std::optional<std::vector<double>> active_volume;

        bool m_circle = false;
        size_t zcorn_fixed = 0;
        bool m_useActnumFromGdfile = false;

        std::optional<MapAxes> m_mapaxes;

        // Mapping to/from active cells.
        int m_nactive {};
        std::vector<int> m_active_to_global;
        std::vector<int> m_global_to_active;
        // Numerical aquifer cells, needs to be active
        std::unordered_set<size_t> m_aquifer_cells;
        // Keep track of aquifer cell depths and (pvtnum,satnum)
        std::map<size_t, double> m_aquifer_cell_depths;
        std::map<size_t, std::array<int,2>> m_aquifer_cell_tabnums;

        // Radial grids need this for volume calculations.
        std::optional<std::vector<double>> m_thetav;
        std::optional<std::vector<double>> m_rv;
        void parseGlobalReferenceToChildren(void);
        int initializeLGRObjectIndices(int);
        void initializeLGRTreeIndices(void);
        void propagateParentIndicesToLGRChildren(int);
        void updateNumericalAquiferCells(const Deck&);
        double computeCellGeometricDepth(size_t globalIndex) const;

        void initGridFromEGridFile(Opm::EclIO::EclFile& egridfile,
                                   const std::string& fileName);
        void resetACTNUM( const int* actnum);

        void initBinaryGrid(const Deck& deck);

        void initCornerPointGrid(const std::vector<double>& coord ,
                                 const std::vector<double>& zcorn ,
                                 const int * actnum);

        bool keywInputBeforeGdfile(const Deck& deck, const std::string& keyword) const;

        void initCylindricalGrid(const Deck&);
        void initSpiderwebGrid(const Deck&);
        void initSpiderwebOrCylindricalGrid(const Deck&, const bool);
        void initCartesianGrid(const Deck&);
        void initDTOPSGrid(const Deck&);
        void initDVDEPTHZGrid(const Deck&);
        void initGrid(const Deck&, const int* actnum);
        void initCornerPointGrid(const Deck&);
        void assertCornerPointKeywords(const Deck&);
        void save_all_lgr_labels(const LgrCollection& );
        static bool hasDTOPSKeywords(const Deck&);
        static void assertVectorSize(const std::vector<double>& vector, size_t expectedSize, const std::string& msg);

        static std::vector<double> createTOPSVector(const std::array<int, 3>& dims, const std::vector<double>& DZ, const Deck&);
        static std::vector<double> createDVector(const std::array<int, 3>& dims, std::size_t dim, const std::string& DKey, const std::string& DVKey, const Deck&);
        static void scatterDim(const std::array<int, 3>& dims , size_t dim , const std::vector<double>& DV , std::vector<double>& D);


        std::vector<double> makeCoordDxDyDzTops(const std::vector<double>& dx, const std::vector<double>& dy, const std::vector<double>& dz, const std::vector<double>& tops) const;
        std::vector<double> makeZcornDzTops(const std::vector<double>& dz, const std::vector<double>& tops) const ;
        std::vector<double> makeZcornDzvDepthz(const std::vector<double>& dzv, const std::vector<double>& depthz) const;
        std::vector<double> makeCoordDxvDyvDzvDepthz(const std::vector<double>& dxv, const std::vector<double>& dyv, const std::vector<double>& dzv, const std::vector<double>& depthz) const;

        void getCellCorners(const std::array<int, 3>& ijk, const std::array<int, 3>& dims, std::array<double,8>& X, std::array<double,8>& Y, std::array<double,8>& Z) const;
        void getCellCorners(const std::size_t globalIndex,
                            std::array<double,8>& X,
                            std::array<double,8>& Y,
                            std::array<double,8>& Z) const;
        std::array<int,3> getCellSubdivisionRatioLGR(const std::string&  lgr_tag, 
                                                     std::array<int,3>   acum = {1,1,1}) const;
    
    
    };

    /// Specialized Class to describe LGR refined cells.
    class EclipseGridLGR : public EclipseGrid
    {
    public:
        using vec_size_t = std::vector<std::size_t>;
        EclipseGridLGR() = default;
        EclipseGridLGR(const std::string& self_label,
                       const std::string& father_label_,
                       size_t nx,
                       size_t ny,
                       size_t nz,
                       const vec_size_t& father_lgr_index,
                       const std::array<int, 3>& low_fatherIJK_,
                       const std::array<int, 3>& up_fatherIJK_);
        const vec_size_t& getFatherGlobalID() const;
        void save(Opm::EclIO::EclOutput&, const std::vector<Opm::NNCdata>&, const Opm::UnitSystem&) const;
        void save_nnc(Opm::EclIO::EclOutput&) const;
        void set_lgr_global_counter(std::size_t counter)
        {
            lgr_global_counter = counter;
        }
        const vec_size_t& get_father_global() const
        {
            return father_global;
        }
        std::optional<std::reference_wrapper<EclipseGridLGR>>
        get_child_LGR_cell(const std::string& lgr_tag) const;
        std::vector<int> save_hostnum(void) const;
        int get_hostnum(std::size_t global_index) const {return(m_hostnum[global_index]);};
        void get_label_child_to_top_father(std::vector<std::reference_wrapper<const std::string>>& list) const;
        void get_global_index_child_to_top_father(std::vector<std::size_t> & list, 
                                                  std::size_t i, std::size_t j, std::size_t k) const;
        void get_global_index_child_to_top_father(std::vector<std::size_t> & list, std::size_t global_ind) const;
  
        
        
        void set_hostnum(std::vector<int>&);
        const std::array<int,3>& get_low_fatherIJK() const{
          return low_fatherIJK;
        }
        const std::string& get_father_label() const{
          return father_label;
        }

        const std::array<int,3>& get_up_fatherIJK() const{
          return up_fatherIJK;
        }


        /**
         * @brief Sets Local Grid Refinement for the EclipseGridLGR.
         *
         * @param lgr_tag The string that contains the name of a given LGR cell.
         * @param coords The coordinates of a given LGR cell in  CPG COORDSformat.
         * @param zcorn The z-coordinates values of a given LGR cell in CPG ZCORN format.
         */
        void set_lgr_refinement(const std::string& lgr_tag,
                                const std::vector<double>& coord,
                                const std::vector<double>& zcorn) override;

        void set_lgr_refinement(const std::vector<double>&, const std::vector<double>&);

    private:
        void init_father_global();
        std::string father_label;
        // references global on the father label
        vec_size_t father_global;
        std::array<int, 3> low_fatherIJK {};
        std::array<int, 3> up_fatherIJK {};
        std::vector<int> m_hostnum;
    };


    class CoordMapper {
    public:
        CoordMapper(size_t nx, size_t ny);
        size_t size() const;


        /*
          dim = 0,1,2 for x, y and z coordinate respectively.
          layer = 0,1 for k=0 and k=nz layers respectively.
        */

        size_t index(size_t i, size_t j, size_t dim, size_t layer) const;
    private:
        size_t nx;
        size_t ny;
    };



    class ZcornMapper {
    public:
        ZcornMapper(size_t nx, size_t ny, size_t nz);
        size_t index(size_t i, size_t j, size_t k, int c) const;
        size_t index(size_t g, int c) const;
        size_t size() const;

        /*
          The fixupZCORN method will take a zcorn vector as input and
          run through it. If the following situation is detected:

                      /|                     /|
                     / |                    / |
                    /  |                  /   |
                   /   |                 /    |
                  /    |     ==>       /      |
                 /     |             /        |
             ---/------x            /---------x
             | /
             |/

        */
        size_t fixupZCORN( std::vector<double>& zcorn);
        bool validZCORN( const std::vector<double>& zcorn) const;
    private:
        std::array<size_t,3> dims;
        std::array<size_t,3> stride;
        std::array<size_t,8> cell_shift;
    };
}
#endif // OPM_PARSER_ECLIPSE_GRID_HPP
