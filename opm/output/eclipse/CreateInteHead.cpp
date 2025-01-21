/*
  Copyright (c) 2021 Equinor ASA
  Copyright (c) 2018 Statoil ASA

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

#include <opm/output/eclipse/WriteRestartHelpers.hpp>

#include <opm/output/eclipse/InteHEAD.hpp>
#include <opm/output/eclipse/VectorItems/intehead.hpp>

#include <opm/input/eclipse/EclipseState/Aquifer/AquiferConfig.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/EclipseState/Runspec.hpp>
#include <opm/input/eclipse/EclipseState/Tables/Regdims.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>

#include <opm/input/eclipse/Schedule/Action/ActionX.hpp>
#include <opm/input/eclipse/Schedule/Action/Actions.hpp>
#include <opm/input/eclipse/Schedule/ArrayDimChecker.hpp>
#include <opm/input/eclipse/Schedule/GasLiftOpt.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRateConfig.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRateModel.hpp>
#include <opm/input/eclipse/Schedule/Network/Balance.hpp>
#include <opm/input/eclipse/Schedule/Network/ExtNetwork.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Tuning.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQActive.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQConfig.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQEnums.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQInput.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <iterator>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

    using nph_enum = Opm::GuideRateModel::Target;
    const std::map<nph_enum, int> nph_enumToECL = {
        {nph_enum::NONE, 0},
        {nph_enum::OIL,  1},
        {nph_enum::GAS,  3},
        {nph_enum::LIQ,  4},
        {nph_enum::RES,  6},
        {nph_enum::COMB, 9},
    };

    using prod_cmode = Opm::Well::ProducerCMode;
    const std::map<prod_cmode, int> prod_cmodeToECL = {
        {prod_cmode::NONE,  0},
        {prod_cmode::ORAT,  1},
        {prod_cmode::WRAT,  2},
        {prod_cmode::GRAT,  3},
        {prod_cmode::LRAT,  4},
        {prod_cmode::RESV,  5},
        {prod_cmode::BHP,   7},
    };

    int maxConnPerWell(const Opm::Schedule& sched,
                       const std::size_t    report_step,
                       const std::size_t    lookup_step)
    {
        if (report_step == std::size_t{0}) {
            return 0;
        }

        auto ncwmax = 0;
        for (const auto& well : sched.getWells(lookup_step)) {
            const auto ncw = well.getConnections().size();

            ncwmax = std::max(ncwmax, static_cast<int>(ncw));
        }

        return ncwmax;
    }

    int numGroupsInField(const Opm::Schedule& sched,
                         const std::size_t    lookup_step)
    {
        const auto ngmax = sched[lookup_step].groups.size();

        if (ngmax < 1) {
            throw std::invalid_argument {
                "Simulation run must include at least FIELD group"
            };
        }

        // Number of non-FIELD groups.
        return ngmax - 1;
    }

    int GroupControl(const Opm::Schedule& sched,
                     const std::size_t    report_step,
                     const std::size_t    lookup_step)
    {
        int gctrl = 0;
        if (report_step == std::size_t{0}) {
            return gctrl;
        }
        bool have_gconprod = false;
        bool have_gconinje = false;

        for (const auto& group_name : sched.groupNames(lookup_step)) {
            const auto& group = sched.getGroup(group_name, lookup_step);
            if (group.isProductionGroup()) {
                have_gconprod = true;
            }
            if (group.isInjectionGroup()) {
                have_gconinje = true;
            }
        }
        if (have_gconinje)
            gctrl = 2;
        else if (have_gconprod)
            gctrl = 1;

        // Index for group control
        return gctrl;
    }

    int noIuads(const Opm::Schedule& sched,
                const std::size_t    rptStep,
                const std::size_t    simStep)
    {
        if (rptStep == std::size_t{0}) {
            return 0;
        }

        return static_cast<int>
            (sched[simStep].udq_active().iuad().size());
    }

    int noIuaps(const Opm::Schedule& sched,
                const std::size_t    rptStep,
                const std::size_t    simStep)
    {
        if (rptStep == std::size_t{0}) {
            return 0;
        }

        // UDQActive::iuap() returns a vector<> by value.
        const auto iuap = sched[simStep].udq_active().iuap();

        return std::accumulate(iuap.begin(), iuap.end(), 0,
            [](const int n, const auto& rec)
        {
            const auto kw = Opm::UDQ::keyword(rec.control);

            const auto is_field_uda =
                ((kw == Opm::UDAKeyword::GCONPROD) ||
                 (kw == Opm::UDAKeyword::GCONINJE))
                && (rec.wgname == "FIELD");

            // One IUAP entry for each "regular" UDA in WCON* or GCON*.  Two
            // IUAP entries for each field level UDA in GCON*.
            return n + 1 + static_cast<int>(is_field_uda);
        });
    }

    int numMultiSegWells(const ::Opm::Schedule& sched,
                         const std::size_t      report_step,
                         const std::size_t      lookup_step)
    {
        if (report_step == 0) { return 0; }

        const auto& wnames = sched.wellNames(lookup_step);

        return std::count_if(std::begin(wnames), std::end(wnames),
            [&sched, lookup_step](const std::string& wname) -> bool
        {
            return sched.getWell(wname, lookup_step).isMultiSegment();
        });
    }

    int maxNumSegments(const ::Opm::Schedule& sched,
                       const std::size_t      report_step,
                       const std::size_t      lookup_step)
    {
        if (report_step == 0) { return 0; }

        const auto& wnames = sched.wellNames(lookup_step);

        return std::accumulate(std::begin(wnames), std::end(wnames), 0,
            [&sched, lookup_step](const int m, const std::string& wname) -> int
        {
            // maxSegmentID() returns 0 for standard (non-MS) wells.
            return std::max(m, sched.getWell(wname, lookup_step).maxSegmentID());
        });
    }

    int maxNumLateralBranches(const ::Opm::Schedule& sched,
                              const std::size_t      report_step,
                              const std::size_t      lookup_step)
    {
        if (report_step == 0) { return 0; }

        const auto& wnames = sched.wellNames(lookup_step);

        return std::accumulate(std::begin(wnames), std::end(wnames), 0,
            [&sched, lookup_step](const int m, const std::string& wname) -> int
        {
            // maxBranchID() returns 0 for standard (non-MS) wells.
            return std::max(m, sched.getWell(wname, lookup_step).maxBranchID());
        });
    }

    Opm::RestartIO::InteHEAD::WellTableDim
    getWellTableDims(const int              nwgmax,
                     const int              ngmax,
                     const ::Opm::Runspec&  rspec,
                     const ::Opm::Schedule& sched,
                     const std::size_t      report_step,
                     const std::size_t      lookup_step)
    {
        const auto& wd = rspec.wellDimensions();

        const auto numWells = static_cast<int>(sched.numWells(lookup_step));

        const auto maxPerf =
            std::max(wd.maxConnPerWell(),
                     maxConnPerWell(sched, report_step, lookup_step));

        const auto maxWellInGroup =
            std::max(wd.maxWellsPerGroup(), nwgmax);

        const auto maxGroupInField =
            std::max(wd.maxGroupsInField(), ngmax);

        const auto nWMaxz = wd.maxWellsInField();

        return {
            (report_step > 0) ? numWells : 0,
            maxPerf,
            maxWellInGroup,
            maxGroupInField,
            (report_step > 0) ? std::max(nWMaxz, numWells) : nWMaxz,
            wd.maxWellListsPrWell(),
            wd.maxDynamicWellLists()
        };
    }

    std::array<int, 4>
    getNGRPZ(const int             grpsz,
             const int             ngrp,
             const int             num_water_tracer,
             const ::Opm::Runspec& rspec)
    {
        const auto& wd = rspec.wellDimensions();

        const auto nwgmax = std::max(grpsz, wd.maxWellsPerGroup());
        const auto ngmax  = std::max(ngrp , wd.maxGroupsInField());

        const int nigrpz = 97 + std::max(nwgmax, ngmax);
        const int nsgrpz = 112;
        const int nxgrpz = 181 + 4*num_water_tracer;
        const int nzgrpz = 5;

        return {{
            nigrpz,
            nsgrpz,
            nxgrpz,
            nzgrpz,
        }};
    }

    Opm::RestartIO::InteHEAD::Phases
    getActivePhases(const ::Opm::Runspec& rspec)
    {
        auto phases = ::Opm::RestartIO::InteHEAD::Phases{};

        const auto& phasePred = rspec.phases();

        phases.oil   = phasePred.active(Opm::Phase::OIL);
        phases.water = phasePred.active(Opm::Phase::WATER);
        phases.gas   = phasePred.active(Opm::Phase::GAS);

        return phases;
    }

    Opm::RestartIO::InteHEAD::TuningPar
    getTuningPars(const ::Opm::Tuning& tuning)
    {
        return {
            tuning.NEWTMX,
            tuning.NEWTMN,
            tuning.LITMAX,
            tuning.LITMIN,
            tuning.MXWSIT,
            tuning.MXWPIT,
            tuning.WSEG_MAX_RESTART
        };
    }

    Opm::RestartIO::InteHEAD::UdqParam
    getUdqParam(const ::Opm::Runspec& rspec,
                const Opm::Schedule&  sched,
                const std::size_t     rptStep,
                const std::size_t     simStep)
    {
        auto param = Opm::RestartIO::InteHEAD::UdqParam{};

        if (rptStep == std::size_t{0}) {
            return param;
        }

        param.udqParam_1 = rspec.udqParams().rand_seed();

        param.num_iuads = noIuads(sched, rptStep, simStep);
        param.num_iuaps = noIuaps(sched, rptStep, simStep);

        sched[simStep].udq().exportTypeCount(param.numUDQs);

        return param;
    }

    Opm::RestartIO::InteHEAD::ActionParam
    getActionParam(const ::Opm::Runspec&         rspec,
                   const ::Opm::Action::Actions& acts,
                   const std::size_t             rptStep)
    {
        if (rptStep == std::size_t{0}) {
            return { 0, 0, 0, 0 };
        }

        const auto& no_act = acts.ecl_size();
        const auto max_lines_pr_action = acts.max_input_lines();
        const auto max_cond_per_action = rspec.actdims().max_conditions();
        const auto max_characters_per_line = rspec.actdims().max_characters();
        
        return {
            static_cast<int>(no_act),
            max_lines_pr_action,
            static_cast<int>(max_cond_per_action),
            static_cast<int>(max_characters_per_line)
        };
    }


    Opm::RestartIO::InteHEAD::WellSegDims
    getWellSegDims(const int              num_water_tracer,
                   const ::Opm::Runspec&  rspec,
                   const ::Opm::Schedule& sched,
                   const std::size_t      report_step,
                   const std::size_t      lookup_step)
    {
        const auto& wsd = rspec.wellSegmentDimensions();

        const auto numMSW = numMultiSegWells(sched, report_step, lookup_step);
        const auto maxNumSeg = maxNumSegments(sched, report_step, lookup_step);
        const auto maxNumBr = maxNumLateralBranches(sched, report_step, lookup_step);

        return {
            numMSW,
            std::max(numMSW, wsd.maxSegmentedWells()),
            std::max(maxNumSeg, wsd.maxSegmentsPerWell()),
            std::max(maxNumBr, wsd.maxLateralBranchesPerWell()),
            22,           // #ISEG elems per segment
            Opm::RestartIO::InteHEAD::numRsegElem(rspec.phases())
               + 8*num_water_tracer, // #RSEG elems per segment
            10            // #ILBR elems per branch
        };
    }

    Opm::RestartIO::InteHEAD::RegDims
    getRegDims(const ::Opm::TableManager& tdims,
               const ::Opm::Regdims&      rdims)
    {
        const auto ntfip  = tdims.numFIPRegions();
        const auto nmfipr = rdims.getNMFIPR();
        const auto nrfreg = rdims.getNRFREG();
        const auto ntfreg = rdims.getNTFREG();
        const auto nplmix = rdims.getNPLMIX();

        return {
            static_cast<int>(ntfip),
            static_cast<int>(nmfipr),
            static_cast<int>(nrfreg),
            static_cast<int>(ntfreg),
            static_cast<int>(nplmix),
        };
    }

    Opm::RestartIO::InteHEAD::RockOpts
    getRockOpts(const ::Opm::RockConfig& rckCfg, const Opm::Regdims& reg_dims)
    {
        int nttyp  = 1;   // Default value (PVTNUM)
        if (rckCfg.rocknum_property() == "SATNUM") nttyp = 2;
        if (rckCfg.rocknum_property() == "ROCKNUM") nttyp = 4 + reg_dims.getNMFIPR();

        return {
            nttyp
        };
    }

    Opm::RestartIO::InteHEAD::GuideRateNominatedPhase
    setGuideRateNominatedPhase(const ::Opm::Schedule& sched,
                               const std::size_t      report_step,
                               const std::size_t      lookup_step)
    {
        int nom_phase = 0;
        if (report_step == std::size_t{0}) {
            return { nom_phase };
        }

        const auto& guideCFG = sched[lookup_step].guide_rate();
        if (guideCFG.has_model()) {
            const auto& guideRateModel = guideCFG.model();
            
            const auto& targPhase = guideRateModel.target();
            const auto& allow_incr = guideRateModel.allow_increase();
            
            const auto it_nph = nph_enumToECL.find(targPhase);
            if (it_nph != nph_enumToECL.end()) {
                nom_phase = it_nph->second;
            }

            //nominated phase has negative sign for allow increment set to 'NO'
            if (!allow_incr) nom_phase *= -1;
        }

        return {nom_phase};
    }

    int getWhistctlMode(const ::Opm::Schedule& sched,
                        const std::size_t      report_step,
                        const std::size_t      lookup_step)
    {
        int mode = 0;
        if (report_step == std::size_t{0}) {
            return mode;
        }

        const auto& w_hist_ctl_mode = sched.getGlobalWhistctlMmode(lookup_step);
        const auto it_ctl = prod_cmodeToECL.find(w_hist_ctl_mode);
        if (it_ctl != prod_cmodeToECL.end()) {
            mode = it_ctl->second;
        }

        return mode;
    }

    int getLiftOptPar(const ::Opm::Schedule& sched,
                      const std::size_t      report_step,
                      const std::size_t      lookup_step)
    {
        using Value = ::Opm::RestartIO::Helpers::VectorItems::InteheadValues::LiftOpt;

        if (report_step == std::size_t{0}) {
            return Value::NotActive;
        }

        const auto& gasLiftOpt = sched.glo(lookup_step);
        if (! gasLiftOpt.active()) {
            return Value::NotActive;
        }

        return gasLiftOpt.all_newton()
            ? Value::EachNupCol
            : Value::FirstIterationOnly;
    }

    Opm::RestartIO::InteHEAD::ActiveNetwork
    getActiveNetwork(const Opm::Schedule&   sched,
                  const std::size_t      lookup_step)
    {
        const auto&  netwrk = sched[lookup_step].network();
        const auto actntwrk = netwrk.active() ? 2 : 0;
        return {
            actntwrk
        };
    }

    Opm::RestartIO::InteHEAD::NetworkDims
    getNetworkDims(const Opm::Schedule&   sched,
                  const std::size_t      lookup_step,
                  const ::Opm::Runspec& rspec)
    {
        const int noactnod = sched[lookup_step].network().node_names().size();
        const int noactbr  = sched[lookup_step].network().NoOfBranches();
        const int nodmax = std::max(rspec.networkDimensions().maxNONodes(), sched[lookup_step].network().NoOfNodes());
        const int nbrmax = std::max(rspec.networkDimensions().maxNoBranches(), sched[lookup_step].network().NoOfBranches());

        //the following dimensions are fixed
        const int nibran = 14;
        const int nrbran = 11;
        const int ninode = 10;
        const int nrnode = 17;
        const int nznode = 2;
        const int ninobr = 2*nbrmax;

        return {
            noactnod,
            noactbr,
            nodmax,
            nbrmax,
            nibran,
            nrbran,
            ninode,
            nrnode,
            nznode,
            ninobr
        };
    }

    Opm::RestartIO::InteHEAD::NetBalanceDims
    getNetworkBalanceParameters(const Opm::Schedule&   sched,
                  const std::size_t      report_step)
    {
        int maxNoItNBC = 0;
        int maxNoItTHP = 10;
        if (report_step > 0) {
            const auto& sched_state = sched[report_step];
            if (sched_state.network().active()) {
                const auto lookup_step = report_step - 1;
                const auto& netbal = sched[lookup_step].network_balance();
                maxNoItNBC = netbal.pressure_max_iter();
                maxNoItTHP  = netbal.thp_max_iter();
            }
        }
        return {
            maxNoItNBC,
            maxNoItTHP
        };
    }

} // Anonymous

// #####################################################################
// Public Interface (createInteHead()) Below Separator
// ---------------------------------------------------------------------

std::vector<int>
Opm::RestartIO::Helpers::
createInteHead(const EclipseState& es,
               const EclipseGrid&  grid,
               const Schedule&     sched,
               const double        simTime,
               const int           num_solver_steps,
               const int           report_step,
               const int           lookup_step)
{
    const auto nwgmax = (report_step == 0)
        ? 0 : maxGroupSize(sched, lookup_step);

    const auto ngmax  = (report_step == 0)
        ? 0 : numGroupsInField(sched, lookup_step);

    const auto& acts  = sched[lookup_step].actions.get();
    const auto& rspec = es.runspec();
    const auto& tdim  = es.getTableManager();
    const auto& rdim  = tdim.getRegdims();
    const auto& rckcfg = es.getSimulationConfig().rock_config();
    auto num_water_tracer = es.runspec().tracers().water_tracers();
    int nxwelz_tracer_shift = num_water_tracer*5 + 2 * (num_water_tracer > 0);

    const auto ih = InteHEAD{}
        .dimensions         (grid.getNXYZ())
        .numActive          (static_cast<int>(grid.getNumActive()))
        .unitConventions    (es.getDeckUnitSystem())
        .wellTableDimensions(getWellTableDims(nwgmax, ngmax, rspec, sched,
                                              report_step, lookup_step))
        .calendarDate       (getSimulationTimePoint(sched.posixStartTime(), simTime))
        .activePhases       (getActivePhases(rspec))
        .drsdt              (sched, lookup_step)
             // -----------------------------------------------------------------------------------
             //              NIWELZ                | NSWELZ                  | NXWELZ                   | NZWELZ
             //              #IWEL elems per well  | #SWEL elems per well    | #XWEL elems per well     | #ZWEL elems per well
        .params_NWELZ       (155 + num_water_tracer, 122 + 2*num_water_tracer, 131 + nxwelz_tracer_shift, 3)
             // -----------------------------------------------------------------------------------
             //              NICONZ               | NSCONZ               | NXCONZ
             //              #ICON elems per conn | #SCON elems per conn | #XCON elems per conn
        .params_NCON        (26,                    42,                    58 + 5*num_water_tracer)
        .params_GRPZ        (getNGRPZ(nwgmax, ngmax, num_water_tracer, rspec))
        .aquiferDimensions  (inferAquiferDimensions(es, sched[lookup_step]))
        .stepParam          (num_solver_steps, report_step)
        .tuningParam        (getTuningPars(sched[lookup_step].tuning()))
        .liftOptParam       (getLiftOptPar(sched, report_step, lookup_step))
        .wellSegDimensions  (getWellSegDims(num_water_tracer, rspec, sched, report_step, lookup_step))
        .regionDimensions   (getRegDims(tdim, rdim))
        .ngroups            ({ ngmax })
        .params_NGCTRL      (GroupControl(sched, report_step, lookup_step))
        .variousParam       (202204, 100)
        .udqParam_1         (getUdqParam(rspec, sched, report_step, lookup_step))
        .actionParam        (getActionParam(rspec, acts, report_step))
        .variousUDQ_ACTIONXParam()
        .nominatedPhaseGuideRate(setGuideRateNominatedPhase(sched, report_step, lookup_step))
        .whistControlMode   (getWhistctlMode(sched, report_step, lookup_step))
        .activeNetwork      (getActiveNetwork(sched, lookup_step))
        .networkDimensions  (getNetworkDims(sched, lookup_step, rspec))
        .netBalanceData     (getNetworkBalanceParameters(sched, report_step))
        .rockOpts           (getRockOpts(rckcfg, rdim))
        ;

    return ih.data();
}
