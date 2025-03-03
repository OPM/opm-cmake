/*
  Copyright 2020 Equinor ASA.

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

#include <opm/io/eclipse/rst/state.hpp>

#include <opm/io/eclipse/RestartFileView.hpp>
#include <opm/io/eclipse/PaddedOutputString.hpp>

#include <opm/io/eclipse/rst/aquifer.hpp>
#include <opm/io/eclipse/rst/action.hpp>
#include <opm/io/eclipse/rst/connection.hpp>
#include <opm/io/eclipse/rst/group.hpp>
#include <opm/io/eclipse/rst/header.hpp>
#include <opm/io/eclipse/rst/network.hpp>
#include <opm/io/eclipse/rst/udq.hpp>
#include <opm/io/eclipse/rst/well.hpp>

#include <opm/common/utility/String.hpp>
#include <opm/common/utility/TimeService.hpp>

#include <opm/output/eclipse/UDQDims.hpp>

#include <opm/output/eclipse/VectorItems/connection.hpp>
#include <opm/output/eclipse/VectorItems/doubhead.hpp>
#include <opm/output/eclipse/VectorItems/intehead.hpp>
#include <opm/output/eclipse/VectorItems/well.hpp>

#include <opm/output/eclipse/WriteRestartHelpers.hpp>

#include <opm/input/eclipse/Schedule/Action/Actdims.hpp>
#include <opm/input/eclipse/Schedule/Action/Condition.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>

#include <opm/input/eclipse/Parser/Parser.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <functional>
#include <iterator>
#include <memory>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <fmt/format.h>
#include <fmt/ranges.h>

#include <boost/range.hpp>

namespace VI = ::Opm::RestartIO::Helpers::VectorItems;

namespace {
    std::string
    udq_define(const std::vector<std::string>& zudl,
               const std::size_t               udq_index)
    {
        auto zudl_begin = zudl.begin();
        auto zudl_end = zudl_begin;
        std::advance(zudl_begin, (udq_index + 0) * Opm::UDQDims::entriesPerZUDL());
        std::advance(zudl_end  , (udq_index + 1) * Opm::UDQDims::entriesPerZUDL());

        // Note: We inject each ZUDL substring into a field of exactly 8
        // characters ({:<8}) in order to preserve any requisite trailing
        // whitespace.  In particular, we don't want to risk consecutive
        // slots of the form
        //
        //   ")/(GEFF "
        //   "TEST)*(1"
        //
        // being fused into a string fragment of the form
        //
        //   )/(GEFFTEST)*(1
        //
        // This would change the meaning of the summary vector (scalar value
        // 'GEFF TEST' -> group level UDQ set 'GEFFTEST').
        //
        // We do nevertheless use rtrim_copy() on the final string in order
        // to ensure that a UDQ without a defining expression--i.e., an
        // assignment--comes back as an empty string.  The 8 character
        // explicit field width in fmt::format() generates a result string
        // which is always of size 128 regardless of the actual character
        // data (e.g., 128 space characters).  If we were to return that
        // string value unchanged, we would break some of the higher-level
        // logic.
        auto define = Opm::rtrim_copy(fmt::format("{:<8}", fmt::join(zudl_begin, zudl_end, "")));
        if (!define.empty() && (define.front() == '~')) {
            define.front() = '-';
        }

        return define;
    }

    Opm::UDQUpdate udq_update(const std::vector<int>& iudq,
                              const std::size_t       udq_index)
    {
        return Opm::UDQ::updateType(iudq[udq_index * Opm::UDQDims::entriesPerIUDQ()]);
    }

    template <typename T>
    boost::iterator_range<typename std::vector<T>::const_iterator>
    getDataWindow(const std::vector<T>& arr,
                  const std::size_t     windowSize,
                  const std::size_t     entity,
                  const std::size_t     subEntity               = 0,
                  const std::size_t     maxSubEntitiesPerEntity = 1)
    {
        const auto off =
            windowSize * (subEntity + maxSubEntitiesPerEntity*entity);

        auto begin = arr.begin() + off;
        auto end   = begin       + windowSize;

        return { begin, end };
    }

    class UDQVectors
    {
    public:
        explicit UDQVectors(std::shared_ptr<Opm::EclIO::RestartFileView> rst_view)
            : rstView_{ std::move(rst_view) }
        {
            const auto& intehead = this->rstView_->getKeyword<int>("INTEHEAD");

            this->maxNumMsWells_  = intehead[VI::intehead::NSWLMX];
            this->maxNumSegments_ = intehead[VI::intehead::NSEGMX];
            this->numGroups_      = intehead[VI::intehead::NGMAXZ]; // Including FIELD
            this->numWells_       = intehead[VI::intehead::NWMAXZ];
        }

        enum Type : std::size_t {
            Field, Group, Segment, Well,
            // ---- Implementation Helper ----
            NumTypes
        };

        void prepareNext(const Type t) { ++this->varIx_[t]; }

        const std::vector<std::string>& zudn() const
        {
            return this->rstView_->getKeyword<std::string>("ZUDN");
        }

        bool hasGroup()   const { return this->rstView_->hasKeyword<double>("DUDG"); }
        bool hasSegment() const { return this->rstView_->hasKeyword<double>("DUDS"); }
        bool hasWell()    const { return this->rstView_->hasKeyword<double>("DUDW"); }

        double currentFieldUDQValue() const
        {
            return this->rstView_->getKeyword<double>("DUDF")[ this->varIx_[Type::Field] ];
        }

        auto currentGroupUDQValue() const
        {
            return getDataWindow(this->rstView_->getKeyword<double>("DUDG"),
                                 this->numGroups_, this->varIx_[Type::Group]);
        }

        auto currentSegmentUDQValue(const std::size_t msWellIx) const
        {
            return getDataWindow(this->rstView_->getKeyword<double>("DUDS"),
                                 this->maxNumSegments_, this->varIx_[Type::Segment],
                                 msWellIx, this->maxNumMsWells_);
        }

        auto currentWellUDQValue() const
        {
            return getDataWindow(this->rstView_->getKeyword<double>("DUDW"),
                                 this->numWells_, this->varIx_[Type::Well]);
        }

    private:
        std::shared_ptr<Opm::EclIO::RestartFileView> rstView_;

        std::size_t maxNumMsWells_{0};
        std::size_t maxNumSegments_{0};
        std::size_t numGroups_{0};
        std::size_t numWells_{0};

        std::array<std::size_t, Type::NumTypes> varIx_{};
    };

    bool isDefaultedUDQ(const double x)
    {
        return x == Opm::UDQ::restart_default;
    }

    void restoreFieldUDQValue(const UDQVectors&       udqs,
                              Opm::RestartIO::RstUDQ& udq)
    {
        const auto dudf = udqs.currentFieldUDQValue();

        if (! isDefaultedUDQ(dudf)) {
            udq.assignScalarValue(dudf);
        }
    }

    void restoreGroupUDQValue(const UDQVectors&                            udqs,
                              const std::vector<Opm::RestartIO::RstGroup>& groups,
                              Opm::RestartIO::RstUDQ&                      udq)
    {
        const auto nGrp = groups.size();
        const auto dudg = udqs.currentGroupUDQValue();

        const auto subEntity = 0;
        auto entity = 0;
        for (auto iGrp = 0*nGrp; iGrp < nGrp; ++iGrp) {
            if (isDefaultedUDQ(dudg[iGrp])) {
                continue;
            }

            udq.addValue(entity++, subEntity, dudg[iGrp]);
            udq.addEntityName(groups[iGrp].name);
        }
    }

    void restoreSegmentUDQValue(const UDQVectors&                           udqs,
                                const std::vector<Opm::RestartIO::RstWell>& wells,
                                Opm::RestartIO::RstUDQ&                     udq)
    {
        // Counter for MS wells with a non-defaulted UDQ value for at least
        // one segment.
        auto activeMsWellID = 0;

        for (const auto& well : wells) {
            if (well.msw_index == 0) {
                // Not a multi-segmented well.
                continue;
            }

            // Subtract one for zero-based indexing.
            const auto duds = udqs.currentSegmentUDQValue(well.msw_index - 1);
            const auto nSeg = duds.size();

            auto isActiveMsWell = false;
            for (auto iSeg = 0*nSeg; iSeg < nSeg; ++iSeg) {
                if (isDefaultedUDQ(duds[iSeg])) {
                    continue;
                }

                // If we get here, there is at least one non-defaulted SU*
                // value for this well.  Record the well as such, in order
                // to properly associate the entity to a well name.
                isActiveMsWell = true;

                // Note: Use 'iSeg' directly as the sub-entity, as this
                // simplifies ordering.  Trust clients to make this into a
                // one-based segment number when needed.
                udq.addValue(activeMsWellID, iSeg, duds[iSeg]);
            }

            if (isActiveMsWell) {
                // The current well has a non-default UDQ value for at least
                // one segment.  Associate the current entity to the well's
                // name and prepare for handling the next MS well.
                udq.addEntityName(well.name);
                ++activeMsWellID;
            }
        }
    }

    void restoreWellUDQValue(const UDQVectors&                           udqs,
                             const std::vector<Opm::RestartIO::RstWell>& wells,
                             Opm::RestartIO::RstUDQ&                     udq)
    {
        const auto nWell = wells.size();
        const auto dudw  = udqs.currentWellUDQValue();

        const auto subEntity = 0;
        auto entity = 0;
        for (auto iWell = 0*nWell; iWell < nWell; ++iWell) {
            if (isDefaultedUDQ(dudw[iWell])) {
                continue;
            }

            udq.addValue(entity++, subEntity, dudw[iWell]);
            udq.addEntityName(wells[iWell].name);
        }
    }

    void restoreSingleUDQ(const Opm::RestartIO::RstState& rst,
                          UDQVectors&                     udqValues,
                          Opm::RestartIO::RstUDQ&         udq)
    {
        udq.prepareValues();

        // Note: Categories ordered by enumerator values in UDQEnums.hpp.
        switch (udq.category) {
        case Opm::UDQVarType::FIELD_VAR:
            restoreFieldUDQValue(udqValues, udq);
            udqValues.prepareNext(UDQVectors::Type::Field);
            break;

        case Opm::UDQVarType::SEGMENT_VAR:
            restoreSegmentUDQValue(udqValues, rst.wells, udq);
            udqValues.prepareNext(UDQVectors::Type::Segment);
            break;

        case Opm::UDQVarType::WELL_VAR:
            restoreWellUDQValue(udqValues, rst.wells, udq);
            udqValues.prepareNext(UDQVectors::Type::Well);
            break;

        case Opm::UDQVarType::GROUP_VAR:
            restoreGroupUDQValue(udqValues, rst.groups, udq);
            udqValues.prepareNext(UDQVectors::Type::Group);
            break;

        default:
            break;
        }

        udq.commitValues();
    }
}

namespace Opm::RestartIO {

RstState::RstState(std::shared_ptr<EclIO::RestartFileView> rstView,
                   const Runspec&                          runspec,
                   const ::Opm::EclipseGrid*               grid)
    : unit_system(rstView->intehead()[VI::intehead::UNIT])
    , header(runspec, unit_system, rstView->intehead(), rstView->logihead(), rstView->doubhead())
    , aquifers(rstView, grid, unit_system)
    , netbalan(rstView->intehead(), rstView->doubhead(), unit_system)
    , network(rstView, unit_system)
    , oilvap(runspec.tabdims().getNumPVTTables())
{
    this->load_tuning(rstView->intehead(), rstView->doubhead());
    this->load_oil_vaporization(rstView->intehead(), rstView->doubhead());
}

void RstState::load_oil_vaporization(const std::vector<int>& intehead,
                                     const std::vector<double>& doubhead)
{
    const std::size_t numPvtRegions = this->oilvap.numPvtRegions();
    const auto tconv = this->unit_system.to_si(::Opm::UnitSystem::measure::time, 1.0);
    std::vector<double> maximums(numPvtRegions, doubhead[VI::doubhead::dRsDt]/tconv);
    std::vector<std::string> options(numPvtRegions, intehead[VI::intehead::DRSDT_FREE]==1 ? "FREE" : "ALL");
    OilVaporizationProperties::updateDRSDT(this->oilvap, maximums, options);
}

void RstState::load_tuning(const std::vector<int>& intehead,
                           const std::vector<double>& doubhead)
{
    using M  = ::Opm::UnitSystem::measure;

    this->tuning.NEWTMX  = intehead[ VI::intehead::NEWTMX ];
    this->tuning.NEWTMN  = intehead[ VI::intehead::NEWTMN ];
    this->tuning.LITMAX  = intehead[ VI::intehead::LITMAX ];
    this->tuning.LITMIN  = intehead[ VI::intehead::LITMIN ];
    this->tuning.MXWSIT  = intehead[ VI::intehead::MXWSIT ];
    this->tuning.MXWPIT  = intehead[ VI::intehead::MXWPIT ];

    this->tuning.WSEG_MAX_RESTART = intehead[ VI::intehead::WSEGITR_IT2 ];

    {
        const auto tsinit = this->unit_system
            .to_si(M::time, doubhead[VI::doubhead::TsInit]);

        tuning.TSINIT = VI::DoubHeadValue::TSINITHasNoValue(tsinit)
            ? std::nullopt
            : std::optional<double>{ tsinit };
    }

    tuning.TSMAXZ = this->unit_system.to_si(M::time, doubhead[VI::doubhead::TsMaxz]);
    tuning.TSMINZ = this->unit_system.to_si(M::time, doubhead[VI::doubhead::TsMinz]);
    tuning.TSMCHP = this->unit_system.to_si(M::time, doubhead[VI::doubhead::TsMchp]);

    tuning.TSFMAX = doubhead[VI::doubhead::TsFMax];
    tuning.TSFMIN = doubhead[VI::doubhead::TsFMin];
    tuning.TSFCNV = doubhead[VI::doubhead::TsFcnv];
    tuning.THRUPT = doubhead[VI::doubhead::ThrUPT];
    tuning.TFDIFF = doubhead[VI::doubhead::TfDiff];
    tuning.TRGTTE = doubhead[VI::doubhead::TrgTTE];
    tuning.TRGCNV = doubhead[VI::doubhead::TrgCNV];
    tuning.TRGMBE = doubhead[VI::doubhead::TrgMBE];
    tuning.TRGLCV = doubhead[VI::doubhead::TrgLCV];
    tuning.XXXTTE = doubhead[VI::doubhead::XxxTTE];
    tuning.XXXCNV = doubhead[VI::doubhead::XxxCNV];
    tuning.XXXMBE = doubhead[VI::doubhead::XxxMBE];
    tuning.XXXLCV = doubhead[VI::doubhead::XxxLCV];
    tuning.XXXWFL = doubhead[VI::doubhead::XxxWFL];
    tuning.TRGFIP = doubhead[VI::doubhead::TrgFIP];
    tuning.TRGSFT = doubhead[VI::doubhead::TrgSFT];
    tuning.TRGDPR = doubhead[VI::doubhead::TrgDPR];
    tuning.XXXDPR = doubhead[VI::doubhead::XxxDPR];
    tuning.DDPLIM = doubhead[VI::doubhead::DdpLim];
    tuning.DDSLIM = doubhead[VI::doubhead::DdsLim];

    this->tuning.WSEG_REDUCTION_FACTOR = doubhead[VI::doubhead::WsegRedFac];
    this->tuning.WSEG_INCREASE_FACTOR = doubhead[VI::doubhead::WsegIncFac];
}

void RstState::add_groups(const std::vector<std::string>& zgrp,
                          const std::vector<int>& igrp,
                          const std::vector<float>& sgrp,
                          const std::vector<double>& xgrp)
{
    auto load_group = [this, &zgrp, &igrp, &sgrp, &xgrp](const int ig)
    {
        std::size_t zgrp_offset = ig * this->header.nzgrpz;
        std::size_t igrp_offset = ig * this->header.nigrpz;
        std::size_t sgrp_offset = ig * this->header.nsgrpz;
        std::size_t xgrp_offset = ig * this->header.nxgrpz;

        this->groups.emplace_back(this->unit_system,
                                  this->header,
                                  zgrp.data() + zgrp_offset,
                                  igrp.data() + igrp_offset,
                                  sgrp.data() + sgrp_offset,
                                  xgrp.data() + xgrp_offset);
    };

    // Load active named/user-defined groups.
    for (int ig=0; ig < this->header.ngroup; ig++)
        load_group(ig);

    // Load FIELD group from zero-based window index NGMAX in the *GRP
    // arrays.  Needed to reconstruct any field-wide constraints (e.g.,
    // GCONINJE and/or GCONPROD applied to the FIELD itself).
    //
    // Recall that 'max_groups_in_field' is really NGMAX + 1 here as FIELD
    // is also included in this value in the restart file.  Subtract one to
    // get the actual NGMAX value.
    load_group(this->header.max_groups_in_field - 1);
}

void RstState::add_wells(const std::vector<std::string>& zwel,
                         const std::vector<int>& iwel,
                         const std::vector<float>& swel,
                         const std::vector<double>& xwel,
                         const std::vector<int>& icon,
                         const std::vector<float>& scon,
                         const std::vector<double>& xcon)
{

    for (int iw = 0; iw < this->header.num_wells; iw++) {
        std::size_t zwel_offset = iw * this->header.nzwelz;
        std::size_t iwel_offset = iw * this->header.niwelz;
        std::size_t swel_offset = iw * this->header.nswelz;
        std::size_t xwel_offset = iw * this->header.nxwelz;
        std::size_t icon_offset = iw * this->header.niconz * this->header.ncwmax;
        std::size_t scon_offset = iw * this->header.nsconz * this->header.ncwmax;
        std::size_t xcon_offset = iw * this->header.nxconz * this->header.ncwmax;
        int group_index = iwel[ iwel_offset + VI::IWell::Group ] - 1;
        const std::string group = this->groups[group_index].name;

        this->wells.emplace_back(this->unit_system,
                                 this->header,
                                 group,
                                 zwel.data() + zwel_offset,
                                 iwel.data() + iwel_offset,
                                 swel.data() + swel_offset,
                                 xwel.data() + xwel_offset,
                                 icon.data() + icon_offset,
                                 scon.data() + scon_offset,
                                 xcon.data() + xcon_offset);

        if (this->wells.back().msw_index)
            throw std::logic_error("MSW data not accounted for in this constructor");
    }
}

void RstState::add_msw(const std::vector<std::string>& zwel,
                       const std::vector<int>& iwel,
                       const std::vector<float>& swel,
                       const std::vector<double>& xwel,
                       const std::vector<int>& icon,
                       const std::vector<float>& scon,
                       const std::vector<double>& xcon,
                       const std::vector<int>& iseg,
                       const std::vector<double>& rseg)
{

    for (int iw = 0; iw < this->header.num_wells; iw++) {
        std::size_t zwel_offset = iw * this->header.nzwelz;
        std::size_t iwel_offset = iw * this->header.niwelz;
        std::size_t swel_offset = iw * this->header.nswelz;
        std::size_t xwel_offset = iw * this->header.nxwelz;
        std::size_t icon_offset = iw * this->header.niconz * this->header.ncwmax;
        std::size_t scon_offset = iw * this->header.nsconz * this->header.ncwmax;
        std::size_t xcon_offset = iw * this->header.nxconz * this->header.ncwmax;
        int group_index = iwel[ iwel_offset + VI::IWell::Group ] - 1;
        const std::string group = this->groups[group_index].name;

        this->wells.emplace_back(this->unit_system,
                                 this->header,
                                 group,
                                 zwel.data() + zwel_offset,
                                 iwel.data() + iwel_offset,
                                 swel.data() + swel_offset,
                                 xwel.data() + xwel_offset,
                                 icon.data() + icon_offset,
                                 scon.data() + scon_offset,
                                 xcon.data() + xcon_offset,
                                 iseg,
                                 rseg);
    }
}

void RstState::add_udqs(std::shared_ptr<EclIO::RestartFileView> rstView)
{
    if (rstView->hasKeyword<int>("IUAD")) {
        const auto rstFileVersion = rstView->
            getKeyword<int>("INTEHEAD")[VI::intehead::VERSION];

        const auto& iuad = rstView->getKeyword<int>("IUAD");
        const auto& iuap = rstView->getKeyword<int>("IUAP");
        const auto& igph = rstView->getKeyword<int>("IGPH");

        this->udq_active.emplace(rstFileVersion, iuad, iuap, igph);
    }

    const auto& iudq = rstView->getKeyword<int>("IUDQ");
    const auto& zudn = rstView->getKeyword<std::string>("ZUDN");
    const auto& zudl = rstView->getKeyword<std::string>("ZUDL");

    auto udqValues = UDQVectors { std::move(rstView) };

    for (auto udq_index = 0*this->header.num_udq();
         udq_index < this->header.num_udq(); ++udq_index)
    {
        const auto& name = zudn[udq_index*UDQDims::entriesPerZUDN() + 0];
        const auto& unit = zudn[udq_index*UDQDims::entriesPerZUDN() + 1];

        const auto define = udq_define(zudl, udq_index);

        auto& udq = define.empty()
            ? this->udqs.emplace_back(name, unit)
            : this->udqs.emplace_back(name, unit, define, udq_update(iudq, udq_index));

        restoreSingleUDQ(*this, udqValues, udq);
    }
}

void RstState::add_actions(const Parser& parser,
                           const Runspec& runspec,
                           std::time_t sim_time,
                           const std::vector<std::string>& zact,
                           const std::vector<int>& iact,
                           const std::vector<float>& sact,
                           const std::vector<std::string>& zacn,
                           const std::vector<int>& iacn,
                           const std::vector<double>& sacn,
                           const std::vector<std::string>& zlact)
{
    const auto& actdims = runspec.actdims();
    auto zact_action_size  = RestartIO::Helpers::entriesPerZACT();
    auto iact_action_size  = RestartIO::Helpers::entriesPerIACT();
    auto sact_action_size  = RestartIO::Helpers::entriesPerSACT();
    auto zacn_action_size  = RestartIO::Helpers::entriesPerZACN(actdims);
    auto iacn_action_size  = RestartIO::Helpers::entriesPerIACN(actdims);
    auto sacn_action_size  = RestartIO::Helpers::entriesPerSACN(actdims);
    auto zlact_action_size = zlact.size() / this->header.num_action;

    auto zacn_cond_size = 13;
    auto iacn_cond_size = 26;
    auto sacn_cond_size = 16;

    for (std::size_t index=0; index < static_cast<std::size_t>(this->header.num_action); index++) {
        std::vector<RstAction::Condition> conditions;
        for (std::size_t icond = 0; icond < actdims.max_conditions(); icond++) {
            const auto zacn_offset = index * zacn_action_size + icond * zacn_cond_size;
            const auto iacn_offset = index * iacn_action_size + icond * iacn_cond_size;

            if (RstAction::Condition::valid(&zacn[zacn_offset], &iacn[iacn_offset])) {
                const auto sacn_offset = index * sacn_action_size + icond * sacn_cond_size;
                conditions.emplace_back(&zacn[zacn_offset], &iacn[iacn_offset], &sacn[sacn_offset]);
            }
        }

        const auto& name = zact[index * zact_action_size + 0];
        const auto& max_run = iact[index * iact_action_size + 5];
        const auto& run_count = iact[index * iact_action_size + 2] - 1;
        const auto& min_wait = this->unit_system.to_si(UnitSystem::measure::time, sact[index * sact_action_size + 3]);
        const auto& last_run_elapsed = this->unit_system.to_si(UnitSystem::measure::time, sact[index * sact_action_size + 4]);

        auto last_run_time = TimeService::advance( runspec.start_time(), last_run_elapsed );
        this->actions.emplace_back(name, max_run, run_count, min_wait, sim_time, last_run_time, conditions );


        std::string action_deck;
        auto zlact_offset = index * zlact_action_size;
        while (true) {
            std::string line;
            for (std::size_t item_index = 0; item_index < actdims.line_size(); item_index++) {
                const auto padded = EclIO::PaddedOutputString<8> { zlact[zlact_offset + item_index] };
                line += padded.c_str();
            }

            line = trim_copy(line);
            if (line.empty())
                continue;

            if (line == "ENDACTIO")
                break;

            action_deck += line + "\n";
            zlact_offset += actdims.line_size();
        }
        const auto& deck = parser.parseString( action_deck );
        for (auto keyword : deck)
            this->actions.back().keywords.push_back(std::move(keyword));
    }
}

void RstState::add_wlist(const std::vector<std::string>& zwls,
                         const std::vector<int>& iwls)
{
    for (auto well_index = 0*this->header.num_wells; well_index < this->header.num_wells; well_index++) {
        const auto zwls_offset = this->header.max_wlist * well_index;
        const auto iwls_offset = this->header.max_wlist * well_index;
        const auto& well_name = this->wells[well_index].name;

        for (auto wlist_index = 0*this->header.max_wlist; wlist_index < this->header.max_wlist; wlist_index++) {
            const auto well_order = iwls[iwls_offset + wlist_index];
            if (well_order < 1) {
                continue;
            }

            auto& wlist = this->wlists[zwls[zwls_offset + wlist_index]];
            if (wlist.size() < static_cast<std::vector<std::string>::size_type>(well_order)) {
                wlist.resize(well_order);
            }

            wlist[well_order - 1] = well_name;
        }
    }
}

const RstWell& RstState::get_well(const std::string& wname) const
{
    const auto well_iter = std::find_if(this->wells.begin(),
                                        this->wells.end(),
                                        [&wname] (const auto& well) {
                                            return well.name == wname;
                                        });
    if (well_iter == this->wells.end())
        throw std::out_of_range("No such well: " + wname);

    return *well_iter;
}

RstState RstState::load(std::shared_ptr<EclIO::RestartFileView> rstView,
                        const Runspec&                          runspec,
                        const Parser&                           parser,
                        const ::Opm::EclipseGrid*               grid)
{
    RstState state(rstView, runspec, grid);

    // At minimum we need any applicable constraint data for FIELD.  Load
    // groups unconditionally.
    {
        const auto& zgrp = rstView->getKeyword<std::string>("ZGRP");
        const auto& igrp = rstView->getKeyword<int>("IGRP");
        const auto& sgrp = rstView->getKeyword<float>("SGRP");
        const auto& xgrp = rstView->getKeyword<double>("XGRP");

        state.add_groups(zgrp, igrp, sgrp, xgrp);
    }

    if (state.header.num_wells > 0) {
        const auto& zwel = rstView->getKeyword<std::string>("ZWEL");
        const auto& iwel = rstView->getKeyword<int>("IWEL");
        const auto& swel = rstView->getKeyword<float>("SWEL");
        const auto& xwel = rstView->getKeyword<double>("XWEL");

        const auto& icon = rstView->getKeyword<int>("ICON");
        const auto& scon = rstView->getKeyword<float>("SCON");
        const auto& xcon = rstView->getKeyword<double>("XCON");

        if (rstView->hasKeyword<int>("ISEG")) {
            // Multi-segmented wells in restart file.
            const auto& iseg = rstView->getKeyword<int>("ISEG");
            const auto& rseg = rstView->getKeyword<double>("RSEG");

            state.add_msw(zwel, iwel, swel, xwel,
                          icon, scon, xcon,
                          iseg, rseg);
        }
        else {
            // Standard wells only.
            state.add_wells(zwel, iwel, swel, xwel,
                            icon, scon, xcon);
        }

        if (rstView->hasKeyword<int>("IWLS")) {
            const auto& iwls = rstView->getKeyword<int>("IWLS");
            const auto& zwls = rstView->getKeyword<std::string>("ZWLS");

            state.add_wlist(zwls, iwls);
        }
    }

    if (state.header.num_udq() > 0) {
        state.add_udqs(rstView);
    }

    if (state.header.num_action > 0) {
        const auto& zact = rstView->getKeyword<std::string>("ZACT");
        const auto& iact = rstView->getKeyword<int>("IACT");
        const auto& sact = rstView->getKeyword<float>("SACT");
        const auto& zacn = rstView->getKeyword<std::string>("ZACN");
        const auto& iacn = rstView->getKeyword<int>("IACN");
        const auto& sacn = rstView->getKeyword<double>("SACN");
        const auto& zlact= rstView->getKeyword<std::string>("ZLACT");
        state.add_actions(parser, runspec, state.header.sim_time(), zact, iact, sact, zacn, iacn, sacn, zlact);
    }

    return state;
}

} // namespace Opm::RestartIO
