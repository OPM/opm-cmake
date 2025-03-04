/*
  Copyright 2016  Statoil ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify it under the terms
  of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

  OPM is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
  A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along with
  OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <opm/input/eclipse/EclipseState/Runspec.hpp>

#include <opm/input/eclipse/EclipseState/Tables/Regdims.hpp>
#include <opm/input/eclipse/EclipseState/Tables/Tabdims.hpp>

#include <opm/common/utility/TimeService.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Deck/DeckSection.hpp>

#include <opm/input/eclipse/Parser/ParserKeywords/A.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/B.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/C.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/F.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/G.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/H.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/M.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/N.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/O.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/P.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/S.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/T.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/W.hpp>

#include <algorithm>
#include <ctime>
#include <stdexcept>
#include <string>

#include <fmt/format.h>

namespace {

    Opm::Phases inferActivePhases(const Opm::Deck& deck)
    {
        return {
            deck.hasKeyword<Opm::ParserKeywords::OIL>(),
            deck.hasKeyword<Opm::ParserKeywords::GAS>() || deck.hasKeyword<Opm::ParserKeywords::GASWAT>(),
            deck.hasKeyword<Opm::ParserKeywords::WATER>() || deck.hasKeyword<Opm::ParserKeywords::GASWAT>(),
            deck.hasKeyword<Opm::ParserKeywords::SOLVENT>(),
            deck.hasKeyword<Opm::ParserKeywords::POLYMER>(),
            deck.hasKeyword<Opm::ParserKeywords::THERMAL>() || deck.hasKeyword<Opm::ParserKeywords::TEMP>(),
            deck.hasKeyword<Opm::ParserKeywords::POLYMW>(),
            deck.hasKeyword<Opm::ParserKeywords::FOAM>(),
            deck.hasKeyword<Opm::ParserKeywords::BRINE>(),
            deck.hasKeyword<Opm::ParserKeywords::PVTSOL>()
        };
    }

    Opm::SatFuncControls::ThreePhaseOilKrModel
    inferThreePhaseOilKrModel(const Opm::Deck& deck)
    {
        using KroModel = Opm::SatFuncControls::ThreePhaseOilKrModel;

        if (deck.hasKeyword<Opm::ParserKeywords::STONE1>())
            return KroModel::Stone1;

        if (deck.hasKeyword<Opm::ParserKeywords::STONE>() ||
            deck.hasKeyword<Opm::ParserKeywords::STONE2>())
            return KroModel::Stone2;

        return KroModel::Default;
    }

    Opm::SatFuncControls::KeywordFamily
    inferKeywordFamily(const Opm::Deck& deck)
    {
        const auto phases = inferActivePhases(deck);
        const auto wat    = phases.active(Opm::Phase::WATER);
        const auto oil    = phases.active(Opm::Phase::OIL);
        const auto gas    = phases.active(Opm::Phase::GAS);

        const auto threeP = gas && oil && wat;
        const auto twoP = (!gas && oil && wat) || (gas && oil && !wat);

        const auto family1 =       // SGOF/SLGOF and/or SWOF
            (gas && (deck.hasKeyword<Opm::ParserKeywords::SGOF>() ||
                     deck.hasKeyword<Opm::ParserKeywords::SLGOF>())) ||
            (wat && deck.hasKeyword<Opm::ParserKeywords::SWOF>()) ||
            (wat && deck.hasKeyword<Opm::ParserKeywords::SWOFLET>());

        // note: we allow for SOF2 to be part of family1 for threeP +
        // solvent simulations.

        const auto family2 =      // SGFN, SOF{2,3}, SWFN, SGWFN
            (gas && (deck.hasKeyword<Opm::ParserKeywords::SGFN>() ||
                     deck.hasKeyword<Opm::ParserKeywords::SGWFN>())) ||
            (oil && ((threeP && deck.hasKeyword<Opm::ParserKeywords::SOF3>()) ||
                     (twoP && deck.hasKeyword<Opm::ParserKeywords::SOF2>()))) ||
            (wat && (deck.hasKeyword<Opm::ParserKeywords::SWFN>() ||
                     deck.hasKeyword<Opm::ParserKeywords::SGWFN>()));
        const auto family3 = //WSF, GSF gas-water CO2STORE case
            deck.hasKeyword<Opm::ParserKeywords::GSF>() &&
            deck.hasKeyword<Opm::ParserKeywords::WSF>();

        if (family1)
            return Opm::SatFuncControls::KeywordFamily::Family_I;

        if (family2)
            return Opm::SatFuncControls::KeywordFamily::Family_II;

        if (family3)
            return Opm::SatFuncControls::KeywordFamily::Family_III;

        return Opm::SatFuncControls::KeywordFamily::Undefined;
    }

    std::time_t create_start_time(const Opm::Deck& deck)
    {
        if (deck.hasKeyword("START")) {
            const auto& keyword = deck["START"].back();
            return Opm::TimeService::timeFromEclipse(keyword.getRecord(0));
        }
        else {
            // Default to START keyword's default date (1st January 1983).
            return Opm::TimeService::mkdate(1983, 1, 1);
        }
    }

    unsigned long long phaseBit(const bool isSet, const Opm::Phase p)
    {
        return isSet
            ? (1ULL << static_cast<unsigned long long>(p))
            : 0ULL;
    }

}

namespace Opm {

Phases::Phases(bool oil, bool gas, bool wat,
               bool sol, bool pol, bool energy,
               bool polymw,
               bool foam,
               bool brine,
               bool zfraction) noexcept
    : bits { phaseBit(oil       , Phase::OIL)     |
             phaseBit(gas       , Phase::GAS)     |
             phaseBit(wat       , Phase::WATER)   |
             phaseBit(sol       , Phase::SOLVENT) |
             phaseBit(pol       , Phase::POLYMER) |
             phaseBit(energy    , Phase::ENERGY)  |
             phaseBit(polymw    , Phase::POLYMW)  |
             phaseBit(foam      , Phase::FOAM)    |
             phaseBit(brine     , Phase::BRINE)   |
             phaseBit(zfraction , Phase::ZFRACTION) }
{}

Phases Phases::serializationTestObject()
{
    return Phases(true, true, true, false, true, false, true, false);
}

bool Phases::active( Phase p ) const noexcept {
    return this->bits[ static_cast< int >( p ) ];
}

size_t Phases::size() const noexcept {
    return this->bits.count();
}

bool Phases::operator==(const Phases& data) const {
    return bits == data.bits;
}

Welldims::Welldims(const Deck& deck)
{
    using WD = ParserKeywords::WELLDIMS;
    if (deck.hasKeyword<WD>()) {
        const auto& keyword = deck.get<WD>().front();
        const auto& wd = keyword.getRecord(0);

        this->nCWMax = wd.getItem<WD::MAXCONN>().get<int>(0);
        this->nWGMax = wd.getItem<WD::MAX_GROUPSIZE>().get<int>(0);

        // Note: nGMax uses the E100 definition.  E300 instead uses
        //
        //   Max{ "MAXGROUPS", "MAXWELLS" }
        //
        // i.e., the maximum of item 1 and item 4 here.
        this->nGMax = wd.getItem<WD::MAXGROUPS>().get<int>(0);
        this->nWMax = wd.getItem<WD::MAXWELLS>().get<int>(0);

        // Maximum number of well lists pr well.  Always at least 1.
        this->nWlistPrWellMax =
            std::max(WD::MAX_WELLIST_PR_WELL::defaultValue,
                     wd.getItem<WD::MAX_WELLIST_PR_WELL>().get<int>(0));

        // Maximum number of dynamic well lists
        this->nDynWlistMax = wd.getItem<WD::MAX_DYNAMIC_WELLIST>().get<int>(0);


        this->m_location = keyword.location();
    }
}

Welldims Welldims::serializationTestObject()
{
    Welldims result;
    result.nWMax = 1;
    result.nCWMax = 2;
    result.nWGMax = 3;
    result.nGMax = 4;
    result.nWlistPrWellMax = 5;
    result.nDynWlistMax = 6;
    result.m_location = KeywordLocation::serializationTestObject();
    return result;
}

WellSegmentDims::WellSegmentDims()
    : nSegWellMax   { ParserKeywords::WSEGDIMS::NSWLMX::defaultValue }
    , nSegmentMax   { ParserKeywords::WSEGDIMS::NSEGMX::defaultValue }
    , nLatBranchMax { ParserKeywords::WSEGDIMS::NLBRMX::defaultValue }
{}

WellSegmentDims::WellSegmentDims(const Deck& deck)
    : WellSegmentDims {}
{
    using WSD = ParserKeywords::WSEGDIMS;

    if (deck.hasKeyword<WSD>()) {
        const auto& keyword = deck.get<WSD>().front();
        const auto& wsd = keyword.getRecord(0);

        if (const auto& maxMSW = wsd.getItem<WSD::NSWLMX>();
            ! maxMSW.defaultApplied(0))
        {
            this->nSegWellMax = maxMSW.get<int>(0);
        }

        if (const auto& maxSeg = wsd.getItem<WSD::NSEGMX>();
            ! maxSeg.defaultApplied(0))
        {
            this->nSegmentMax = maxSeg.get<int>(0);
        }

        if (const auto& maxBranch = wsd.getItem<WSD::NLBRMX>();
            ! maxBranch.defaultApplied(0))
        {
            this->nLatBranchMax = maxBranch.get<int>(0);
        }

        this->location_ = keyword.location();
    }
}

WellSegmentDims WellSegmentDims::serializationTestObject()
{
    WellSegmentDims result;

    result.nSegWellMax = 1;
    result.nSegmentMax = 2;
    result.nLatBranchMax = 3;
    result.location_ = KeywordLocation::serializationTestObject();

    return result;
}

bool WellSegmentDims::operator==(const WellSegmentDims& data) const
{
    return (this->maxSegmentedWells() == data.maxSegmentedWells())
        && (this->maxSegmentsPerWell() == data.maxSegmentsPerWell())
        && (this->maxLateralBranchesPerWell() == data.maxLateralBranchesPerWell())
        && (this->location() == data.location());
}

NetworkDims::NetworkDims()
    : nMaxNoNodes(0)
    , nMaxNoBranches(0)
    , nMaxNoBranchesConToNode(ParserKeywords::NETWORK::NBCMAX::defaultValue)
{}

NetworkDims::NetworkDims(const Deck& deck)
    : NetworkDims{}
{
    if (deck.hasKeyword<ParserKeywords::NETWORK>()) {
        const auto& wsd = deck.get<ParserKeywords::NETWORK>()[0].getRecord(0);

        this->nMaxNoNodes    = wsd.getItem<ParserKeywords::NETWORK::NODMAX>().get<int>(0);
        this->nMaxNoBranches = wsd.getItem<ParserKeywords::NETWORK::NBRMAX>().get<int>(0);
        this->nMaxNoBranchesConToNode = wsd.getItem<ParserKeywords::NETWORK::NBCMAX>().get<int>(0);

        this->type_ = Type::Extended;
    }
    else if (deck.hasKeyword<ParserKeywords::GRUPNET>()) {
        this->type_ = Type::Standard;
    }
}

NetworkDims NetworkDims::serializationTestObject()
{
    NetworkDims result;
    result.nMaxNoNodes = 1;
    result.nMaxNoBranches = 2;
    result.nMaxNoBranchesConToNode = 3;

    return result;
}

bool NetworkDims::operator==(const NetworkDims& data) const
{
    return this->maxNONodes() == data.maxNONodes() &&
           this->maxNoBranches() == data.maxNoBranches() &&
           this->maxNoBranchesConToNode() == data.maxNoBranchesConToNode();
}

AquiferDimensions::AquiferDimensions()
    : maxNumAnalyticAquifers   { ParserKeywords::AQUDIMS::NANAQU::defaultValue }
    , maxNumAnalyticAquiferConn{ ParserKeywords::AQUDIMS::NCAMAX::defaultValue }
{}

AquiferDimensions::AquiferDimensions(const Deck& deck)
    : AquiferDimensions{}
{
    using AD = ParserKeywords::AQUDIMS;

    if (deck.hasKeyword<AD>()) {
        const auto& keyword = deck.get<AD>().front();
        const auto& ad = keyword.getRecord(0);

        this->maxNumAnalyticAquifers    = ad.getItem<AD::NANAQU>().get<int>(0);
        this->maxNumAnalyticAquiferConn = ad.getItem<AD::NCAMAX>().get<int>(0);
    }
}

AquiferDimensions AquiferDimensions::serializationTestObject()
{
    auto dim = AquiferDimensions{};

    dim.maxNumAnalyticAquifers = 3;
    dim.maxNumAnalyticAquiferConn = 10;

    return dim;
}

bool operator==(const AquiferDimensions& lhs, const AquiferDimensions& rhs)
{
    return (lhs.maxAnalyticAquifers() == rhs.maxAnalyticAquifers())
        && (lhs.maxAnalyticAquiferConnections() == rhs.maxAnalyticAquiferConnections());
}

EclHysterConfig::EclHysterConfig(const Opm::Deck& deck)
    {

        if (!deck.hasKeyword("SATOPTS"))
            return;

        const auto& satoptsItem = deck["SATOPTS"].back().getRecord(0).getItem(0);
        for (unsigned i = 0; i < satoptsItem.data_size(); ++i) {
            std::string satoptsValue = satoptsItem.get< std::string >(0);
            std::transform(satoptsValue.begin(),
                           satoptsValue.end(),
                           satoptsValue.begin(),
                           ::toupper);

            if (satoptsValue == "HYSTER")
                activeHyst = true;
        }

        // check for the (deprecated) HYST keyword
        if (deck.hasKeyword("HYST"))
            activeHyst = true;

        if (!activeHyst)
	      return;

        if (deck.hasKeyword("WAGHYSTR")) {
            if ( !(deck.hasKeyword<Opm::ParserKeywords::OIL>() &&
                   deck.hasKeyword<Opm::ParserKeywords::GAS>() &&
                   deck.hasKeyword<Opm::ParserKeywords::WATER>()) )
                throw std::runtime_error("WAG hysteresis (kw 'WAGHYSTR') requires 'OIL', 'WATER' and 'GAS' present in model. ");

            activeWagHyst = true;
        }

        if (!deck.hasKeyword("EHYSTR")) {
            std::string msg = "Hysteresis is enabled via the HYST parameter for SATOPTS, but the EHYSTR " 
                              "keyword is not present in the deck. \n"
                              "Default values are used for the EHYSTR keyword.";

            OpmLog::warning(msg);
            krHystMod = 0;
            pcHystMod = 0;
            return;
        }
	    /*!
       * \brief Set the type of the hysteresis model which is used for relative permeability.
       *
       * -1: relperm hysteresis is disabled
       * 0: use the Carlson model for relative permeability hysteresis of the non-wetting
       *    phase and the drainage curve for the relperm of the wetting phase
       * 1: use the Carlson model for relative permeability hysteresis of the non-wetting
       *    phase and the imbibition curve for the relperm of the wetting phase
       */
        const auto& ehystrKeyword = deck["EHYSTR"].back();
        std::string whereFlag = ehystrKeyword.getRecord(0).getItem("limiting_hyst_flag").getTrimmedString(0);

        if (deck.hasKeyword("NOHYKR") || whereFlag == "PC")
            krHystMod = -1;
        else {
            krHystMod = ehystrKeyword.getRecord(0).getItem("relative_perm_hyst").get<int>(0);
        }

        // this is slightly screwed: it is possible to specify contradicting hysteresis
        // models with HYPC/NOHYPC and the fifth item of EHYSTR. Let's ignore that for
        // now.
            /*!
       * \brief Return the type of the hysteresis model which is used for capillary pressure.
       *
       * -1: capillary pressure hysteresis is disabled
       * 0: use the Killough model for capillary pressure hysteresis
       */
        if (deck.hasKeyword("NOHYPC") || whereFlag == "KR")
            pcHystMod = -1;
        else {
            // if capillary pressure hysteresis is enabled, Eclipse always uses the
            // Killough model
            pcHystMod = 0;
            curvatureCapPrsValue = ehystrKeyword.getRecord(0).getItem("curvature_capillary_pressure_hyst").get<double>(0);
            if (curvatureCapPrsValue <= 0.0)
                throw std::runtime_error(
                    "Only positive values allowed for the 'capillary pressure curvature parameter' "
                    "(the first item of the 'EHYSTR' keyword).");
        }

        // Killough model: Regularisation for trapped critical saturation calculations
        if (pcHystMod == 0 || krHystMod == 2 || krHystMod == 3)
            modParamTrappedValue = ehystrKeyword.getRecord(0).getItem("mod_param_trapped").get<double>(0);
    }

EclHysterConfig EclHysterConfig::serializationTestObject()
{
    EclHysterConfig result;
    result.activeHyst = true;
    result.pcHystMod = 1;
    result.krHystMod = 2;
    result.modParamTrappedValue = 3;
    result.curvatureCapPrsValue = 4;
    result.activeWagHyst = true;

    return result;
}

bool EclHysterConfig::active() const
{ return activeHyst; }

int EclHysterConfig::pcHysteresisModel() const
{ return pcHystMod; }

int EclHysterConfig::krHysteresisModel() const
{ return krHystMod; }

double EclHysterConfig::modParamTrapped() const
{ return modParamTrappedValue; }

double EclHysterConfig::curvatureCapPrs() const
{ return curvatureCapPrsValue; }

bool EclHysterConfig::activeWag() const
{ return activeWagHyst; }

bool EclHysterConfig::operator==(const EclHysterConfig& data) const
{
    return (this->active() == data.active())
        && (this->pcHysteresisModel() == data.pcHysteresisModel())
        && (this->krHysteresisModel() == data.krHysteresisModel())
        && (this->modParamTrapped() == data.modParamTrapped())
        && (this->curvatureCapPrs() == data.curvatureCapPrs())
        && (this->activeWag() == data.activeWag())
        ;
}

SatFuncControls::SatFuncControls()
    : tolcrit(ParserKeywords::TOLCRIT::VALUE::defaultValue)
{}

SatFuncControls::SatFuncControls(const Deck& deck)
    : SatFuncControls()
{
    using TolCrit = ParserKeywords::TOLCRIT;

    if (deck.hasKeyword<TolCrit>()) {
        // SIDouble doesn't perform any unit conversions here since
        // TOLCRIT is a pure scalar (Dimension = 1).
        this->tolcrit = deck.get<TolCrit>().front().getRecord(0)
            .getItem<TolCrit::VALUE>().getSIDouble(0);
    }

    this->krmodel = inferThreePhaseOilKrModel(deck);
    this->satfunc_family = inferKeywordFamily(deck);
}

SatFuncControls::SatFuncControls(const double tolcritArg,
                                 const ThreePhaseOilKrModel model,
                                 const KeywordFamily family)
    : tolcrit(tolcritArg)
    , krmodel(model)
    , satfunc_family(family)
{}

SatFuncControls SatFuncControls::serializationTestObject()
{
    return SatFuncControls {
        1.0, ThreePhaseOilKrModel::Stone2, KeywordFamily::Family_I
    };
}

bool SatFuncControls::operator==(const SatFuncControls& rhs) const
{
    return (this->minimumRelpermMobilityThreshold() == rhs.minimumRelpermMobilityThreshold())
        && (this->krModel() == rhs.krModel())
        && (this->family() == rhs.family());
}

Nupcol::Nupcol()
    : Nupcol { ParserKeywords::MINNPCOL::VALUE::defaultValue }
{}

Nupcol::Nupcol(int min_value)
    : min_nupcol(min_value)
{
    this->update(ParserKeywords::NUPCOL::NUM_ITER::defaultValue);
}

void Nupcol::update(int value)
{
    if ((value < this->min_nupcol) &&
        (this->min_nupcol == ParserKeywords::MINNPCOL::VALUE::defaultValue))
    {
        OpmLog::note(fmt::format("Minimum NUPCOL value is {} - see keyword "
                                 "MINNPCOL to adjust the minimum value",
                                 this->min_nupcol));
    }

    this->nupcol_value = std::max(value, this->min_nupcol);
}

Nupcol Nupcol::serializationTestObject() {
    Nupcol nc;
    nc.update(123);
    return nc;
}

int Nupcol::value() const
{
    return this->nupcol_value;
}

bool Nupcol::operator==(const Nupcol& data) const {
    return this->min_nupcol == data.min_nupcol &&
           this->nupcol_value == data.nupcol_value;
}


bool Tracers::operator==(const Tracers& other) const {
    return
        this->m_oil_tracers == other.m_oil_tracers &&
        this->m_water_tracers == other.m_water_tracers &&
        this->m_gas_tracers == other.m_gas_tracers &&
        this->m_env_tracers == other.m_env_tracers &&
        this->diffusion_control == other.diffusion_control &&
        this->max_iter == other.max_iter &&
        this->min_iter == other.min_iter;
}

Tracers Tracers::serializationTestObject() {
    Tracers tracers;
    tracers.m_oil_tracers = 123;
    tracers.m_water_tracers = 5;
    tracers.m_gas_tracers = 77;
    tracers.m_env_tracers = 43;
    tracers.diffusion_control = false;
    tracers.max_iter = 11;
    tracers.min_iter = 9;
    return tracers;
}

int Tracers::water_tracers() const {
    return this->m_water_tracers;
}
int Tracers::gas_tracers() const {
    return this->m_gas_tracers;
}
int Tracers::oil_tracers() const {
    return this->m_oil_tracers;
}


Tracers::Tracers(const Deck& deck) {
    using TR = ParserKeywords::TRACERS;

    if (deck.hasKeyword<TR>()) {
        const auto& keyword = deck.get<TR>().back();
        const auto& record = keyword[0];
        this->m_oil_tracers = record.getItem<TR::MAX_OIL_TRACERS>().get<int>(0);
        this->m_water_tracers = record.getItem<TR::MAX_WATER_TRACERS>().get<int>(0);
        this->m_gas_tracers = record.getItem<TR::MAX_GAS_TRACERS>().get<int>(0);
        this->m_env_tracers = record.getItem<TR::MAX_ENV_TRACERS>().get<int>(0);
        this->max_iter = record.getItem<TR::MAX_ITER>().get<int>(0);
        this->min_iter = record.getItem<TR::MIN_ITER>().get<int>(0);

        const auto& diff_control = record.getItem<TR::NUMERIC_DIFF>().get<std::string>(0);
        this->diffusion_control = (diff_control == "DIFF" || diff_control == "SPECIAL");
    } else {
        this->m_oil_tracers = TR::MAX_OIL_TRACERS::defaultValue;
        this->m_water_tracers = TR::MAX_WATER_TRACERS::defaultValue;
        this->m_gas_tracers = TR::MAX_GAS_TRACERS::defaultValue;
        this->m_env_tracers = TR::MAX_ENV_TRACERS::defaultValue;
        this->diffusion_control = false;
        this->max_iter = TR::MAX_ITER::defaultValue;
        this->min_iter = TR::MIN_ITER::defaultValue;
    }
}

const Tracers& Runspec::tracers() const {
    return this->m_tracers;
}

Runspec::Runspec(const Deck& deck)
    : m_start_time (create_start_time(deck))
    , active_phases(inferActivePhases(deck))
    , m_tabdims    (deck)
    , m_regdims    (deck)
    , endscale     (deck)
    , welldims     (deck)
    , wsegdims     (deck)
    , netwrkdims   (deck)
    , aquiferdims  (deck)
    , udq_params   (deck)
    , hystpar      (deck)
    , m_actdims    (deck)
    , m_sfuncctrl  (deck)
    , m_nupcol     ()
    , m_tracers    (deck)
    , m_co2storage (false)
    , m_co2sol     (false)
    , m_h2sol      (false)
    , m_h2storage  (false)
    , m_micp       (false)
    , m_mech       (false)
    , m_temp       (false)
{
    if (DeckSection::hasRUNSPEC(deck)) {
        const RUNSPECSection runspecSection{deck};

        if (runspecSection.hasKeyword<ParserKeywords::MINNPCOL>()) {
            const auto& min_item = runspecSection.get<ParserKeywords::MINNPCOL>()
                .back().getRecord(0).getItem<ParserKeywords::MINNPCOL::VALUE>();

            auto min_value = min_item.get<int>(0);
            this->m_nupcol = Nupcol(min_value);
        }

        using NC = ParserKeywords::NUPCOL;
        if (runspecSection.hasKeyword<NC>()) {
            const auto& item = runspecSection.get<NC>().back()[0].getItem<NC::NUM_ITER>();

            if (item.defaultApplied(0)) {
                const auto msg = fmt::format("OPM Flow uses {} as default "
                                             "NUPCOL value", NC::NUM_ITER::defaultValue);
                OpmLog::note(msg);
            }

            auto deck_nupcol = item.get<int>(0);
            this->m_nupcol.update(deck_nupcol);
        }

        if (runspecSection.hasKeyword<ParserKeywords::CO2STORE>() ||
            runspecSection.hasKeyword<ParserKeywords::CO2STOR>())
        {
            m_co2storage = true;

            if (phases().active(Phase::GAS) &&
                (phases().active(Phase::OIL) || phases().active(Phase::WATER)))
            {
                const std::string msg = "The CO2 storage option is given. PVT properties "
                    "from the Brine-CO2 system is used\n"
                    "See the OPM manual for details on the used models.";

                OpmLog::note(msg);
            }
            else {
                throw std::runtime_error {
                    "The CO2 storage option is given. Activate GAS, plus WATER or OIL"
                };
            }
        }

        if (runspecSection.hasKeyword<ParserKeywords::CO2SOL>()) {
            m_co2sol = true;

            if (phases().active(Phase::SOLVENT)) {
                const std::string msg = "The CO2SOL option is given together "
                    "with SOLVENT. PVT properties from the CO2-Brine system is used\n"
                    "See the OPM manual for details on the used models.";

                OpmLog::note(msg);
            }
            else {
                throw std::runtime_error("The CO2SOL option is given. Activate SOLVENT.");
            }
        }

        if (runspecSection.hasKeyword<ParserKeywords::COMPS>()) {
            const auto& comps_item = runspecSection.get<ParserKeywords::COMPS>().back().getRecord(0).getItem<ParserKeywords::COMPS::NUM_COMPS>();
            const auto num_comps = comps_item.get<int>(0);
            if (num_comps < 1) {
                throw std::logic_error(fmt::format("non-positive COMPS value {} is specified", num_comps));
            }
            this->m_comps = num_comps;
        }

        if (runspecSection.hasKeyword<ParserKeywords::H2SOL>()) {
            m_h2sol = true;

            if (phases().active(Phase::SOLVENT)) {
                const std::string msg = "The H2SOL option is given together with SOLVENT. "
                    "PVT properties from the H2-Brine system is used\n"
                    "See the OPM manual for details on the used models.";

                OpmLog::note(msg);
            }
            else {
                throw std::runtime_error("The H2SOL option is given. Activate SOLVENT.");
            }
        }

        if (runspecSection.hasKeyword<ParserKeywords::H2STORE>()) {
            m_h2storage = true;

            const std::string msg = "The H2 storage option is given. PVT properties "
                "from the Brine-H2 system is used\n"
                "See the OPM manual for details on the used models.";

            OpmLog::note(msg);
        }

        if (runspecSection.hasKeyword<ParserKeywords::MICP>()) {
            m_micp = true;

            const std::string msg = "The MICP option is given. Single phase (WATER) "
                "+ 3 transported components\n + 3 solid phases are used. See "
                "https://doi.org/10.1016/j.ijggc.2021.103256\n"
                "for details on the used model.";

            OpmLog::note(msg);
        }

        if (runspecSection.hasKeyword<ParserKeywords::MECH>()) {
            m_mech = true;

            const std::string msg = "Simulation will solve for mechanical quantities";

            OpmLog::note(msg);
        }

        if (runspecSection.hasKeyword<ParserKeywords::TEMP>()) {
            m_temp = true;

            const std::string msg = "TEMP is used to start a thermal simulation. "
                                    "Conservation equation for internal energy is used instead of enthalpy "
                                    "(i.e., no work terms are included).";

            OpmLog::note(msg);

        }
    }
}

Runspec Runspec::serializationTestObject()
{
    Runspec result;
    result.m_start_time = 1234;
    result.active_phases = Phases::serializationTestObject();
    result.m_tabdims = Tabdims::serializationTestObject();
    result.m_regdims = Regdims::serializationTestObject();
    result.endscale = EndpointScaling::serializationTestObject();
    result.welldims = Welldims::serializationTestObject();
    result.wsegdims = WellSegmentDims::serializationTestObject();
    result.netwrkdims = NetworkDims::serializationTestObject();
    result.aquiferdims = AquiferDimensions::serializationTestObject();
    result.udq_params = UDQParams::serializationTestObject();
    result.hystpar = EclHysterConfig::serializationTestObject();
    result.m_actdims = Actdims::serializationTestObject();
    result.m_sfuncctrl = SatFuncControls::serializationTestObject();
    result.m_nupcol = Nupcol::serializationTestObject();
    result.m_tracers = Tracers::serializationTestObject();
    result.m_comps = 3;
    result.m_co2storage = true;
    result.m_co2sol = true;
    result.m_h2sol = true;
    result.m_h2storage = true;
    result.m_micp = true;
    result.m_mech = true;
    result.m_temp = true;

    return result;
}

const Phases& Runspec::phases() const noexcept {
    return this->active_phases;
}

const Tabdims& Runspec::tabdims() const noexcept {
    return this->m_tabdims;
}

const Regdims& Runspec::regdims() const noexcept {
    return this->m_regdims;
}

const Actdims& Runspec::actdims() const noexcept {
    return this->m_actdims;
}

const EndpointScaling& Runspec::endpointScaling() const noexcept {
    return this->endscale;
}

const Welldims& Runspec::wellDimensions() const noexcept
{
    return this->welldims;
}

const WellSegmentDims& Runspec::wellSegmentDimensions() const noexcept
{
    return this->wsegdims;
}

const NetworkDims& Runspec::networkDimensions() const noexcept
{
    return this->netwrkdims;
}

const AquiferDimensions& Runspec::aquiferDimensions() const noexcept
{
    return this->aquiferdims;
}

const EclHysterConfig& Runspec::hysterPar() const noexcept
{
    return this->hystpar;
}

const SatFuncControls& Runspec::saturationFunctionControls() const noexcept
{
    return this->m_sfuncctrl;
}

const Nupcol& Runspec::nupcol() const noexcept
{
    return this->m_nupcol;
}

bool Runspec::compositionalMode() const
{
    return this->m_comps > 0;
}

size_t Runspec::numComps() const
{
    return this->m_comps;
}

bool Runspec::co2Storage() const noexcept
{
    return this->m_co2storage;
}

bool Runspec::co2Sol() const noexcept
{
    return this->m_co2sol;
}

bool Runspec::h2Sol() const noexcept
{
    return this->m_h2sol;
}

bool Runspec::h2Storage() const noexcept
{
    return this->m_h2storage;
}

bool Runspec::micp() const noexcept
{
    return this->m_micp;
}

bool Runspec::mech() const noexcept
{
    return this->m_mech;
}

bool Runspec::temp() const noexcept
{
    return this->m_temp;
}

bool Runspec::compositional() const noexcept
{
    // Note: co2store and h2store are only in blackoil setting for now
    return (this->m_comps > 0) && !this->m_co2storage && !this->m_h2storage;
}

std::time_t Runspec::start_time() const noexcept
{
    return this->m_start_time;
}

// Returns an integer in the range 0...7 which can be used to indicate
// available phases in Eclipse restart and init files.
int Runspec::eclPhaseMask() const noexcept
{
    const int water = 1 << 2;
    const int oil   = 1 << 0;
    const int gas   = 1 << 1;

    return ( active_phases.active( Phase::WATER ) ? water : 0 )
         | ( active_phases.active( Phase::OIL ) ? oil : 0 )
         | ( active_phases.active( Phase::GAS ) ? gas : 0 );
}

const UDQParams& Runspec::udqParams() const noexcept
{
    return this->udq_params;
}

bool Runspec::rst_cmp(const Runspec& full_spec, const Runspec& rst_spec)
{
    return
        // full_spec.start_time() == rst_spec.start_time() &&  // Can be different between base and restart.
        full_spec.phases() == rst_spec.phases() &&
        full_spec.tabdims() == rst_spec.tabdims() &&
        full_spec.regdims() == rst_spec.regdims() &&
        full_spec.endpointScaling() == rst_spec.endpointScaling() &&
        full_spec.wellSegmentDimensions() == rst_spec.wellSegmentDimensions() &&
        full_spec.networkDimensions() == rst_spec.networkDimensions() &&
        full_spec.aquiferDimensions() == rst_spec.aquiferDimensions() &&
        // full_spec.udqParams() == rst_spec.udqParams() &&   // Can be different between base and restart.
        full_spec.hysterPar() == rst_spec.hysterPar() &&
        full_spec.actdims() == rst_spec.actdims() &&
        full_spec.saturationFunctionControls() == rst_spec.saturationFunctionControls() &&
        full_spec.m_nupcol == rst_spec.m_nupcol &&
        full_spec.m_tracers == rst_spec.m_tracers &&
        full_spec.m_comps == rst_spec.m_comps &&
        full_spec.m_co2storage == rst_spec.m_co2storage &&
        full_spec.m_co2sol == rst_spec.m_co2sol &&
        full_spec.m_h2sol == rst_spec.m_h2sol &&
        full_spec.m_h2storage == rst_spec.m_h2storage &&
        full_spec.m_micp == rst_spec.m_micp &&
        full_spec.m_mech == rst_spec.m_mech &&
        full_spec.m_temp == rst_spec.m_temp &&
        Welldims::rst_cmp(full_spec.wellDimensions(), rst_spec.wellDimensions());
}

bool Runspec::operator==(const Runspec& data) const
{
    return (this->start_time() == data.start_time())
        && (this->phases() == data.phases())
        && (this->tabdims() == data.tabdims())
        && (this->regdims() == data.regdims())
        && (this->endpointScaling() == data.endpointScaling())
        && (this->wellDimensions() == data.wellDimensions())
        && (this->wellSegmentDimensions() == data.wellSegmentDimensions())
        && (this->networkDimensions() == data.networkDimensions())
        && (this->aquiferDimensions() == data.aquiferDimensions())
        && (this->udqParams() == data.udqParams())
        && (this->hysterPar() == data.hysterPar())
        && (this->actdims() == data.actdims())
        && (this->saturationFunctionControls() == data.saturationFunctionControls())
        && (this->m_nupcol == data.m_nupcol)
        && (this->m_tracers == data.m_tracers)
        && (this->m_comps == data.m_comps)
        && (this->m_co2storage == data.m_co2storage)
        && (this->m_co2sol == data.m_co2sol)
        && (this->m_h2sol == data.m_h2sol)
        && (this->m_h2storage == data.m_h2storage)
        && (this->m_micp == data.m_micp)
        && (this->m_mech == data.m_mech)
        && (this->m_temp == data.m_temp)
        ;
}

std::size_t declaredMaxRegionID(const Runspec& rspec)
{
    return std::max(rspec.tabdims().getNumFIPRegions(),
                    rspec.regdims().getNTFIP());
}

} // namespace Opm
