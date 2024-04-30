/*
  Copyright 2019 Equinor ASA.

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


#ifndef WELL_PRODUCTION_CONTROLS_HPP
#define WELL_PRODUCTION_CONTROLS_HPP

#include <opm/input/eclipse/Schedule/ScheduleTypes.hpp>
#include <opm/input/eclipse/Schedule/Well/WellEnums.hpp>

#include <cmath>

namespace Opm {

struct WellProductionControls {
public:
    explicit WellProductionControls(int controls_arg) :
        controls(controls_arg)
    {
    }

    bool hasControl(WellProducerCMode cmode_arg) const
    {
        return (this->controls & static_cast<int>(cmode_arg)) != 0;
    }

    void skipControl(WellProducerCMode cmode_arg) {
        auto int_arg = static_cast<int>(cmode_arg);
        if ((this->controls & int_arg) != 0)
            this->controls -= int_arg;
    }

    void addControl(WellProducerCMode cmode_arg) {
        auto int_arg = static_cast<int>(cmode_arg);
        if ((this->controls & int_arg) == 0)
            this->controls += int_arg;
    }

    void clearControls(){
        this->controls = 0;
    }

    bool anyZeroRateConstraint() const {
        auto is_zero = [](const double x)
        {
            return std::isfinite(x) && !std::isnormal(x);
        };

        if (this->hasControl(WellProducerCMode::ORAT) && is_zero(this->oil_rate)) {
            return true;
        }

        if (this->hasControl(WellProducerCMode::WRAT) && is_zero(this->water_rate)) {
            return true;
        }

        if (this->hasControl(WellProducerCMode::GRAT) && is_zero(this->gas_rate)) {
            return true;
        }

        if (this->hasControl(WellProducerCMode::LRAT) && is_zero(this->liquid_rate)) {
            return true;
        }

        if (this->hasControl(WellProducerCMode::RESV)) {
            // currently, the RESV constraint for a history matching producers will be calculated based on the reservoir condition
            // at the moment it is used, and we do not have a fixed resv_rate to check against zero. as a result,
            // we check whether the historic rates for oil, gas and water are all zero to determine whether it is under zero
            // RESV constraint for history matching producers.
            const auto zero_rate = this->prediction_mode
                ? is_zero(this->resv_rate)
                : is_zero(this->oil_rate) && is_zero(this->water_rate) && is_zero(this->gas_rate);

            if (zero_rate) {
                return true;
            }
        }

        return false;
    }

    bool operator==(const WellProductionControls& other) const
    {
        return this->cmode == other.cmode &&
               this->oil_rate == other.oil_rate &&
               this->water_rate == other.water_rate &&
               this->gas_rate == other.gas_rate &&
               this->liquid_rate == other.liquid_rate &&
               this->resv_rate == other.resv_rate &&
               this->bhp_history == other.bhp_history &&
               this->thp_history == other.thp_history &&
               this->bhp_limit == other.bhp_limit &&
               this->thp_limit == other.thp_limit &&
               this->alq_value == other.alq_value &&
               this->vfp_table_number == other.vfp_table_number &&
               this->prediction_mode == other.prediction_mode;
    }

    WellProducerCMode cmode = WellProducerCMode::NONE;
    double oil_rate{0};
    double water_rate{0};
    double gas_rate{0};
    double liquid_rate{0};
    double resv_rate{0};
    double bhp_history{0};
    double thp_history{0};
    double bhp_limit{0};
    double thp_limit{0};
    double alq_value{0};
    int    vfp_table_number{0};
    bool   prediction_mode{0};

private:
    int controls;
};

}

#endif
