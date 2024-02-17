/*
  Copyright 2023 Equinor.

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

#ifndef WVFPDP_HPP_HEADER_INCLUDED
#define WVFPDP_HPP_HEADER_INCLUDED

namespace Opm {
    class DeckRecord;
} // namespace Opm

namespace Opm { namespace RestartIO {
    struct RstWell;
}} // namespace Opm::RestartIO

namespace Opm {

    class WVFPDP
    {
    public:
        static WVFPDP serializationTestObject();

        double getPressureLoss (double bhp_tab, double thp_limit) const;
        double getPressureAdjustment() const { return m_dp; }
        double getPLossScalingFactor() const { return m_fp; }
        void update(const DeckRecord& record);
        void update(const RestartIO::RstWell& rst_well);

        bool active() const;

        bool operator==(const WVFPDP& other) const;
        bool operator!=(const WVFPDP& other) const;

        template<class Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(m_dp);
            serializer(m_fp);
            serializer(m_active);
        }

    private:
        double m_dp{0.0};
        double m_fp{1.0};
        bool   m_active{false};
    };

} // namespace Opm

#endif  // WVFPDP_HPP_HEADER_INCLUDED
