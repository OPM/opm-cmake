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

#ifndef ACTIONX_CONDITION_HPP
#define ACTIONX_CONDITION_HPP

#include <string>
#include <vector>

#include <opm/common/OpmLog/KeywordLocation.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Action/Enums.hpp>

namespace Opm {

namespace Action {


class Quantity {
public:
    Quantity() = default;

    Quantity(const std::string& arg) :
        quantity(arg)
    {}

    void add_arg(const std::string& arg);
    std::string quantity;
    std::vector<std::string> args;
    bool date() const;

    bool operator==(const Quantity& data) const {
        return quantity == data.quantity &&
               args == data.args;
    }

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(quantity);
        serializer(args);
    }
};



class Condition {
public:

    Condition() = default;
    Condition(const std::vector<std::string>& tokens, const KeywordLocation& location);


    Quantity lhs;
    Quantity rhs;
    Logical logic = Logical::END;
    Comparator cmp = Comparator::INVALID;
    bool left_paren = false;
    bool right_paren = false;

    std::string cmp_string;

    static Logical logic_from_int(int);
    int logic_as_int() const;
    static Comparator comparator_from_int(int);
    int comparator_as_int() const;

    bool open_paren() const;
    bool close_paren() const;
    bool operator==(const Condition& data) const;

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        lhs.serializeOp(serializer);
        rhs.serializeOp(serializer);
        serializer(logic);
        serializer(cmp);
        serializer(cmp_string);
        serializer(left_paren);
        serializer(right_paren);
    }
};


}
}

#endif
