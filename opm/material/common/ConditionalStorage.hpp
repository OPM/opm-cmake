// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \copydoc Opm::ConditionalStorage
 */
#ifndef OPM_CONDITIONAL_STORAGE_HH
#define OPM_CONDITIONAL_STORAGE_HH

#include "opm/common/ErrorMacros.hpp"
#include <cassert>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include <opm/common/utility/gpuDecorators.hpp>
namespace Opm {
/*!
 * \ingroup Common
 *
 * \brief A simple class which only stores a given member attribute if a boolean
 *        condition is true
 *
 * If the condition is false, nothing is stored and an exception is thrown when trying to
 * access the object.
 */
template <bool cond, class T>
class ConditionalStorage
{
public:
    typedef T type;
    static constexpr bool condition = cond;

    OPM_HOST_DEVICE ConditionalStorage()
    {}

    OPM_HOST_DEVICE explicit ConditionalStorage(const T& v)
        : data_(v)
    {}

    OPM_HOST_DEVICE explicit ConditionalStorage(T&& v)
        : data_(std::move(v))
    {}

    template <class ...Args>
    OPM_HOST_DEVICE ConditionalStorage(Args... args)
        : data_(args...)
    {}

    OPM_HOST_DEVICE ConditionalStorage(const ConditionalStorage& t)
        : data_(t.data_)
    {};

    OPM_HOST_DEVICE ConditionalStorage(ConditionalStorage&& t)
        : data_(std::move(t.data_))
    {};

    OPM_HOST_DEVICE ConditionalStorage& operator=(const ConditionalStorage& v)
    {
        data_ = v.data_;
        return *this;
    }

    OPM_HOST_DEVICE ConditionalStorage& operator=(ConditionalStorage&& v)
    {
        data_ = std::move(v.data_);
        return *this;
    }

    OPM_HOST_DEVICE const T& operator*() const
    { return data_; }
    OPM_HOST_DEVICE T& operator*()
    { return data_; }

    OPM_HOST_DEVICE const T* operator->() const
    { return &data_; }
    OPM_HOST_DEVICE T* operator->()
    { return &data_; }

private:
    T data_{};
};

template <class T>
class ConditionalStorage<false, T>
{
public:
    typedef T type;
    static constexpr bool condition = false;

    OPM_HOST_DEVICE ConditionalStorage()
    {
        static_assert(std::is_default_constructible_v<T>);
    }

    OPM_HOST_DEVICE explicit ConditionalStorage(const T&)
    {
        static_assert(std::is_copy_constructible_v<T>);
    }

    OPM_HOST_DEVICE ConditionalStorage(const ConditionalStorage &)
    {
        // copying an empty conditional storage object does not do anything.
    };

    template <class ...Args>
    OPM_HOST_DEVICE ConditionalStorage(Args...)
    {
        static_assert(std::is_constructible_v<T, Args...>);
    }

    OPM_HOST_DEVICE ConditionalStorage& operator=(const ConditionalStorage&)
    {
        static_assert(std::is_copy_assignable_v<T>);
        return *this;
    }

    OPM_HOST_DEVICE const T& operator*() const
    { 
        OPM_THROW(std::logic_error, "data member deactivated"); 
    }
    OPM_HOST_DEVICE T& operator*()
    { 
        OPM_THROW(std::logic_error, "data member deactivated"); 
    }

    OPM_HOST_DEVICE const T* operator->() const
    { 
        OPM_THROW(std::logic_error, "data member deactivated"); 
    }
    OPM_HOST_DEVICE T* operator->()
    { 
        OPM_THROW(std::logic_error, "data member deactivated"); 
    }

};

} // namespace Opm

#endif
