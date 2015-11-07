/************************************************************************/
/*                                                                      */
/*               Copyright 2014-2015 by Ullrich Koethe                  */
/*                                                                      */
/*    This file is part of the VIGRA2 computer vision library.          */
/*    The VIGRA2 Website is                                             */
/*        http://ukoethe.github.io/vigra2                               */
/*    Please direct questions, bug reports, and contributions to        */
/*        ullrich.koethe@iwr.uni-heidelberg.de    or                    */
/*        vigra@informatik.uni-hamburg.de                               */
/*                                                                      */
/*    Permission is hereby granted, free of charge, to any person       */
/*    obtaining a copy of this software and associated documentation    */
/*    files (the "Software"), to deal in the Software without           */
/*    restriction, including without limitation the rights to use,      */
/*    copy, modify, merge, publish, distribute, sublicense, and/or      */
/*    sell copies of the Software, and to permit persons to whom the    */
/*    Software is furnished to do so, subject to the following          */
/*    conditions:                                                       */
/*                                                                      */
/*    The above copyright notice and this permission notice shall be    */
/*    included in all copies or substantial portions of the             */
/*    Software.                                                         */
/*                                                                      */
/*    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    */
/*    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   */
/*    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          */
/*    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       */
/*    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      */
/*    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      */
/*    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     */
/*    OTHER DEALINGS IN THE SOFTWARE.                                   */
/*                                                                      */
/************************************************************************/

#pragma once

#ifndef VIGRA2_SHAPE_HXX_HXX
#define VIGRA2_SHAPE_HXX_HXX

#include "config.hxx"
#include "tinyarray.hxx"
#include <vector>
#include <type_traits>

namespace vigra {

template <class T>
class RangeIter
{
  public:
    RangeIter(T begin, T end, T step)
    : begin_(begin)
    , end_(end)
    , step_(step)
    {
        vigra_precondition(step != 0,
            "RangeIter(): step must be non-zero.");
        vigra_precondition((step > 0 && begin <= end) || (step < 0 && begin >= end),
            "RangeIter(): sign mismatch between step and (end-begin).");
    }

    RangeIter(RangeIter const & other, ReverseCopyTag)
    : RangeIter(other.end_, other.begin_, -other.step_)
    {}

  public:
    RangeIter begin() const
    {
        return *this;
    }

    RangeIter end() const
    {
        return RangeIter(*this, ReverseCopy);
    }

    T const & operator*() const
    {
        return begin_;
    }

    RangeIter & operator++()
    {
        begin_ += step_;
        return *this;
    }

    bool operator!=(RangeIter const & other) const
    {
        return begin_ != other.begin_;
    }

    bool operator<(RangeIter const & other) const
    {
        return (other.begin_ - begin_)*step_ > 0;
    }

    ArrayIndex size() const
    {
        // FIXME: test RangeIter::size() for non-integer ranges
        return floor((abs(end_-begin_+step_)-1)/abs(step_));
    }

  private:
    T begin_, end_, step_;
};

template <class T1, class T2, class T3>
RangeIter<T1>
range(T1 begin, T2 end, T3 step)
{
    return RangeIter<T1>(begin, end, step);
}

template <class T1, class T2>
RangeIter<T1>
range(T1 begin, T2 end)
{
    return RangeIter<T1>(begin, end, 1);
}

template <class T>
RangeIter<T>
range(T end)
{
    return RangeIter<T>(0, end, 1);
}

template <int N=runtime_size>
using Shape = TinyArray<ArrayIndex, N>;

} // namespace vigra

#endif // VIGRA2_SHAPE_HXX_HXX
