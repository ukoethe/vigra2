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

#ifndef VIGRA2_LOOP_HXX_HXX
#define VIGRA2_LOOP_HXX_HXX

#include "config.hxx"
#include "error.hxx"
#include "tinyarray.hxx"
#include "shape.hxx"
#include <vector>
#include <string>
#include <initializer_list>
#include <type_traits>
#include <tuple>
#include <algorithm>

namespace vigra {

class LoopSlice;

class LoopIndex
{
  public:
    typedef ArrayIndex               IndexType;
    typedef TinyArray<IndexType, 3>  Range;

    LoopIndex (IndexType begin, IndexType end, IndexType step=1)
    : id_(this - (LoopIndex *)0)
    , begin_(begin)
    , end_(end)
    , step_(step)
    , has_range_(begin!=end)
    {
        if(has_range_)
            vigra_precondition(
                (step_ > 0 && begin_ < end_) ||
                (step_ < 0 && begin_ > end_),
                "LoopIndex(): sign mismatch between 'step' and '(end-begin)'.");
    }

    LoopIndex (IndexType end=0)
    : LoopIndex(0, end)
    {}

  protected:
    LoopIndex (size_t id, IndexType begin, IndexType end, IndexType step)
    : id_(id)
    , begin_(begin)
    , end_(end)
    , step_(step)
    , has_range_(begin!=end)
    {
        if(has_range_)
            vigra_precondition(
                (step_ > 0 && begin_ < end_) ||
                (step_ < 0 && begin_ > end_),
                "LoopIndex(): sign mismatch between 'step' and '(end-begin)'.");
    }

  public:

    LoopIndex & setRange(LoopIndex const & i)
    {
        begin_     = i.begin_;
        end_       = i.end_;
        step_      = i.step_;
        has_range_ = i.has_range_;
        return *this;
    }

    LoopIndex & setRange(IndexType begin, IndexType end, IndexType step=1)
    {
        begin_     = begin;
        end_       = end;
        step_      = step;
        has_range_ = begin != end;
        if(has_range_)
            vigra_precondition(
                (step_ > 0 && begin_ < end_) ||
                (step_ < 0 && begin_ > end_),
                "LoopIndex::setRange(): sign mismatch between 'step' and '(end-begin)'.");
        return *this;
    }

    LoopIndex & setRange(IndexType end)
    {
        return setRange(0, end);
    }

    LoopSlice slice(IndexType begin, IndexType end, IndexType step=1);

    // LoopSlice sliceStep(IndexType s)
    // {
        // step_ = s;
        // return *this;
    // }

    Range range() const
    {
        return Range(begin_, end_, step_);
    }

    IndexType distance() const
    {
        return end_ - begin_;
    }

    size_t  id_;
    IndexType  begin_, end_, step_;
    bool has_range_;
};

class LoopSlice
: public LoopIndex
{
  public:
    typedef ArrayIndex               IndexType;
    typedef TinyArray<IndexType, 3>  Range;
    friend class LoopIndex;

    LoopSlice(LoopIndex const & i)
    : LoopIndex(i)
    , left_(0)
    , right_(0)
    , has_window_(false)
    {}

  protected:

    LoopSlice(IndexType id, IndexType begin, IndexType end, IndexType step)
    : LoopIndex(id, begin, end, step)
    , left_(0)
    , right_(0)
    , has_window_(false)
    {}

  public:

    LoopSlice & window(IndexType left, IndexType right)
    {
        left_ = left;
        right_ = right;
        has_window_ = true;
        return *this;
    }

    LoopSlice & left(IndexType l)
    {
        return window(l, right_);
    }

    LoopSlice & right(IndexType r)
    {
        return window(left_, r);
    }

    IndexType  left_, right_;
    bool has_window_;
};

inline LoopSlice
LoopIndex::slice(IndexType begin, IndexType end, IndexType step)
{
    return LoopSlice(id_, begin, end, step);
}


template <class T>
class MArrayView
{
  public:
    using SliceVector = std::vector<LoopSlice>;
    using MultiIndex = Shape<>;

    MArrayView(MultiIndex const & shape)
    : shape_(shape)
    , strides_(shapeToStrides(shape_))
    , data_(0)
    {}

    MArrayView(MultiIndex const & shape, MultiIndex const & strides, T * data)
    : shape_(shape)
    , strides_(strides)
    , data_(data)
    {}

    MArrayView(std::initializer_list<ArrayIndex> shape)
    : MArrayView(MultiIndex(shape.begin(), shape.end()))
    {}

    std::tuple<MArrayView, SliceVector>
    operator()(SliceVector indices) const
    {
        vigra_precondition(indices.size() == ndim(),
            "MArrayView::operator(): Expected " + std::to_string(ndim()) +
            " arguments, " + std::to_string(indices.size()) + " given.");
        return std::make_tuple(*this, std::move(indices));
    }

    template <class ... REST>
    std::tuple<MArrayView, SliceVector>
    operator()(LoopSlice i, REST... indices) const
    {
        return operator()({i, indices...});
    }

    std::tuple<MArrayView, SliceVector>
    operator[](SliceVector indices) const
    {
        return operator()(std::move(indices));
    }

    T & operator[](MultiIndex i)
    {
        return data_[dot(i, strides_)];
    }

    T const & operator[](MultiIndex i) const
    {
        return data_[dot(i, strides_)];
    }

    ArrayIndex ndim() const
    {
        return shape_.size();
    }

    ArrayIndex size() const
    {
        return prod(shape_);
    }

    void swap(MArrayView & v)
    {
        std::swap(data_, v.data_);
        shape_.swap(v.shape_);
        strides_.swap(v.strides_);
    }

    MArrayView bind(ArrayIndex dim, ArrayIndex where) const
    {
        T * data = const_cast<T *>(data_);
        return MArrayView(shape_.erase(dim), strides_.erase(dim),
                          data + where*strides_[dim]);
    }

    MultiIndex shape_, strides_;
    T * data_;
};

template <class T>
class MArray
: public MArrayView<T>
{
    using BaseType = MArrayView<T>;

  public:

    using typename BaseType::MultiIndex;
    using Data = std::vector<T>;

    MArray(MultiIndex const & shape)
    : BaseType(shape)
    , alloc_data_(prod(shape), 0)
    {
        this->data_ = alloc_data_.data();
        std::iota(alloc_data_.begin(), alloc_data_.end(), 1);
    }

    MArray(std::initializer_list<ArrayIndex> shape)
    : MArray(MultiIndex(shape.begin(), shape.end()))
    {}

    Data alloc_data_;
};

// template <class T>
// class LoopArray

// FIXME: should NDIM be a template parameter?
class Loop
{
    typedef std::vector<LoopIndex>  IndexVector;
    typedef std::vector<LoopSlice>  SliceVector;
    typedef Shape<>                 MultiIndex;

    template <int N>
    using Level = std::integral_constant<int, N>;

  public:

    friend class LoopTest;

    template <class ... REST>
    Loop(LoopIndex i, REST... indices)
    : Loop{i, indices...}
    {}

    Loop(std::initializer_list<LoopIndex> indices)
    : indices_{indices}
    , nesting_(indices_.size()-1)
    {}

    Loop(IndexVector indices)
    : indices_{std::move(indices)}
    , nesting_(indices_.size()-1)
    {}

    int ndim() const
    {
        return indices_.size();
    }

    Loop & add(std::tuple<MArrayView<int>, SliceVector> const & data)
    {
        MArrayView<int> const & array = std::get<0>(data);
        SliceVector const & slices = std::get<1>(data);
        MArrayView<int> newarray(array);
        MultiIndex strides(ndim());

        Shape<> axes_to_bind(array.ndim());
        Shape<> bind_where(array.ndim());
        for(int i=0; i<ndim(); ++i)
        {
            int j=0;
            for(; j<array.ndim(); ++j)
            {
                if(slices[j].id_ == indices_[i].id_)
                    break;
            }
            if(j < array.ndim()) // index found
            {
                // check range
                if(indices_[i].has_range_)
                {
                    // FIXME: also consider step size
                    vigra_precondition(
                       !slices[j].has_range_  ||
                       (slices[j].distance() == indices_[i].distance()),
                       "Loop::add(): range mismatch in dimension " + std::to_string(i) + ".");
                }
                else if(slices[j].has_range_)
                {
                    indices_[i].setRange(slices[j]);
                }
                else
                {
                    indices_[i].setRange(array.shape_[j]);
                }
                strides[i] = array.strides_[j];
                axes_to_bind[j] = 1;
                bind_where[j] = indices_[i].begin_;
            }
            else
            {
                strides[i] = 0;
            }
        }
        for(int j: range(array.ndim()-1, -1l, -1l))
            if(axes_to_bind[j])
                newarray.bind(j, bind_where[j]).swap(newarray);
        data_.push_back(newarray);
        strides_.push_back(strides);
        return *this;
    }

    template <class FUNC>
    void run(FUNC f) const
    {
        for(int i; i<= nesting_; ++i)
        {
            vigra_precondition(indices_[i].has_range_,
                "Loop::run(): missing range for dimension " + std::to_string(i) + ".");
            if(indices_[i].begin_ == indices_[i].end_)
                return;   // empty range, nothing to do
            vigra_precondition(
                (indices_[i].step_ > 0 && indices_[i].begin_ < indices_[i].end_) ||
                (indices_[i].step_ < 0 && indices_[i].begin_ > indices_[i].end_),
                "Loop::run(): sign mismatch between 'step' and '(end-begin)' in dimension " + std::to_string(i) + ".");
        }
        switch(nesting_)
        {
          // case 0:
          // {
            // Shape<1> index;
            // runImpl(f, index, Level<0>());
            // break;
          // }
          // case 1:
          // {
            // Shape<2> index;
            // runImpl(f, index, Level<0>());
            // break;
          // }
          // case 2:
          // {
            // Shape<3> index;
            // runImpl(f, index, Level<0>());
            // break;
          // }
          default:
          {
            MultiIndex index(ndim());
            runImpl(f, index, 0);
          }
        }
    }

    template<class A1, class A2>
    static void savePointers(A1 & data, A2 & p)
    {
        for(int i=0; i<data.size(); ++i)
            p[i] = data[i].data_;
    }

    template<class A1, class A2>
    static void restorePointers(A1 & data, A2 & p)
    {
        for(int i=0; i<data.size(); ++i)
            data[i].data_ = p[i];
    }

    template<class A1, class A2>
    static void incrementPointers(A1 & data, A2 const & strides, ArrayIndex step, int level)
    {
        for(int i=0; i<data.size(); ++i)
            data[i].data_ += strides[i][level]*step;
    }

    template <class FUNC>
    void runImpl(FUNC f, MultiIndex & index, int level) const
    {
        LoopIndex const & li = indices_[level];
        TinyArray<int *> savedPtrs(data_.size());
        if (level==nesting_)
        {
            if(li.step_ > 0)
            {
                savePointers(data_, savedPtrs);
                for(index[level] =  li.begin_;
                    index[level] <  li.end_;
                    index[level] += li.step_,
                    incrementPointers(data_, strides_, li.step_, level))
                {
                    f(index, data_);
                }
                restorePointers(data_, savedPtrs);
            }
            else
            {
                savePointers(data_, savedPtrs);
                for(index[level] =  li.begin_;
                    index[level] >  li.end_;
                    index[level] += li.step_,
                    incrementPointers(data_, strides_, li.step_, level))
                {
                    f(index, data_);
                }
                restorePointers(data_, savedPtrs);
            }
        }
        else
        {
            if(li.step_ > 0)
            {
                savePointers(data_, savedPtrs);
                for(index[level] =  li.begin_;
                    index[level] <  li.end_;
                    index[level] += li.step_,
                    incrementPointers(data_, strides_, li.step_, level))
                {
                    runImpl(f, index, level+1);
                }
                restorePointers(data_, savedPtrs);
            }
            else
            {
                savePointers(data_, savedPtrs);
                for(index[level] =  li.begin_;
                    index[level] >  li.end_;
                    index[level] += li.step_,
                    incrementPointers(data_, strides_, li.step_, level))
                {
                    runImpl(f, index, level+1);
                }
                restorePointers(data_, savedPtrs);
            }
        }
    }

    template <class FUNC, int N, int LEVEL>
    typename std::enable_if<LEVEL != N-1, void>::type
    runImpl(FUNC f, Shape<N> & index, Level<LEVEL>) const
    {
        LoopIndex const & li = indices_[LEVEL];
        if(li.step_ > 0)
        {
            auto data_ptr = data_[0].data_;
            for(index[LEVEL] =  li.begin_;
                index[LEVEL] <  li.end_;
                index[LEVEL] += li.step_,
                data_[0].data_ += li.step_*strides_[0][LEVEL])
            {
                runImpl(f, index, Level<LEVEL+1>());
            }
            data_[0].data_ = data_ptr;
        }
        else
        {
            auto data_ptr = data_[0].data_;
            for(index[LEVEL] =  li.begin_;
                index[LEVEL] >  li.end_;
                index[LEVEL] += li.step_,
                data_[0].data_ += li.step_*strides_[0][LEVEL])
            {
                runImpl(f, index, Level<LEVEL+1>());
            }
            data_[0].data_ = data_ptr;
        }
    }

    template <class FUNC, int N>
    void runImpl(FUNC f, Shape<N> & index, Level<N-1>) const
    {
        LoopIndex const & li = indices_[N-1];
        if(li.step_ > 0)
        {
            auto data_ptr = data_[0].data_;
            for(index[N-1] =  li.begin_;
                index[N-1] <  li.end_;
                index[N-1] += li.step_,
                data_[0].data_ += li.step_*strides_[0][N-1])
            {
                f(index, data_);
            }
            data_[0].data_ = data_ptr;
        }
        else
        {
            auto data_ptr = data_[0].data_;
            for(index[N-1] =  li.begin_;
                index[N-1] >  li.end_;
                index[N-1] += li.step_,
                data_[0].data_ += li.step_*strides_[0][N-1])
            {
                f(index, data_);
            }
            data_[0].data_ = data_ptr;
        }
    }

    IndexVector indices_;
    mutable std::vector<MArrayView<int> >  data_;
    mutable std::vector<int *>             save_ptrs_;
    mutable std::vector<MultiIndex >       strides_;
    int nesting_;
};

} // namespace vigra

#endif // VIGRA2_LOOP_HXX_HXX
