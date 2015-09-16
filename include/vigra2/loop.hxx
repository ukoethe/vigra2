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

#ifndef VIGRA_LOOP_HXX
#define VIGRA_LOOP_HXX

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
    

// FIXME: Restricted to 2D for now
template <class T>
class MArrayView
{
    typedef std::vector<T>  Data;
    typedef std::vector<LoopSlice>  SliceVector;
    typedef Shape<2> MultiIndex;
    typedef Shape<> TestMultiIndex;
    
  public:
    
    MArrayView(MultiIndex const & shape)
    : shape_(shape)
    , stride_(shapeToStride(shape_))
    , data_(prod(shape_), 0)
    {
        std::iota(data_.begin(), data_.end(), 1);
    }
    
    MArrayView(std::initializer_list<ArrayIndex> shape)
    : MArrayView(MultiIndex(shape.begin()))
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
    
    int ndim() const
    {
        return shape_.size();
    }
    
    MultiIndex shape_, stride_;
    Data data_;
};

// template <class T>
// class LoopArray

// FIXME: should NDIM be a template parameter?
class Loop
{
    typedef std::vector<LoopIndex>  IndexVector;
    typedef std::vector<LoopSlice>  SliceVector;
    typedef std::vector<ArrayIndex> MultiIndex;
    
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
    
    Loop & add(std::tuple<MArrayView<int>, SliceVector> const & data)
    {
        MArrayView<int> const & array = std::get<0>(data);
        SliceVector const & slices = std::get<1>(data);
        
        for(int i; i< array.ndim(); ++i)
        {
            // FIXME: find corresponding index
            
            // check range
            if(indices_[i].has_range_)
            {
                std::cerr << 1 << std::endl;
                // FIXME: also consider step size
                vigra_precondition(
                   !slices[i].has_range_  ||
                   (slices[i].distance() == indices_[i].distance()),
                   "Loop::add(): range mismatch in dimension " + std::to_string(i) + ".");
            }
            else if(slices[i].has_range_)
            {
                std::cerr << 2 << std::endl;
                indices_[i].setRange(slices[i]);
            }
            else
            {
                std::cerr << 3 << std::endl;
                indices_[i].setRange(array.shape_[i]);
            }
        }
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
          case 0:
          {
            Shape<1> index;
            runImpl(f, index, Level<0>());
            break;
          }
          case 1:
          {
            Shape<2> index;
            runImpl(f, index, Level<0>());
            break;
          }
          case 2:
          {
            Shape<3> index;
            runImpl(f, index, Level<0>());
            break;
          }
          default:
          {
            MultiIndex index(indices_.size(), 0);
            runImpl(f, index, 0);
          }
        }
    }
    
    template <class FUNC>
    void runImpl(FUNC f, MultiIndex & index, int level) const
    {
        if (level==nesting_)
        {
            if(indices_[level].step_ > 0)
            {
                for(index[level] = indices_[level].begin_;
                    index[level] < indices_[level].end_;
                    index[level] += indices_[level].step_)
                {
                    f(index);
                }
            }
            else
            {
                for(index[level] = indices_[level].begin_;
                    index[level] > indices_[level].end_;
                    index[level] += indices_[level].step_)
                {
                    f(index);
                }
            }
        }
        else
        {
            if(indices_[level].step_ > 0)
            {
                for(index[level] = indices_[level].begin_;
                    index[level] < indices_[level].end_;
                    index[level] += indices_[level].step_)
                {
                    runImpl(f, index, level+1);
                }
            }
            else
            {
                for(index[level] = indices_[level].begin_;
                    index[level] > indices_[level].end_;
                    index[level] += indices_[level].step_)
                {
                    runImpl(f, index, level+1);
                }
            }
        }
    }
    
    template <class FUNC, int N, int LEVEL>
    typename std::enable_if<LEVEL != N-1, void>::type
    runImpl(FUNC f, Shape<N> & index, Level<LEVEL>) const
    {
        if(indices_[LEVEL].step_ > 0)
        {
            for(index[LEVEL] =  indices_[LEVEL].begin_;
                index[LEVEL] <  indices_[LEVEL].end_;
                index[LEVEL] += indices_[LEVEL].step_)
            {
                runImpl(f, index, Level<LEVEL+1>());
            }
        }
        else
        {
            for(index[LEVEL] =  indices_[LEVEL].begin_;
                index[LEVEL] >  indices_[LEVEL].end_;
                index[LEVEL] += indices_[LEVEL].step_)
            {
                runImpl(f, index, Level<LEVEL+1>());
            }
        }
    }

    template <class FUNC, int N>
    void runImpl(FUNC f, Shape<N> & index, Level<N-1>) const
    {
        if(indices_[N-1].step_ > 0)
        {
            for(index[N-1] =  indices_[N-1].begin_;
                index[N-1] <  indices_[N-1].end_;
                index[N-1] += indices_[N-1].step_)
            {
                f(index);
            }
        }
        else
        {
            for(index[N-1] =  indices_[N-1].begin_;
                index[N-1] >  indices_[N-1].end_;
                index[N-1] += indices_[N-1].step_)
            {
                f(index);
            }
        }
    }
    
    IndexVector indices_;
    std::vector<MArrayView<int> >  data_;
    int nesting_;
};

} // namespace vigra

#endif // VIGRA_LOOP_HXX
