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

#ifndef VIGRA_TINYARRAY_DYNAMIC_HXX
#define VIGRA_TINYARRAY_DYNAMIC_HXX

#include "config.hxx"
#include "tinyarray.hxx"
#include <memory>

namespace vigra {

template <class VALUETYPE, class DERIVED>
class TinyArrayBase<VALUETYPE, DERIVED, 0>
{
  public:
  
    template <class NEW_VALUETYPE>
    using AsType = TinyArray<NEW_VALUETYPE, 0>;
    
    using value_type             = VALUETYPE;
    using reference              = value_type &;
    using const_reference        = value_type const &;
    using pointer                = value_type *;
    using const_pointer          = value_type const *;
    using iterator               = value_type *;
    using const_iterator         = value_type const *;
    using reverse_iterator       = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;
    using size_type              = std::size_t;
    using difference_type        = std::ptrdiff_t;
    using index_type             = ArrayIndex;
    
    static const ArrayIndex static_ndim  = 1;
    
  protected:
  
    template <int LEVEL, class ... V2>
    void initImpl(VALUETYPE v1, V2... v2)
    {
        data_[LEVEL] = v1;
        initImpl<LEVEL+1>(v2...);
    }
    
    template <int LEVEL>
    void initImpl(VALUETYPE v1)
    {
        data_[LEVEL] = v1;
    }
    
    TinyArrayBase(ArrayIndex size=0, pointer data=0)
    : size_(size)
    , data_(data)
    {}

  public:
    
    TinyArrayBase(TinyArrayBase const &) = default;
    
    // assignment
    
    TinyArrayBase & operator=(value_type v)
    {
        for(int i=0; i<size_; ++i)
            data_[i] = v;
        return *this;
    }
    
    TinyArrayBase & operator=(TinyArrayBase const & rhs)
    {
        vigra_precondition(size_ == rhs.size(),
            "TinyArrayBase::operator=(): size mismatch.");
        for(int i=0; i<size_; ++i)
            data_[i] = rhs[i];
        return *this;
    }
    
    template <class OTHER, class OTHER_DERIVED, int N>
    TinyArrayBase & operator=(TinyArrayBase<OTHER, OTHER_DERIVED, N> const & other)
    {
        vigra_precondition(size_ == other.size(),
            "TinyArrayBase::operator=(): size mismatch.");
        for(int i=0; i<size_; ++i)
            data_[i] = detail::RequiresExplicitCast<value_type>::cast(other[i]);
        return *this;
    }
    
    DERIVED & init(value_type v = value_type())
    {
        for(int i=0; i<size_; ++i)
            data_[i] = v;
        return static_cast<DERIVED &>(*this);
    }
    
    template <class ... V>
    DERIVED & init(value_type v0, value_type v1, V... v)
    {
        vigra_precondition(sizeof...(V)+2 == size_, 
                      "TinyArrayBase::init(): wrong number of arguments.");
        initImpl<0>(v0, v1, v...);
        return static_cast<DERIVED &>(*this);
    }
    
    template <class Iterator>
    DERIVED & init(Iterator first, Iterator end)
    {
        int range = std::distance(first, end);
        if(size_ < range)
            range = size_;
        for(int i=0; i<range; ++i, ++first)
            data_[i] = detail::RequiresExplicitCast<value_type>::cast(*first);
        return static_cast<DERIVED &>(*this);
    }
    
    // arithmetic assignment
    
    DERIVED & operator+=(value_type v)
    {
        for(int i=0; i<size_; ++i)
            data_[i] += v;
        return static_cast<DERIVED &>(*this);
    }
    
    template <class OTHER, class OTHER_DERIVED, int N>
    DERIVED & operator+=(TinyArrayBase<OTHER, OTHER_DERIVED, N> const & other)
    {
        vigra_precondition(size_ == other.size(),
            "TinyArrayBase::operator+=(): size mismatch.");
        for(int i=0; i<size_; ++i)
            data_[i] += other[i];
        return static_cast<DERIVED &>(*this);
    }
    
    DERIVED & operator-=(value_type v)
    {
        for(int i=0; i<size_; ++i)
            data_[i] -= v;
        return static_cast<DERIVED &>(*this);
    }
    
    template <class OTHER, class OTHER_DERIVED, int N>
    DERIVED & operator-=(TinyArrayBase<OTHER, OTHER_DERIVED, N> const & other)
    {
        vigra_precondition(size_ == other.size(),
            "TinyArrayBase::operator-=(): size mismatch.");
        for(int i=0; i<size_; ++i)
            data_[i] -= other[i];
        return static_cast<DERIVED &>(*this);
    }
    
    DERIVED & operator*=(value_type v)
    {
        for(int i=0; i<size_; ++i)
            data_[i] *= v;
        return static_cast<DERIVED &>(*this);
    }
    
    template <class OTHER, class OTHER_DERIVED, int N>
    DERIVED & operator*=(TinyArrayBase<OTHER, OTHER_DERIVED, N> const & other)
    {
        vigra_precondition(size_ == other.size(),
            "TinyArrayBase::operator*=(): size mismatch.");
        for(int i=0; i<size_; ++i)
            data_[i] *= other[i];
        return static_cast<DERIVED &>(*this);
    }
    
    DERIVED & operator/=(value_type v)
    {
        for(int i=0; i<size_; ++i)
            data_[i] /= v;
        return static_cast<DERIVED &>(*this);
    }
    
    template <class OTHER, class OTHER_DERIVED, int N>
    DERIVED & operator/=(TinyArrayBase<OTHER, OTHER_DERIVED, N> const & other)
    {
        vigra_precondition(size_ == other.size(),
            "TinyArrayBase::operator/=(): size mismatch.");
        for(int i=0; i<size_; ++i)
            data_[i] /= other[i];
        return static_cast<DERIVED &>(*this);
    }
    
    DERIVED & operator%=(value_type v)
    {
        for(int i=0; i<size_; ++i)
            data_[i] %= v;
        return static_cast<DERIVED &>(*this);
    }
    
    template <class OTHER, class OTHER_DERIVED, int N>
    DERIVED & operator%=(TinyArrayBase<OTHER, OTHER_DERIVED, N> const & other)
    {
        vigra_precondition(size_ == other.size(),
            "TinyArrayBase::operator%=(): size mismatch.");
        for(int i=0; i<size_; ++i)
            data_[i] %= other[i];
        return static_cast<DERIVED &>(*this);
    }
    
    // vector functions
    
    SquaredNormType<value_type> squaredNormImpl() const
    {
        SquaredNormType<value_type> result = SquaredNormType<value_type>();
        for(int i=0; i<size_; ++i)
            result += squaredNorm(data_[i]);
        return result;
    }

        /** Return the minimal element.
        */
    const_reference minimum() const
    {
        int m = 0;
        for(int i=1; i<size_; ++i)
            if(data_[i] < data_[m])
                m = i;
        return data_[m];
    }

        /** Return the maximal element.
        */
    const_reference maximum() const
    {
        int m = 0;
        for(int i=1; i<size_; ++i)
            if(data_[m] < data_[i])
                m = i;
        return data_[m];
    }
    
        /** Check that all elements of this vector are zero (or 'false' if T is bool).
        */
    bool allZero() const
    {
        for(int i=0; i<size_; ++i)
            if(data_[i] != value_type())
                return false;
        return true;
    }
    
        /** Check that all elements of this vector are non-zero (or 'true' if T is bool).
        */
    bool all() const
    {
        for(int i=0; i<size_; ++i)
            if(data_[i] == value_type())
                return false;
        return true;
    }
    
        /** Check that at least one element of this vector is non-zero (or 'true' if T is bool).
        */
    bool any() const
    {
        for(int i=0; i<size_; ++i)
            if(data_[i] != value_type())
                return true;
        return false;
    }
    
    // index access
    
    reference operator[](ArrayIndex i)
    {
        return data_[i];
    }
    
    constexpr const_reference operator[](ArrayIndex i) const
    {
        return data_[i];
    }
    
    reference at(ArrayIndex i)
    {
        if(i < 0 || i >= size_)
            throw std::out_of_range("TinyArrayBase::at()");
        return data_[i];
    }
    
    const_reference at(ArrayIndex i) const
    {
        if(i < 0 || i >= size_)
            throw std::out_of_range("TinyArrayBase::at()");
        return data_[i];
    }
        
        /** Get a view to the subarray with length <tt>(TO-FROM)</tt> starting at <tt>FROM</tt>.
            The bounds must fullfill <tt>0 <= FROM < TO <= size_</tt>.
        */
    template <int FROM, int TO>
    TinyArrayBase
    subarray() const
    {
        vigra_precondition(FROM >= 0 && FROM < TO && TO <= size_, 
                      "TinyArrayBase::subarray(): range out of bounds.");
        return TinyArrayBase(TO-FROM, data_+FROM);
    }
    
    TinyArray<value_type, 0>
    dropIndex(ArrayIndex m) const
    {
        TinyArray<value_type, 0> res(size_-1, DontInit);
        for(int k=0; k<m; ++k)
            res[k] = data_[k];
        for(int k=m; k<size_-1; ++k)
            res[k] = data_[k+1];
        return res;
    }

    // boiler plate
    
    iterator begin() { return data_; }
    iterator end()   { return data_ + size_; }
    const_iterator begin() const { return data_; }
    const_iterator end()   const { return data_ + size_; }
    const_iterator cbegin() const { return data_; }
    const_iterator cend()   const { return data_ + size_; }
    
    reverse_iterator rbegin() { return reverse_iterator(data_ + size_); }
    reverse_iterator rend()   { return reverse_iterator(data_); }
    const_reverse_iterator rbegin() const { return const_reverse_iterator(data_ + size_); }
    const_reverse_iterator rend()   const { return const_reverse_iterator(data_); }
    const_reverse_iterator crbegin() const { return const_reverse_iterator(data_ + size_); }
    const_reverse_iterator crend()   const { return const_reverse_iterator(data_); }
    
    pointer data() { return data_; }
    const_pointer data() const { return data_; }
    
    reference front() { return data_[0]; }
    reference back()  { return data_[size_-1]; }
    const_reference front() const { return data_[0]; }
    const_reference back()  const { return data_[size_-1]; }
    
    bool       empty() const { return size_ == 0; }
    ArrayIndex size()  const { return size_; }
    ArrayIndex max_size()  const { return size_; }
    ArrayIndex ndim()  const { return static_ndim; }

    TinyArrayBase & reverse()
    {
        ArrayIndex i=0, j=size_-1;
        while(i < j)
             std::swap(data_[i++], data_[j--]);
        return *this;
    }

    void swap(TinyArrayBase & other)
    {
        std::swap(size_, other.size_);
        std::swap(data_, other.data_);
    }

        /// factory function for the fixed-size k-th unit vector 
    static inline 
    TinyArray<value_type, 0>
    unitVector(ArrayIndex size, ArrayIndex k)
    {
        TinyArray<value_type, 0> res(size);
        res[k] = 1;
        return res;
    }

        // /// factory function for fixed-size linear sequence starting at <tt>start</tt> with stepsize <tt>step</tt>
    // static inline 
    // TinyArray<value_type, N...>
    // linearSequence(value_type start = value_type(), value_type step = value_type(1))
    // {
        // TinyArray<value_type, N...> res;
        // for(int k=0; k < static_size; ++k, start += step)
            // res[k] = start;
        // return res;
    // }
    
  protected:
    ArrayIndex size_;
    pointer data_;
};

template <class VALUETYPE>
class TinyArray<VALUETYPE, 0>
: public TinyArrayBase<VALUETYPE, TinyArray<VALUETYPE, 0>, 0>
{
    using BaseType = TinyArrayBase<VALUETYPE, TinyArray<VALUETYPE, 0>, 0>;
  public:
  
    using value_type = VALUETYPE;
    
    TinyArray()
    : BaseType()
    {}

    TinyArray(ArrayIndex size, SkipInitialization)
    : BaseType(size)
    {
        this->data_ = alloc_.allocate(this->size_);
    }

    explicit TinyArray(ArrayIndex size)
    : TinyArray(size, DontInit)
    {
        std::uninitialized_fill(this->begin(), this->end(), value_type());
    }

    TinyArray( ArrayIndex size, value_type const & initial)
    : TinyArray(size)
    {
        std::uninitialized_fill(this->begin(), this->end(), initial);
    }

    TinyArray(ArrayIndex size, lemon::Invalid const &)
    : TinyArray(size, value_type(-1))
    {}
    
    TinyArray( TinyArray const & rhs )
    : TinyArray(rhs.size(), DontInit)
    {
        std::uninitialized_copy(rhs.begin(), rhs.end(), this->begin());
    }

    template <class U, class D, int N>
    explicit TinyArray(TinyArrayBase<U, D, N> const & rhs)
    : TinyArray(rhs.size(), DontInit)
    {
        std::uninitialized_copy(rhs.begin(), rhs.end(), this->begin());
    }

    template <class U>
    TinyArray(std::initializer_list<U> rhs)
    : TinyArray(rhs.size(), DontInit)
    {
        std::uninitialized_copy(rhs.begin(), rhs.end(), this->begin());
    }
    
    template <class InputIterator>
    TinyArray(InputIterator i, InputIterator end)
    : TinyArray(std::distance(i, end), DontInit)
    {
        std::uninitialized_copy(i, end, this->begin());
    }

    TinyArray & operator=(value_type v)
    {
        BaseType::operator=(v);
        return *this;
    }

    TinyArray & operator=( TinyArray const & rhs )
    {
        if(this == &rhs)
            return *this;
        if(this->size_ == 0)
            TinyArray(rhs).swap(*this);
        else
            BaseType::operator=(rhs);
        return *this;
    }

    template <class U, class D, int N>
    TinyArray & operator=(TinyArrayBase<U, D, N> const & rhs)
    {
        if(this->size_ == 0)
            TinyArray(rhs).swap(*this);
        else
            BaseType::operator=(rhs);
        return *this;
    }

    ~TinyArray()
    {
        if(this->data_)
        {
            for(ArrayIndex i=0; i<this->size_; ++i)
                (this->data_+i)->~value_type();
            alloc_.deallocate(this->data_, this->size_);
        }
    }

  private:
    // FIXME: implement an optimized allocator
    std::allocator<value_type> alloc_;
};

} // namespace vigra

#endif // VIGRA_TINYARRAY_DYNAMIC_HXX
