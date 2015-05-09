/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
/*                                                                      */
/*    This file is part of the MULI computer vision library.           */
/*    The MULI Website is                                              */
/*        http://hci.iwr.uni-heidelberg.de/muli/                       */
/*    Please direct questions, bug reports, and contributions to        */
/*        ullrich.koethe@iwr.uni-heidelberg.de    or                    */
/*        muli@informatik.uni-hamburg.de                               */
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


#ifndef MULI_TINYARRAY_HXX
#define MULI_TINYARRAY_HXX

namespace lemon {

struct Invalid;

} // namespace lemon

#include <cmath>    // abs(double)
#include <cstdlib>  // abs(int)
#include <iosfwd>   // ostream
#include <algorithm>
#include <iterator>
#include <array>
#include "config.hxx"
#include "promote.hxx"
// #include "error.hxx"
// #include "mathutil.hxx"

#ifdef MULI_CHECK_BOUNDS
#define MULI_ASSERT_INSIDE(diff) \
  muli_precondition(diff >= 0, "Index out of bounds");\
  muli_precondition(diff < SIZE, "Index out of bounds");
#else
#define MULI_ASSERT_INSIDE(diff)
#endif

namespace muli {

// mask cl.exe shortcomings [begin]
#if defined(_MSC_VER)
#pragma warning( push )
#pragma warning( disable : 4503 )
#endif

using std::abs;
using std::ceil;
using std::floor;
using std::sqrt;

using ArrayIndex = std::ptrdiff_t;

enum SkipInitialization { DontInit };

template <ArrayIndex LEVEL, ArrayIndex ... N>
struct TinyShapeImpl;

template <ArrayIndex LEVEL, ArrayIndex N, ArrayIndex ... REST>
struct TinyShapeImpl<LEVEL, N, REST...>
{
    using NextType = TinyShapeImpl<LEVEL+1, REST...>;
    
    static const ArrayIndex level      = LEVEL;
    static const ArrayIndex stride     = NextType::total_size;
    static const ArrayIndex total_size = N * stride;
    
    static ArrayIndex offset(ArrayIndex const * coord)
    {
        return stride*coord[level] + NextType::offset(coord);
    }
    
    template <class ... V>
    static ArrayIndex offset(ArrayIndex i, V...rest)
    {
        return stride*i + NextType::offset(rest...);
    }
};

template <ArrayIndex LEVEL, ArrayIndex N>
struct TinyShapeImpl<LEVEL, N>
{
    static const ArrayIndex level      = LEVEL;
    static const ArrayIndex stride     = 1;
    static const ArrayIndex total_size = N;
    
    static ArrayIndex offset(ArrayIndex const * coord)
    {
        return coord[level];
    }
    
    static ArrayIndex offset(ArrayIndex i)
    {
        return i;
    }
};

template <class VALUETYPE, ArrayIndex ... N>
class TinyArray;

template <class VALUETYPE, ArrayIndex ... N>
class TinyArrayView;

template <class VALUETYPE, class DERIVED, ArrayIndex ... N>
class TinyArrayBase
{
  protected:
    using ShapeHelper = TinyShapeImpl<0, N...>;
    
    static const bool derived_is_view = std::is_same<DERIVED, TinyArrayView<VALUETYPE, N...> >::value;
    using data_array_type = typename std::conditional<derived_is_view, 
                                                VALUETYPE *, 
                                                VALUETYPE[ShapeHelper::total_size]>::type;
    
  public:
  
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
    using index_type             = TinyArray<ArrayIndex, sizeof...(N)>;
    
    static constexpr ArrayIndex static_ndim  = sizeof...(N);
    static constexpr ArrayIndex static_size  = ShapeHelper::total_size;
    static constexpr index_type static_shape = index_type(N...);
    
    // constructors
        
    constexpr TinyArrayBase(TinyArrayBase const &) = default;
    
    TinyArrayBase(SkipInitialization)
    {}
    
    TinyArrayBase(value_type v = value_type())
    {
        for(ArrayIndex i=0; i<static_size; ++i)
            data_[i] = v;
    }
    
    constexpr TinyArrayBase(value_type const (&v)[static_ndim])
    : data_{v}
    {}
    
    // template <class ... V>
    // constexpr TinyArrayBase(V... v)
    // : data_{v...}
    // {
        // static_assert(sizeof...(V) == static_size, 
                      // "TinyArrayBase(): wrong number of arguments.");
    // }
    
    template <class OTHER, class OTHER_DERIVED>
    TinyArrayBase(TinyArrayBase<OTHER, OTHER_DERIVED, N...> const & other)
    {
        for(ArrayIndex i=0; i<static_size; ++i)
            data_[i] = other[i];
    }
    
    // assignment
    
    TinyArrayBase & operator=(TinyArrayBase const &) = default;

    TinyArrayBase & operator=(value_type v)
    {
        for(ArrayIndex i=0; i<static_size; ++i)
            data_[i] = v;
        return *this;
    }
    
    TinyArrayBase & operator=(value_type const (&v)[static_ndim])
    {
        for(ArrayIndex i=0; i<static_size; ++i)
            data_[i] = v[i];
        return *this;
    }
    
    template <class OTHER, class OTHER_DERIVED>
    TinyArrayBase & operator=(TinyArrayBase<OTHER, OTHER_DERIVED, N...> const & other)
    {
        for(ArrayIndex i=0; i<static_size; ++i)
            data_[i] = other[i];
        return *this;
    }
    
    // arithmetic assignment
    
    DERIVED & operator+=(value_type v)
    {
        for(ArrayIndex i=0; i<static_size; ++i)
            data_[i] += v;
        return static_cast<DERIVED &>(*this);
    }
    
    template <class OTHER, class OTHER_DERIVED>
    DERIVED & operator+=(TinyArrayBase<OTHER, OTHER_DERIVED, N...> const & other)
    {
        for(ArrayIndex i=0; i<static_size; ++i)
            data_[i] += other[i];
        return static_cast<DERIVED &>(*this);
    }
    
    DERIVED & operator-=(value_type v)
    {
        for(ArrayIndex i=0; i<static_size; ++i)
            data_[i] -= v;
        return static_cast<DERIVED &>(*this);
    }
    
    template <class OTHER, class OTHER_DERIVED>
    DERIVED & operator-=(TinyArrayBase<OTHER, OTHER_DERIVED, N...> const & other)
    {
        for(ArrayIndex i=0; i<static_size; ++i)
            data_[i] -= other[i];
        return static_cast<DERIVED &>(*this);
    }
    
    DERIVED & operator*=(value_type v)
    {
        for(ArrayIndex i=0; i<static_size; ++i)
            data_[i] *= v;
        return static_cast<DERIVED &>(*this);
    }
    
    template <class OTHER, class OTHER_DERIVED>
    DERIVED & operator*=(TinyArrayBase<OTHER, OTHER_DERIVED, N...> const & other)
    {
        for(ArrayIndex i=0; i<static_size; ++i)
            data_[i] *= other[i];
        return static_cast<DERIVED &>(*this);
    }
    
    DERIVED & operator/=(value_type v)
    {
        for(ArrayIndex i=0; i<static_size; ++i)
            data_[i] /= v;
        return static_cast<DERIVED &>(*this);
    }
    
    template <class OTHER, class OTHER_DERIVED>
    DERIVED & operator/=(TinyArrayBase<OTHER, OTHER_DERIVED, N...> const & other)
    {
        for(ArrayIndex i=0; i<static_size; ++i)
            data_[i] /= other[i];
        return static_cast<DERIVED &>(*this);
    }
    
    DERIVED & operator%=(value_type v)
    {
        for(ArrayIndex i=0; i<static_size; ++i)
            data_[i] %= v;
        return static_cast<DERIVED &>(*this);
    }
    
    template <class OTHER, class OTHER_DERIVED>
    DERIVED & operator%=(TinyArrayBase<OTHER, OTHER_DERIVED, N...> const & other)
    {
        for(ArrayIndex i=0; i<static_size; ++i)
            data_[i] %= other[i];
        return static_cast<DERIVED &>(*this);
    }
    
    // vector functions
    
    Promote<value_type> squaredNorm() const
    {
        Promote<value_type> result = Promote<value_type>();
        for(ArrayIndex i=0; i<static_size; ++i)
            result += data_[i]*data_[i];
        return result;
    }
    
    double norm() const
    {
         return sqrt(this->squaredNorm());
    }


        /** Return the minimal element.
        */
    const_reference minimum() const
    {
        ArrayIndex m = 0;
        for(ArrayIndex i=1; i<static_size; ++i)
            if(data_[i] < data_[m])
                m = i;
        return data_[m];
    }

        /** Return the maximal element.
        */
    const_reference maximum() const
    {
        ArrayIndex m = 0;
        for(ArrayIndex i=1; i<static_size; ++i)
            if(data_[m] < data_[i])
                m = i;
        return data_[m];
    }
    
        /** Check that all elements of this vector are non-zero (or 'true' if T is bool).
        */
    bool all() const
    {
        for(ArrayIndex i=0; i<static_size; ++i)
            if(data_[i] == value_type())
                return false;
        return true;
    }
    
        /** Check that at least one element of this vector is non-zero (or 'true' if T is bool).
        */
    bool any() const
    {
        for(ArrayIndex i=0; i<static_size; ++i)
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
        if(i < 0 || i >= static_size)
            throw std::out_of_range("TinyArrayBase::at()");
        return data_[i];
    }
    
    const_reference at(ArrayIndex i) const
    {
        if(i < 0 || i >= static_size)
            throw std::out_of_range("TinyArrayBase::at()");
        return data_[i];
    }
    
    reference operator[](ArrayIndex const (&i)[static_ndim])
    {
        return data_[ShapeHelper::offset(i)];
    }
    
    constexpr const_reference operator[](ArrayIndex const (&i)[static_ndim]) const
    {
        return data_[ShapeHelper::offset(i)];
    }
    
    reference at(ArrayIndex const (&i)[static_ndim])
    {
        return at(ShapeHelper::offset(i));
    }
    
    const_reference at(ArrayIndex const (&i)[static_ndim]) const
    {
        return at(ShapeHelper::offset(i));
    }
    
    reference operator[](index_type const & i)
    {
        return data_[ShapeHelper::offset(i.data())];
    }
    
    constexpr const_reference operator[](index_type const & i) const
    {
        return data_[ShapeHelper::offset(i.data())];
    }
    
    reference at(index_type const & i)
    {
        return at(ShapeHelper::offset(i.data()));
    }
    
    const_reference at(index_type const & i) const
    {
        return at(ShapeHelper::offset(i.data()));
    }
    
    template <class ... V>
    reference operator()(V...v)
    {
        static_assert(sizeof...(V) == static_ndim, 
                      "TinyArrayBase::operator(): wrong number of arguments.");
        return data_[ShapeHelper::offset(v...)];
    }
    
    template <class ... V>
    constexpr const_reference operator()(V...v) const
    {
        static_assert(sizeof...(V) == static_ndim, 
                      "TinyArrayBase::operator(): wrong number of arguments.");
        return data_[ShapeHelper::offset(v...)];
    }
    
    // boiler plate
    
    iterator begin() { return data_; }
    iterator end()   { return data_ + static_size; }
    const_iterator begin() const { return data_; }
    const_iterator end()   const { return data_ + static_size; }
    const_iterator cbegin() const { return data_; }
    const_iterator cend()   const { return data_ + static_size; }
    
    reverse_iterator rbegin() { return data_ + static_size; }
    reverse_iterator rend()   { return data_; }
    const_reverse_iterator rbegin() const { return data_ + static_size; }
    const_reverse_iterator rend()   const { return data_; }
    const_reverse_iterator crbegin() const { return data_ + static_size; }
    const_reverse_iterator crend()   const { return data_; }
    
    pointer data() { return data_; }
    const_pointer data() const { return data_; }
    
    reference front() { return data_[0]; }
    reference back()  { return data_[static_size-1]; }
    constexpr const_reference front() const { return data_[0]; }
    constexpr const_reference back()  const { return data_[static_size-1]; }
    
    constexpr index_type shape() const { return static_shape; }
    constexpr ArrayIndex size()  const { return static_size; }
    constexpr ArrayIndex max_size()  const { return static_size; }
    constexpr ArrayIndex ndim()  const { return static_ndim; }
    constexpr bool empty() const { return static_size == 0; }
    
  protected:
    data_array_type data_;
};

template <class T, class DERIVED, ArrayIndex ... N>
constexpr 
typename TinyArrayBase<T, DERIVED, N...>::index_type 
TinyArrayBase<T, DERIVED, N...>::static_shape;

template <class T, class DERIVED, ArrayIndex ... N>
std::ostream & operator<<(std::ostream & o, TinyArrayBase<T, DERIVED, N...> const & v)
{
    if(TinyArrayBase<T, DERIVED, N...>::static_size > 0)
        o << v[0];
    for(ArrayIndex i=1; i < TinyArrayBase<T, DERIVED, N...>::static_size; ++i)
        o << " " << v[i];
    return o;
}

template <class VALUETYPE, ArrayIndex ... N>
class TinyArray
: public TinyArrayBase<VALUETYPE, TinyArray<VALUETYPE, N...>, N...>
{
    using BaseType = TinyArrayBase<VALUETYPE, TinyArray<VALUETYPE, N...>, N...>;
    
    template <int LEVEL, class V1, class ... V2>
    void initImpl(V1 v1, V2... v2)
    {
        BaseType::data_[LEVEL] =  v1;
        initImpl<LEVEL+1>(v2...);
    }
    
    template <int LEVEL, class V1>
    void initImpl(V1 v1)
    {
        BaseType::data_[LEVEL] =  v1;
    }
    
  public:
  
    template <class NEW_VALUETYPE>
    using AsType = TinyArray<NEW_VALUETYPE, N...>;
    
    constexpr TinyArray(TinyArray const &) = default;
    
    TinyArray(VALUETYPE v = VALUETYPE())
    : BaseType(v)
    {}
    
    constexpr TinyArray(VALUETYPE const (&v)[BaseType::static_ndim])
    : BaseType(v)
    {}
    
    template <class ... V>
    // constexpr TinyArray(V... v)
    TinyArray(V... v)
    : BaseType(DontInit)
    {
        initImpl<0>(v...);
    }
    
    template <class OTHER>
    TinyArray(TinyArray<OTHER, N...> const & other)
    : BaseType(static_cast<TinyArrayBase<OTHER, TinyArray<OTHER, N...>, N...> const &>(other))
    {}
    
    template <class OTHER, class DERIVED>
    TinyArray(TinyArrayBase<OTHER, DERIVED, N...> const & other)
    : BaseType(other)
    {}
    
    template <class ... V>
    void init(V... v)
    {
        static_assert(sizeof...(V) == BaseType::static_size, 
                      "TinyArray::init(): wrong number of arguments.");
        initImpl<0>(v...);
    }
    
    TinyArray & operator=(TinyArray const &) = default;

    TinyArray & operator=(VALUETYPE v)
    {
        BaseType::operator=(v);
        return *this;
    }
    
    TinyArray & operator=(VALUETYPE const (&v)[BaseType::static_ndim])
    {
        BaseType::operator=(v);
        return *this;
    }
    
    template <class OTHER, class OTHER_DERIVED>
    TinyArray & operator=(TinyArrayBase<OTHER, OTHER_DERIVED, N...> const & other)
    {
        BaseType::operator=(other);
        return *this;
    }
};


    /// component-wise addition
template <class V1, class D1, class V2, class D2, ArrayIndex ... N>
inline
TinyArray<Promote<V1, V2>, N...>
operator+(TinyArrayBase<V1, D1, N...> const & l,
          TinyArrayBase<V2, D2, N...> const & r)
{
    return TinyArray<Promote<V1, V2>, N...>(l) += r;
}

#if 0

/********************************************************/
/*                                                      */
/*                    TinyArrayBase                    */
/*                                                      */
/********************************************************/

/** \brief Base class for fixed size vectors.

    This class contains functionality shared by
    \ref TinyVector and \ref TinyVectorView, and enables these classes
    to be freely mixed within expressions. It is typically not used directly.

    <b>\#include</b> \<muli/tinyvector.hxx\><br>
    Namespace: muli
**/
template <class VALUETYPE, int SIZE, class DATA, class DERIVED>
class TinyArrayBase
{
    TinyArrayBase(TinyArrayBase const &); // do not use

    TinyArrayBase & operator=(TinyArrayBase const & other); // do not use

  protected:

    typedef typename detail::LoopType<SIZE>::type Loop;

    TinyArrayBase()
    {}

  public:
        /** STL-compatible definition of valuetype
        */
    typedef VALUETYPE value_type;

        /** reference (return of operator[]).
        */
    typedef VALUETYPE & reference;

        /** const reference (return of operator[] const).
        */
    typedef VALUETYPE const & const_reference;

        /** pointer (return of operator->).
        */
    typedef VALUETYPE * pointer;

        /** const pointer (return of operator-> const).
        */
    typedef VALUETYPE const * const_pointer;

        /** STL-compatible definition of iterator
        */
    typedef value_type * iterator;

        /** STL-compatible definition of const iterator
        */
    typedef value_type const * const_iterator;

        /** STL-compatible definition of size_type
        */
    typedef unsigned int size_type;

        /** STL-compatible definition of difference_type
        */
    typedef std::ptrdiff_t difference_type;

        /** the scalar type for the outer product
        */
    typedef double scalar_multiplier;

        /** the vector's squared norm type
        */
    typedef typename NormTraits<VALUETYPE>::SquaredNormType SquaredNormType;

        /** the vector's norm type
        */
    typedef typename SquareRootTraits<SquaredNormType>::SquareRootResult NormType;

        /** the vector's size
        */
    enum { static_size = SIZE };

        /** Initialize from another sequence (must have length SIZE!)
        */
    template <class Iterator>
    void init(Iterator i, Iterator end)
    {
        muli_precondition(end-i == SIZE,
            "TinyVector::init(): Sequence has wrong size.");
        Loop::assignCast(data_, i);
    }

        /** Initialize with constant value
        */
    void init(value_type initial)
    {
        Loop::assignScalar(data_, initial);
    }

        /** Component-wise add-assignment
        */
    template <class T1, class D1, class D2>
    DERIVED & operator+=(TinyArrayBase<T1, SIZE, D1, D2> const & r)
    {
        Loop::add(data_, r.begin());
        return static_cast<DERIVED &>(*this);
    }

        /** Component-wise subtract-assignment
        */
    template <class T1, class D1, class D2>
    DERIVED & operator-=(TinyArrayBase<T1, SIZE, D1, D2> const & r)
    {
        Loop::sub(data_, r.begin());
        return static_cast<DERIVED &>(*this);
    }

        /** Component-wise multiply-assignment
        */
    template <class T1, class D1, class D2>
    DERIVED & operator*=(TinyArrayBase<T1, SIZE, D1, D2> const & r)
    {
        Loop::mul(data_, r.begin());
        return static_cast<DERIVED &>(*this);
    }

        /** Component-wise divide-assignment
        */
    template <class T1, class D1, class D2>
    DERIVED & operator/=(TinyArrayBase<T1, SIZE, D1, D2> const & r)
    {
        Loop::div(data_, r.begin());
        return static_cast<DERIVED &>(*this);
    }

        /** Component-wise modulo-assignment
        */
    template <class T1, class D1, class D2>
    DERIVED & operator%=(TinyArrayBase<T1, SIZE, D1, D2> const & r)
    {
        Loop::mod(data_, r.begin());
        return static_cast<DERIVED &>(*this);
    }

        /** Component-wise scalar multiply-assignment
        */
    DERIVED & operator+=(double r)
    {
        Loop::addScalar(data_, r);
        return static_cast<DERIVED &>(*this);
    }

        /** Component-wise scalar divide-assignment
        */
    DERIVED & operator-=(double r)
    {
        Loop::subScalar(data_, r);
        return static_cast<DERIVED &>(*this);
    }

        /** Component-wise scalar multiply-assignment
        */
    DERIVED & operator*=(double r)
    {
        Loop::mulScalar(data_, r);
        return static_cast<DERIVED &>(*this);
    }

        /** Component-wise scalar divide-assignment
        */
    DERIVED & operator/=(double r)
    {
        Loop::divScalar(data_, r);
        return static_cast<DERIVED &>(*this);
    }

        /** Calculate magnitude (i.e. 2-norm / Euclidean norm / length).
         * \see squaredMagnitude()
         */
    NormType magnitude() const
    {
         return sqrt(static_cast<typename
              SquareRootTraits<SquaredNormType>::SquareRootArgument>(squaredMagnitude()));
    }

        /** Calculate squared magnitude (i.e. sum of squared elements).
        */
    SquaredNormType squaredMagnitude() const
    {
        return Loop::squaredNorm(data_);
    }

        /** Return the minimal element.
        */
    VALUETYPE const & minimum() const
    {
        return Loop::minimum(data_);
    }

        /** Return the maximal element.
        */
    VALUETYPE const & maximum() const
    {
        return Loop::maximum(data_);
    }
    
        /** Check that all elements of this vector are non-zero (or 'true' if T is bool).
        */
    bool all() const
    {
        return Loop::all(data_, VALUETYPE());
    }
    
        /** Check that at least one element of this vector is non-zero (or 'true' if T is bool).
        */
    bool any() const
    {
        return Loop::any(data_, VALUETYPE());
    }

        /** Access component by index.
        */
    reference operator[](difference_type i)
    {
        MULI_ASSERT_INSIDE(i);
        return data_[i];
    }

        /** Get component by index.
        */
    const_reference operator[](difference_type i) const
    {
        MULI_ASSERT_INSIDE(i);
        return data_[i];
    }

        /** Get random access iterator to begin of vector.
        */
    iterator begin() { return data_; }
        /** Get random access iterator past-the-end of vector.
        */
    iterator end() { return data_ + SIZE; }

        /** Get const random access iterator to begin of vector.
        */
    const_iterator begin() const { return data_; }

        /** Get const random access iterator past-the-end of vector.
        */
    const_iterator end() const { return data_ + SIZE; }

        /** Get const random access iterator to begin of vector.
        */
    const_iterator cbegin() const { return data_; }

        /** Get const random access iterator past-the-end of vector.
        */
    const_iterator cend() const { return data_ + SIZE; }
    
        /** Get a view to the subarray with length <tt>(TO-FROM)</tt> starting at <tt>FROM</tt>.
            The bounds must fullfill <tt>0 <= FROM < TO <= SIZE</tt>, but this is only
            checked when <tt>MULI_CHECK_BOUNDS</tt> is \#define'd.
        */
    template <int FROM, int TO>
    TinyVectorView<VALUETYPE, TO-FROM> subarray() const
    {
#ifdef MULI_CHECK_BOUNDS
        muli_precondition(FROM >= 0, "Index out of bounds");
        muli_precondition(FROM < TO, "Index out of bounds");
        muli_precondition(TO <=SIZE, "Index out of bounds");
#endif
        return TinyVectorView<VALUETYPE, TO-FROM>(data_+FROM);
    }
    
    TinyVector<VALUETYPE, SIZE-1>
    dropIndex(int m) const
    {
#ifdef MULI_CHECK_BOUNDS
        muli_precondition(0 <= m && m < SIZE, "Dimension out of bounds");
#endif
        TinyVector<VALUETYPE, SIZE-1> res(SkipInitialization);
        for(int k=0; k<m; ++k)
            res[k] = data_[k];
        for(int k=m; k<SIZE-1; ++k)
            res[k] = data_[k+1];
        return res;
    }

        /** Size of TinyVector vector always equals the template parameter SIZE.
        */
    size_type size() const { return SIZE; }

    pointer data() { return data_; }

    const_pointer data() const { return data_; }
    
    reference front()
    {
        return data_[0];
    }
    
    const_reference front() const
    {
        return data_[0];
    }
    
    reference back()
    {
        return data_[SIZE-1];
    }
    
    const_reference back() const
    {
        return data_[SIZE-1];
    }
    
        /** \brief Factory function for a unit vector for dimension \a k.
        */
    static TinyVector<VALUETYPE, SIZE> unitVector(int k)
    {
        MULI_ASSERT_INSIDE(k);
        TinyVector<VALUETYPE, SIZE> ret;
        ret[k] = 1;
        return ret;
    }
    
        /** \brief Factory function for a linear sequence.
        
            The result will be initialized as <tt>res[k] = start + k*step</tt>.
        */
    static TinyVector<VALUETYPE, SIZE> linearSequence(VALUETYPE start=VALUETYPE(), VALUETYPE step=VALUETYPE(1))
    {
        TinyVector<VALUETYPE, SIZE> ret(SkipInitialization);
        for(int k=0; k<SIZE; ++k, start+=step)
            ret[k] = start;
        return ret;
    }

  protected:
  
    DATA data_;
};

#ifndef DOXYGEN

template <int SIZE, int DESIRED_SIZE>
struct TinyVector_constructor_has_wrong_number_of_arguments
: staticAssert::AssertBool<SIZE == DESIRED_SIZE>
{};

#endif /* DOXYGEN */

/** \brief Class for fixed size vectors.
    \ingroup RangesAndPoints

    This class contains an array of size SIZE of the specified VALUETYPE.
    The interface conforms to STL vector, except that there are no functions
    that change the size of a TinyVector.

    \ref TinyVectorOperators "Arithmetic operations"
    on TinyVectors are defined as component-wise applications of these
    operations. Addition and subtraction of two TinyVectors
    (+=, -=, +, -, unary -), multiplication and division of an
    TinyVector with a double, and NumericTraits/PromoteTraits are defined,
    so that TinyVector fulfills the requirements of \ref LinearAlgebraConcept "Linear Algebra".

    MULI algorithms typically use \ref muli::VectorAccessor to access
    TinyVectors as a whole, or specific components of them.

    See also:<br>
    <UL style="list-style-image:url(documents/bullet.gif)">
        <LI> \ref muli::TinyArrayBase
        <LI> \ref muli::TinyVectorView
        <LI> \ref TinyVectorTraits
        <LI> \ref TinyVectorOperators
    </UL>

    <b>\#include</b> \<muli/tinyvector.hxx\><br>
    Namespace: muli
**/
template <class T, int SIZE>
class TinyVector
: public TinyArrayBase<T, SIZE, T[SIZE], TinyVector<T, SIZE> >
{
    typedef TinyArrayBase<T, SIZE, T[SIZE], TinyVector<T, SIZE> > BaseType;
    typedef typename BaseType::Loop Loop;

  public:

    typedef typename BaseType::value_type value_type;
    typedef typename BaseType::reference reference;
    typedef typename BaseType::const_reference const_reference;
    typedef typename BaseType::pointer pointer;
    typedef typename BaseType::const_pointer const_pointer;
    typedef typename BaseType::iterator iterator;
    typedef typename BaseType::const_iterator const_iterator;
    typedef typename BaseType::size_type size_type;
    typedef typename BaseType::difference_type difference_type;
    typedef typename BaseType::scalar_multiplier scalar_multiplier;
    typedef typename BaseType::SquaredNormType SquaredNormType;
    typedef typename BaseType::NormType NormType;
    
    enum ReverseCopyTag { ReverseCopy };

        /** Construction with constant value.
        
            Initializes all vector elements with the given value.
        */
    explicit TinyVector(value_type const & initial)
    : BaseType()
    {
        Loop::assignScalar(BaseType::begin(), initial);
    }

        /** Construction from lemon::Invalid.
        
            Initializes all vector elements with -1.
        */
    TinyVector(lemon::Invalid const &)
    : BaseType()
    {
        Loop::assignScalar(BaseType::begin(), -1);
    }

        /** Construction with Diff2D.
        
            Use only when <tt>SIZE == 2</tt>.
        */
    explicit TinyVector(Diff2D const & initial)
    : BaseType()
    {
        BaseType::data_[0] = detail::RequiresExplicitCast<T>::cast(initial.x);
        BaseType::data_[1] = detail::RequiresExplicitCast<T>::cast(initial.y);
    }

        /** Construction with explicit values.
            Call only if SIZE == 2
        */
    TinyVector(value_type const & i1, value_type const & i2)
    : BaseType()
    {
        MULI_STATIC_ASSERT((TinyVector_constructor_has_wrong_number_of_arguments<SIZE, 2>));
        BaseType::data_[0] = i1;
        BaseType::data_[1] = i2;
    }

        /** Construction with explicit values.
            Call only if SIZE == 3
        */
    TinyVector(value_type const & i1, value_type const & i2, value_type const & i3)
    : BaseType()
    {
        MULI_STATIC_ASSERT((TinyVector_constructor_has_wrong_number_of_arguments<SIZE, 3>));
        BaseType::data_[0] = i1;
        BaseType::data_[1] = i2;
        BaseType::data_[2] = i3;
    }

        /** Construction with explicit values.
            Call only if SIZE == 4
        */
    TinyVector(value_type const & i1, value_type const & i2,
               value_type const & i3, value_type const & i4)
    : BaseType()
    {
        MULI_STATIC_ASSERT((TinyVector_constructor_has_wrong_number_of_arguments<SIZE, 4>));
        BaseType::data_[0] = i1;
        BaseType::data_[1] = i2;
        BaseType::data_[2] = i3;
        BaseType::data_[3] = i4;
    }

        /** Construction with explicit values.
            Call only if SIZE == 5
        */
    TinyVector(value_type const & i1, value_type const & i2,
               value_type const & i3, value_type const & i4,
               value_type const & i5)
    : BaseType()
    {
        MULI_STATIC_ASSERT((TinyVector_constructor_has_wrong_number_of_arguments<SIZE, 5>));
        BaseType::data_[0] = i1;
        BaseType::data_[1] = i2;
        BaseType::data_[2] = i3;
        BaseType::data_[3] = i4;
        BaseType::data_[4] = i5;
    }
    
       /** Default constructor (initializes all elements with zero).
        */
    TinyVector()
    : BaseType()
    {
        Loop::assignScalar(BaseType::data_, value_type());
    }

        /** Construct without initializing the vector elements.
        */
    explicit TinyVector(SkipInitializationTag)
    : BaseType()
    {}

    explicit TinyVector(detail::DontInit)
    : BaseType()
    {}

        /** Copy constructor.
        */
    TinyVector(TinyVector const & r)
    : BaseType()
    {
        Loop::assign(BaseType::data_, r.data_);
    }

        /** Constructor from C array.
        */
    template <class U>
    explicit TinyVector(U const * data)
    : BaseType()
    {
        Loop::assign(BaseType::data_, data);
    }

        /** Constructor by reverse copy from C array.
            
            Usage:
            \code
            TinyVector<int, 3> v(1,2,3);
            TinyVector<int, 3> reversed(v.begin(), TinyVector<int, 3>::ReverseCopy);
            \endcode
        */
    explicit TinyVector(const_pointer data, ReverseCopyTag)
    : BaseType()
    {
        Loop::reverseAssign(BaseType::data_, data+SIZE-1);
    }

        /** Copy with type conversion.
        */
    template <class U, class DATA, class DERIVED>
    TinyVector(TinyArrayBase<U, SIZE, DATA, DERIVED> const & r)
    : BaseType()
    {
        Loop::assignCast(BaseType::data_, r.begin());
    }

        /** Copy assignment.
        */
    TinyVector & operator=(TinyVector const & r)
    {
        Loop::assign(BaseType::data_, r.data_);
        return *this;
    }

        /** Copy assignment with type conversion.
        */
    template <class U, class DATA, class DERIVED>
    TinyVector & operator=(TinyArrayBase<U, SIZE, DATA, DERIVED> const & r)
    {
        Loop::assignCast(BaseType::data_, r.begin());
        return *this;
    }

        /** Assignment from Diff2D.
        
            Use only when <tt>SIZE == 2</tt>.
        */
    TinyVector & operator=(Diff2D const & r)
    {
        BaseType::data_[0] = detail::RequiresExplicitCast<T>::cast(r.x);
        BaseType::data_[1] = detail::RequiresExplicitCast<T>::cast(r.y);
        return *this;
    }

        /** Assignment from scalar. Will set all entries to the given value.
        */
    TinyVector & operator=(value_type const & v)
    {
        Loop::assignScalar(BaseType::begin(), v);
        return *this;
    }

        /** Copy from a TinyVector with a different number of elements.
        
            Only the first <tt>min(SIZE, USIZE)</tt> elements are copied.
        */
    template <class U, int USIZE, class DATA, class DERIVED>
    TinyVector & copy(TinyArrayBase<U, USIZE, DATA, DERIVED> const & r)
    {
        static const int minSize = USIZE < SIZE
                                        ? USIZE
                                        : SIZE;
        
        typedef typename detail::LoopType<minSize>::type MinLoop;
        MinLoop::assignCast(BaseType::data_, r.begin());
        return *this;
    }
};

/** \brief Wrapper for fixed size vectors.

    This class wraps an array of size SIZE of the specified VALUETYPE.
    Thus, the array can be accessed with an interface similar to
    that of std::vector (except that there are no functions
    that change the size of a TinyVectorView). The TinyVectorView
    does <em>not</em> assume ownership of the given memory.

    \ref TinyVectorOperators "Arithmetic operations"
    on TinyVectorViews are defined as component-wise applications of these
    operations. Addition and subtraction of two TinyVectorViews
    (+=, -=, +, -, unary -), multiplication and division of an
    TinyVectorViews with a double, and NumericTraits/PromoteTraits are defined,
    so that TinyVectorView fulfills the requirements of \ref LinearAlgebraConcept "Linear Algebra".

    MULI algorithms typically use \ref muli::VectorAccessor to access
    TinyVectorViews as a whole, or specific components of them.

    <b>See also:</b>
    <ul>
        <li> \ref muli::TinyArrayBase
        <li> \ref muli::TinyVector
        <li> \ref TinyVectorTraits
        <li> \ref TinyVectorOperators
    </ul>

    <b>\#include</b> \<muli/tinyvector.hxx\><br>
    Namespace: muli
**/
template <class T, int SIZE>
class TinyVectorView
: public TinyArrayBase<T, SIZE, T *, TinyVectorView<T, SIZE> >
{
    typedef TinyArrayBase<T, SIZE, T *, TinyVectorView<T, SIZE> > BaseType;
    typedef typename BaseType::Loop Loop;

  public:

    typedef typename BaseType::value_type value_type;
    typedef typename BaseType::reference reference;
    typedef typename BaseType::const_reference const_reference;
    typedef typename BaseType::pointer pointer;
    typedef typename BaseType::const_pointer const_pointer;
    typedef typename BaseType::iterator iterator;
    typedef typename BaseType::const_iterator const_iterator;
    typedef typename BaseType::size_type size_type;
    typedef typename BaseType::difference_type difference_type;
    typedef typename BaseType::scalar_multiplier scalar_multiplier;
    typedef typename BaseType::SquaredNormType SquaredNormType;
    typedef typename BaseType::NormType NormType;

        /** Default constructor
            (pointer to wrapped data is NULL).
        */
    TinyVectorView()
    : BaseType()
    {
        BaseType::data_ = 0;
    }

        /** Construct view for given data array
        */
    TinyVectorView(const_pointer data)
    : BaseType()
    {
        BaseType::data_ = const_cast<pointer>(data);
    }

        /** Copy constructor (shallow copy).
        */
    TinyVectorView(TinyVectorView const & other)
    : BaseType()
    {
        BaseType::data_ = const_cast<pointer>(other.data_);
    }

        /** Construct view from other TinyVector.
        */
    template <class DATA, class DERIVED>
    TinyVectorView(TinyArrayBase<T, SIZE, DATA, DERIVED> const & other)
    : BaseType()
    {
        BaseType::data_ = const_cast<pointer>(other.data());
    }

        /** Copy the data (not the pointer) of the rhs.
        */
   TinyVectorView & operator=(TinyVectorView const & r)
    {
        Loop::assign(BaseType::data_, r.begin());
        return *this;
    }

        /** Copy the data of the rhs with cast.
        */
    template <class U, class DATA, class DERIVED>
    TinyVectorView & operator=(TinyArrayBase<U, SIZE, DATA, DERIVED> const & r)
    {
        Loop::assignCast(BaseType::data_, r.begin());
        return *this;
    }
};

/********************************************************/
/*                                                      */
/*                     TinyVector Comparison            */
/*                                                      */
/********************************************************/

/** \addtogroup TinyVectorOperators Functions for TinyVector

    \brief Implement basic arithmetic and equality for TinyVector.

    These functions fulfill the requirements of a Linear Space (vector space).
    Return types are determined according to \ref TinyVectorTraits.

    <b>\#include</b> \<muli/tinyvector.hxx\><br>
    Namespace: muli
*/
//@{
    /// component-wise equal
template <class V1, int SIZE, class D1, class D2, class V2, class D3, class D4>
inline bool
operator==(TinyArrayBase<V1, SIZE, D1, D2> const & l,
           TinyArrayBase<V2, SIZE, D3, D4> const & r)
{
    return !(l != r);
}

    /// component-wise not equal
template <class V1, int SIZE, class D1, class D2, class V2, class D3, class D4>
inline bool
operator!=(TinyArrayBase<V1, SIZE, D1, D2> const & l,
           TinyArrayBase<V2, SIZE, D3, D4> const & r)
{
    typedef typename detail::LoopType<SIZE>::type ltype;
    return ltype::notEqual(l.begin(), r.begin());
}

    /// lexicographical comparison
template <class V1, int SIZE, class D1, class D2, class V2, class D3, class D4>
inline bool
operator<(TinyArrayBase<V1, SIZE, D1, D2> const & l,
                      TinyArrayBase<V2, SIZE, D3, D4> const & r)
{
    typedef typename detail::LoopType<SIZE>::type ltype;
    return ltype::lexicographicLessThan(l.begin(), r.begin());
}


    /// pointwise less-than
template <class V1, int SIZE, class D1, class D2, class V2, class D3, class D4>
inline bool
allLess(TinyArrayBase<V1, SIZE, D1, D2> const & l,
        TinyArrayBase<V2, SIZE, D3, D4> const & r)
{
    for(int k=0; k < SIZE; ++k)
        if (l[k] >= r[k])
            return false;
    return true;
}

    /// pointwise greater-than
template <class V1, int SIZE, class D1, class D2, class V2, class D3, class D4>
inline bool
allGreater(TinyArrayBase<V1, SIZE, D1, D2> const & l,
           TinyArrayBase<V2, SIZE, D3, D4> const & r)
{
    for(int k=0; k < SIZE; ++k)
        if(l[k] <= r[k])
            return false;
    return true;
}

    /// pointwise less-equal
template <class V1, int SIZE, class D1, class D2, class V2, class D3, class D4>
inline bool
allLessEqual(TinyArrayBase<V1, SIZE, D1, D2> const & l,
             TinyArrayBase<V2, SIZE, D3, D4> const & r)
{
    for(int k=0; k < SIZE; ++k)
        if (l[k] > r[k])
            return false;
    return true;
}

    /// pointwise greater-equal
template <class V1, int SIZE, class D1, class D2, class V2, class D3, class D4>
inline bool
allGreaterEqual(TinyArrayBase<V1, SIZE, D1, D2> const & l,
                TinyArrayBase<V2, SIZE, D3, D4> const & r)
{
    for(int k=0; k < SIZE; ++k)
        if (l[k] < r[k])
            return false;
    return true;
}

template <class V, int SIZE, class D1, class D2, class D3, class D4>
bool 
closeAtTolerance(TinyArrayBase<V, SIZE, D1, D2> const & l,
                 TinyArrayBase<V, SIZE, D3, D4> const & r, 
                 V epsilon = NumericTraits<V>::epsilon())
{
    typedef typename detail::LoopType<SIZE>::type ltype;
    return ltype::closeAtTolerance(l.begin(), r.begin(), epsilon);
}

template <class V, int SIZE>
bool 
closeAtTolerance(TinyVector<V, SIZE> const & l,
                 TinyVector<V, SIZE> const & r, 
                 V epsilon = NumericTraits<V>::epsilon())
{
    typedef typename detail::LoopType<SIZE>::type ltype;
    return ltype::closeAtTolerance(l.begin(), r.begin(), epsilon);
}

/********************************************************/
/*                                                      */
/*                     TinyVector Output                */
/*                                                      */
/********************************************************/

    /// stream output
template <class V1, int SIZE, class DATA, class DERIVED>
std::ostream &
operator<<(std::ostream & out, TinyArrayBase<V1, SIZE, DATA, DERIVED> const & l)
{
    out << "(";
    int i;
    for(i=0; i<SIZE-1; ++i)
        out << l[i] << ", ";
    out << l[i] << ")";
    return out;
}
//@}

/********************************************************/
/*                                                      */
/*                      TinyVector-Traits               */
/*                                                      */
/********************************************************/

/** \page TinyVectorTraits Numeric and Promote Traits of TinyVector
    The numeric and promote traits for TinyVectors follow
    the general specifications for \ref NumericPromotionTraits.
    They are implemented in terms of the traits of the basic types by
    partial template specialization:

    \code

    template <class T, int SIZE>
    struct NumericTraits<TinyVector<T, SIZE> >
    {
        typedef TinyVector<typename NumericTraits<T>::Promote, SIZE> Promote;
        typedef TinyVector<typename NumericTraits<T>::RealPromote, SIZE> RealPromote;

        typedef typename NumericTraits<T>::isIntegral isIntegral;
        typedef VigraFalseType isScalar;
        typedef typename NumericTraits<T>::isSigned isSigned;

        // etc.
    };

    template <class T, int SIZE>
    struct NormTraits<TinyVector<T, SIZE> >
    {
        typedef TinyVector<T, SIZE> Type;
        typedef typename Type::SquaredNormType    SquaredNormType;
        typedef typename Type::NormType           NormType;
    };

    template <class T1, class T2, SIZE>
    struct PromoteTraits<TinyVector<T1, SIZE>, TinyVector<T2, SIZE> >
    {
        typedef TinyVector<typename PromoteTraits<T1, T2>::Promote, SIZE> Promote;
    };
    \endcode

    <b>\#include</b> \<muli/tinyvector.hxx\><br>
    Namespace: muli

    On compilers that don't support partial template specialization (e.g.
    MS VisualC++), the traits classes are explicitly specialized for
    <TT>TinyVector<VALUETYPE, SIZE></TT> with
    <TT>VALUETYPE = unsigned char | int | float | double</TT> and <TT>SIZE = 2 | 3 | 4</TT>.

*/

#if !defined(NO_PARTIAL_TEMPLATE_SPECIALIZATION)

template <class T, int SIZE>
struct NumericTraits<TinyVector<T, SIZE> >
{
    typedef TinyVector<T, SIZE> Type;
    typedef TinyVector<typename NumericTraits<T>::Promote, SIZE> Promote;
    typedef TinyVector<typename NumericTraits<T>::RealPromote, SIZE> RealPromote;
    typedef TinyVector<typename NumericTraits<T>::ComplexPromote, SIZE> ComplexPromote;
    typedef T ValueType;

    typedef typename NumericTraits<T>::isIntegral isIntegral;
    typedef VigraFalseType isScalar;
    typedef typename NumericTraits<T>::isSigned isSigned;
    typedef VigraTrueType isOrdered;
    typedef VigraFalseType isComplex;

    static TinyVector<T, SIZE> zero()
    {
        return TinyVector<T, SIZE>(NumericTraits<T>::zero());
    }
    static TinyVector<T, SIZE> one()
    {
        return TinyVector<T, SIZE>(NumericTraits<T>::one());
    }
    static TinyVector<T, SIZE> nonZero()
    {
        return TinyVector<T, SIZE>(NumericTraits<T>::nonZero());
    }

    static TinyVector<T, SIZE> min()
    {
        return TinyVector<T, SIZE>(NumericTraits<T>::min());
    }
    static TinyVector<T, SIZE> max()
    {
        return TinyVector<T, SIZE>(NumericTraits<T>::max());
    }

    template <class D1, class D2>
    static Promote toPromote(TinyArrayBase<T, SIZE, D1, D2> const & v)
    {
        return Promote(v);
    }

    template <class D1, class D2>
    static RealPromote toRealPromote(TinyArrayBase<T, SIZE, D1, D2> const & v)
    {
        return RealPromote(v);
    }

    template <class D1, class D2>
    static TinyVector<T, SIZE>
    fromPromote(TinyArrayBase<typename NumericTraits<T>::Promote, SIZE, D1, D2> const & v)
    {
        TinyVector<T, SIZE> res(detail::dontInit());
        typedef typename detail::LoopType<SIZE>::type ltype;
        ltype::fromPromote(res.begin(), v.begin());
        return res;
    }

    template <class D1, class D2>
    static TinyVector<T, SIZE>
    fromRealPromote(TinyArrayBase<typename NumericTraits<T>::RealPromote, SIZE, D1, D2> const & v)
    {
        TinyVector<T, SIZE> res(detail::dontInit());
        typedef typename detail::LoopType<SIZE>::type ltype;
        ltype::fromRealPromote(res.begin(), v.begin());
        return res;
    }
};

template <class T, int SIZE>
struct NumericTraits<TinyVectorView<T, SIZE> >
: public NumericTraits<TinyVector<T, SIZE> >
{
    typedef TinyVector<T, SIZE> Type;
    typedef TinyVector<typename NumericTraits<T>::Promote, SIZE> Promote;
    typedef TinyVector<typename NumericTraits<T>::RealPromote, SIZE> RealPromote;
    typedef TinyVector<typename NumericTraits<T>::ComplexPromote, SIZE> ComplexPromote;
    typedef T ValueType;

    typedef typename NumericTraits<T>::isIntegral isIntegral;
    typedef VigraFalseType isScalar;
    typedef typename NumericTraits<T>::isSigned isSigned;
    typedef VigraFalseType isOrdered;
    typedef VigraFalseType isComplex;
};

template <class T, int SIZE>
struct NormTraits<TinyVector<T, SIZE> >
{
    typedef TinyVector<T, SIZE> Type;
    typedef typename Type::SquaredNormType    SquaredNormType;
    typedef typename Type::NormType           NormType;
};

template <class T, int SIZE>
struct NormTraits<TinyVectorView<T, SIZE> >
{
    typedef TinyVector<T, SIZE> Type;
    typedef typename Type::SquaredNormType    SquaredNormType;
    typedef typename Type::NormType           NormType;
};

template <class T1, class T2, int SIZE>
struct PromoteTraits<TinyVector<T1, SIZE>, TinyVector<T2, SIZE> >
{
    typedef TinyVector<typename PromoteTraits<T1, T2>::Promote, SIZE> Promote;
};

template <class T1, class T2, int SIZE>
struct PromoteTraits<TinyVectorView<T1, SIZE>, TinyVectorView<T2, SIZE> >
{
    typedef TinyVector<typename PromoteTraits<T1, T2>::Promote, SIZE> Promote;
};

template <class T1, class T2, int SIZE>
struct PromoteTraits<TinyVectorView<T1, SIZE>, TinyVector<T2, SIZE> >
{
    typedef TinyVector<typename PromoteTraits<T1, T2>::Promote, SIZE> Promote;
};

template <class T1, class T2, int SIZE>
struct PromoteTraits<TinyVector<T1, SIZE>, TinyVectorView<T2, SIZE> >
{
    typedef TinyVector<typename PromoteTraits<T1, T2>::Promote, SIZE> Promote;
};

template <class T, int SIZE>
struct PromoteTraits<TinyVector<T, SIZE>, double >
{
    typedef TinyVector<typename NumericTraits<T>::RealPromote, SIZE> Promote;
};

template <class T, int SIZE>
struct PromoteTraits<double, TinyVector<T, SIZE> >
{
    typedef TinyVector<typename NumericTraits<T>::RealPromote, SIZE> Promote;
};

template <class T, int SIZE>
struct PromoteTraits<TinyVectorView<T, SIZE>, double >
{
    typedef TinyVector<typename NumericTraits<T>::RealPromote, SIZE> Promote;
};

template <class T, int SIZE>
struct PromoteTraits<double, TinyVectorView<T, SIZE> >
{
    typedef TinyVector<typename NumericTraits<T>::RealPromote, SIZE> Promote;
};

template<class T, int SIZE>
struct CanSkipInitialization<TinyVectorView<T, SIZE> >
{
    typedef typename CanSkipInitialization<T>::type type;
    static const bool value = type::asBool;
};

template<class T, int SIZE>
struct CanSkipInitialization<TinyVector<T, SIZE> >
{
    typedef typename CanSkipInitialization<T>::type type;
    static const bool value = type::asBool;
};



#else // NO_PARTIAL_TEMPLATE_SPECIALIZATION


#define TINYVECTOR_NUMTRAITS(T, SIZE) \
template<>\
struct NumericTraits<TinyVector<T, SIZE> >\
{\
    typedef TinyVector<T, SIZE> Type;\
    typedef TinyVector<NumericTraits<T>::Promote, SIZE> Promote;\
    typedef TinyVector<NumericTraits<T>::RealPromote, SIZE> RealPromote;\
    typedef TinyVector<NumericTraits<T>::ComplexPromote, SIZE> ComplexPromote;\
    typedef T ValueType; \
    typedef NumericTraits<T>::isIntegral isIntegral;\
    typedef VigraFalseType isScalar;\
    typedef NumericTraits<T>::isSigned isSigned; \
    typedef VigraFalseType isOrdered;\
    typedef VigraFalseType isComplex;\
    \
    static TinyVector<T, SIZE> zero() { \
        return TinyVector<T, SIZE>(NumericTraits<T>::zero()); \
    }\
    static TinyVector<T, SIZE> one() { \
        return TinyVector<T, SIZE>(NumericTraits<T>::one()); \
    }\
    static TinyVector<T, SIZE> nonZero() { \
        return TinyVector<T, SIZE>(NumericTraits<T>::nonZero()); \
    }\
    \
    static Promote toPromote(TinyVector<T, SIZE> const & v) { \
        return Promote(v); \
    }\
    static RealPromote toRealPromote(TinyVector<T, SIZE> const & v) { \
        return RealPromote(v); \
    }\
    static TinyVector<T, SIZE> fromPromote(Promote const & v) { \
        TinyVector<T, SIZE> res;\
        TinyVector<T, SIZE>::iterator d = res.begin(), dend = res.end();\
        Promote::const_iterator s = v.begin();\
        for(; d != dend; ++d, ++s)\
            *d = NumericTraits<T>::fromPromote(*s);\
        return res;\
    }\
    static TinyVector<T, SIZE> fromRealPromote(RealPromote const & v) {\
        TinyVector<T, SIZE> res;\
        TinyVector<T, SIZE>::iterator d = res.begin(), dend = res.end();\
        RealPromote::const_iterator s = v.begin();\
        for(; d != dend; ++d, ++s)\
            *d = NumericTraits<T>::fromRealPromote(*s);\
        return res;\
    }\
}; \
template<>\
struct NormTraits<TinyVector<T, SIZE> >\
{\
    typedef TinyVector<T, SIZE> Type;\
    typedef Type::SquaredNormType           SquaredNormType; \
    typedef Type::NormType NormType; \
};

#define TINYVECTOR_PROMTRAITS1(type1, SIZE) \
template<> \
struct PromoteTraits<TinyVector<type1, SIZE>, TinyVector<type1, SIZE> > \
{ \
    typedef TinyVector<PromoteTraits<type1, type1>::Promote, SIZE> Promote; \
    static Promote toPromote(TinyVector<type1, SIZE> const & v) { \
        return static_cast<Promote>(v); } \
};

#define TINYVECTOR_PROMTRAITS2(type1, type2, SIZE) \
template<> \
struct PromoteTraits<TinyVector<type1, SIZE>, TinyVector<type2, SIZE> > \
{ \
    typedef TinyVector<PromoteTraits<type1, type2>::Promote, SIZE> Promote; \
    static Promote toPromote(TinyVector<type1, SIZE> const & v) { \
        return static_cast<Promote>(v); } \
    static Promote toPromote(TinyVector<type2, SIZE> const & v) { \
       return static_cast<Promote>(v); } \
};

#define TINYVECTOR_TRAITS(SIZE) \
TINYVECTOR_NUMTRAITS(unsigned char, SIZE)\
TINYVECTOR_NUMTRAITS(int, SIZE)\
TINYVECTOR_NUMTRAITS(float, SIZE)\
TINYVECTOR_NUMTRAITS(double, SIZE)\
TINYVECTOR_PROMTRAITS1(unsigned char, SIZE)\
TINYVECTOR_PROMTRAITS1(int, SIZE)\
TINYVECTOR_PROMTRAITS1(float, SIZE)\
TINYVECTOR_PROMTRAITS1(double, SIZE)\
TINYVECTOR_PROMTRAITS2(float, unsigned char, SIZE)\
TINYVECTOR_PROMTRAITS2(unsigned char, float, SIZE)\
TINYVECTOR_PROMTRAITS2(int, unsigned char, SIZE)\
TINYVECTOR_PROMTRAITS2(unsigned char, int, SIZE)\
TINYVECTOR_PROMTRAITS2(int, float, SIZE)\
TINYVECTOR_PROMTRAITS2(float, int, SIZE)\
TINYVECTOR_PROMTRAITS2(double, unsigned char, SIZE)\
TINYVECTOR_PROMTRAITS2(unsigned char, double, SIZE)\
TINYVECTOR_PROMTRAITS2(int, double, SIZE)\
TINYVECTOR_PROMTRAITS2(double, int, SIZE)\
TINYVECTOR_PROMTRAITS2(double, float, SIZE)\
TINYVECTOR_PROMTRAITS2(float, double, SIZE)

TINYVECTOR_TRAITS(2)
TINYVECTOR_TRAITS(3)
TINYVECTOR_TRAITS(4)

#undef TINYVECTOR_NUMTRAITS
#undef TINYVECTOR_PROMTRAITS1
#undef TINYVECTOR_PROMTRAITS2
#undef TINYVECTOR_TRAITS

#endif // NO_PARTIAL_TEMPLATE_SPECIALIZATION


/********************************************************/
/*                                                      */
/*                      TinyVector-Arithmetic           */
/*                                                      */
/********************************************************/

/** \addtogroup TinyVectorOperators
 */
//@{

    /// component-wise addition
template <class V1, int SIZE, class D1, class D2, class V2, class D3, class D4>
inline
typename PromoteTraits<TinyVector<V1, SIZE>, TinyVector<V2, SIZE> >::Promote
operator+(TinyArrayBase<V1, SIZE, D1, D2> const & l,
          TinyArrayBase<V2, SIZE, D3, D4> const & r)
{
    return typename PromoteTraits<TinyVector<V1, SIZE>, TinyVector<V2 , SIZE> >::Promote(l) += r;
}

    /// component-wise subtraction
template <class V1, int SIZE, class D1, class D2, class V2, class D3, class D4>
inline
typename PromoteTraits<TinyVector<V1, SIZE>, TinyVector<V2, SIZE> >::Promote
operator-(TinyArrayBase<V1, SIZE, D1, D2> const & l,
          TinyArrayBase<V2, SIZE, D3, D4> const & r)
{
    return typename PromoteTraits<TinyVector<V1, SIZE>, TinyVector<V2 , SIZE> >::Promote(l) -= r;
}

    /// component-wise multiplication
template <class V1, int SIZE, class D1, class D2, class V2, class D3, class D4>
inline
typename PromoteTraits<TinyVector<V1, SIZE>, TinyVector<V2, SIZE> >::Promote
operator*(TinyArrayBase<V1, SIZE, D1, D2> const & l,
          TinyArrayBase<V2, SIZE, D3, D4> const & r)
{
    return typename PromoteTraits<TinyVector<V1, SIZE>, TinyVector<V2 , SIZE> >::Promote(l) *= r;
}

    /// component-wise division
template <class V1, int SIZE, class D1, class D2, class V2, class D3, class D4>
inline
typename PromoteTraits<TinyVector<V1, SIZE>, TinyVector<V2, SIZE> >::Promote
operator/(TinyArrayBase<V1, SIZE, D1, D2> const & l,
          TinyArrayBase<V2, SIZE, D3, D4> const & r)
{
    return typename PromoteTraits<TinyVector<V1, SIZE>, TinyVector<V2 , SIZE> >::Promote(l) /= r;
}

    /// component-wise modulo
template <class V1, int SIZE, class D1, class D2, class V2, class D3, class D4>
inline
typename PromoteTraits<TinyVector<V1, SIZE>, TinyVector<V2, SIZE> >::Promote
operator%(TinyArrayBase<V1, SIZE, D1, D2> const & l,
          TinyArrayBase<V2, SIZE, D3, D4> const & r)
{
    return typename PromoteTraits<TinyVector<V1, SIZE>, TinyVector<V2 , SIZE> >::Promote(l) %= r;
}

    /// component-wise left scalar addition
template <class V, int SIZE, class D1, class D2>
inline
typename NumericTraits<TinyVector<V, SIZE> >::RealPromote
operator+(double v, TinyArrayBase<V, SIZE, D1, D2> const & r)
{
    return typename NumericTraits<TinyVector<V, SIZE> >::RealPromote(r) += v;
}

    /// component-wise right scalar addition
template <class V, int SIZE, class D1, class D2>
inline
typename NumericTraits<TinyVector<V, SIZE> >::RealPromote
operator+(TinyArrayBase<V, SIZE, D1, D2> const & l, double v)
{
    return typename NumericTraits<TinyVector<V, SIZE> >::RealPromote(l) += v;
}

    /// component-wise left scalar subtraction
template <class V, int SIZE, class D1, class D2>
inline
typename NumericTraits<TinyVector<V, SIZE> >::RealPromote
operator-(double v, TinyArrayBase<V, SIZE, D1, D2> const & r)
{
    return typename NumericTraits<TinyVector<V, SIZE> >::RealPromote(v) -= r;
}

    /// component-wise right scalar subtraction
template <class V, int SIZE, class D1, class D2>
inline
typename NumericTraits<TinyVector<V, SIZE> >::RealPromote
operator-(TinyArrayBase<V, SIZE, D1, D2> const & l, double v)
{
    return typename NumericTraits<TinyVector<V, SIZE> >::RealPromote(l) -= v;
}

    /// component-wise left scalar multiplication
template <class V, int SIZE, class D1, class D2>
inline
typename NumericTraits<TinyVector<V, SIZE> >::RealPromote
operator*(double v, TinyArrayBase<V, SIZE, D1, D2> const & r)
{
    return typename NumericTraits<TinyVector<V, SIZE> >::RealPromote(r) *= v;
}

    /// component-wise right scalar multiplication
template <class V, int SIZE, class D1, class D2>
inline
typename NumericTraits<TinyVector<V, SIZE> >::RealPromote
operator*(TinyArrayBase<V, SIZE, D1, D2> const & l, double v)
{
    return typename NumericTraits<TinyVector<V, SIZE> >::RealPromote(l) *= v;
}

    /// component-wise left scalar division
template <class V, int SIZE, class D1, class D2>
inline
typename NumericTraits<TinyVector<V, SIZE> >::RealPromote
operator/(double v, TinyArrayBase<V, SIZE, D1, D2> const & r)
{
    return typename NumericTraits<TinyVector<V, SIZE> >::RealPromote(v) /= r;
}

    /// component-wise right scalar division
template <class V, int SIZE, class D1, class D2>
inline
typename NumericTraits<TinyVector<V, SIZE> >::RealPromote
operator/(TinyArrayBase<V, SIZE, D1, D2> const & l, double v)
{
    return typename NumericTraits<TinyVector<V, SIZE> >::RealPromote(l) /= v;
}

    /// component-wise scalar division without type promotion
template <class V, int SIZE, class D1, class D2>
inline
TinyVector<V, SIZE>
div(TinyArrayBase<V, SIZE, D1, D2> const & l, V v)
{
    TinyVector<V, SIZE> result(l);
    typedef typename detail::LoopType<SIZE>::type Loop;
    Loop::divScalar(result.data(), v);
    return result;
}


    /** Unary negation (construct TinyVector with negative values)
    */
template <class V, int SIZE, class D1, class D2>
inline
TinyVector<V, SIZE>
operator-(TinyArrayBase<V, SIZE, D1, D2> const & v)
{
    TinyVector<V, SIZE> res(detail::dontInit());
    typedef typename detail::LoopType<SIZE>::type ltype;
    ltype::neg(res.begin(), v.begin());
    return res;
}

    /// component-wise absolute value
template <class V, int SIZE, class D1, class D2>
inline
TinyVector<V, SIZE>
abs(TinyArrayBase<V, SIZE, D1, D2> const & v)
{
    TinyVector<V, SIZE> res(detail::dontInit());
    typedef typename detail::LoopType<SIZE>::type ltype;
    ltype::abs(res.begin(), v.begin());
    return res;
}

    /** Apply ceil() function to each vector component.
    */
template <class V, int SIZE, class D1, class D2>
inline
TinyVector<V, SIZE>
ceil(TinyArrayBase<V, SIZE, D1, D2> const & v)
{
    TinyVector<V, SIZE> res(detail::dontInit());
    typedef typename detail::LoopType<SIZE>::type ltype;
    ltype::ceil(res.begin(), v.begin());
    return res;
}

    /** Apply floor() function to each vector component.
    */
template <class V, int SIZE, class D1, class D2>
inline
TinyVector<V, SIZE>
floor(TinyArrayBase<V, SIZE, D1, D2> const & v)
{
    TinyVector<V, SIZE> res(detail::dontInit());
    typedef typename detail::LoopType<SIZE>::type ltype;
    ltype::floor(res.begin(), v.begin());
    return res;
}

    /** Apply round() function to each vector component.
    */
template <class V, int SIZE, class D1, class D2>
inline
TinyVector<V, SIZE>
round(TinyArrayBase<V, SIZE, D1, D2> const & v)
{
    TinyVector<V, SIZE> res(detail::dontInit());
    typedef typename detail::LoopType<SIZE>::type ltype;
    ltype::round(res.begin(), v.begin());
    return res;
}

    /** Apply roundi() function to each vector component, i.e. return an integer vector.
    */
template <class V, int SIZE, class D1, class D2>
inline
TinyVector<std::ptrdiff_t, SIZE>
roundi(TinyArrayBase<V, SIZE, D1, D2> const & v)
{
    TinyVector<V, SIZE> res(detail::dontInit());
    for(int k=0; k<SIZE; ++k)
        res[k] = roundi(v[k]);
    return res;
}

    /** Apply sqrt() function to each vector component.
    */
template <class V, int SIZE, class D1, class D2>
inline
TinyVector<V, SIZE>
sqrt(TinyArrayBase<V, SIZE, D1, D2> const & v)
{
    TinyVector<V, SIZE> res(detail::dontInit());
    typedef typename detail::LoopType<SIZE>::type ltype;
    ltype::sqrt(res.begin(), v.begin());
    return res;
}

using std::pow;

    /** Apply pow() function to each vector component.
    */
template <class V, int SIZE, class D1, class D2, class E>
inline
TinyVector<V, SIZE>
pow(TinyArrayBase<V, SIZE, D1, D2> const & v, E exponent)
{
    TinyVector<V, SIZE> res(v);
    typedef typename detail::LoopType<SIZE>::type ltype;
    ltype::power(res.begin(), exponent);
    return res;
}

    /// cross product
template <class V1, class D1, class D2, class V2, class D3, class D4>
inline
TinyVector<typename PromoteTraits<V1, V2>::Promote, 3>
cross(TinyArrayBase<V1, 3, D1, D2> const & r1,
      TinyArrayBase<V2, 3, D3, D4> const & r2)
{
    typedef TinyVector<typename PromoteTraits<V1, V2>::Promote, 3>
            Res;
    return  Res(r1[1]*r2[2] - r1[2]*r2[1],
                r1[2]*r2[0] - r1[0]*r2[2],
                r1[0]*r2[1] - r1[1]*r2[0]);
}

    /// dot product
template <class V1, int SIZE, class D1, class D2, class V2, class D3, class D4>
inline
typename PromoteTraits<V1, V2>::Promote
dot(TinyArrayBase<V1, SIZE, D1, D2> const & l,
    TinyArrayBase<V2, SIZE, D3, D4> const & r)
{
    typedef typename detail::LoopType<SIZE>::type ltype;
    return ltype::dot(l.begin(), r.begin());
}

    /// sum of the vector's elements
template <class V, int SIZE, class D1, class D2>
inline
typename NumericTraits<V>::Promote
sum(TinyArrayBase<V, SIZE, D1, D2> const & l)
{
    typename NumericTraits<V>::Promote res = l[0];
    for(int k=1; k<SIZE; ++k)
        res += l[k];
    return res;
}

    /// cumulative sum of the vector's elements
template <class V, int SIZE, class D1, class D2>
inline
TinyVector<typename NumericTraits<V>::Promote, SIZE>
cumsum(TinyArrayBase<V, SIZE, D1, D2> const & l)
{
    TinyVector<typename NumericTraits<V>::Promote, SIZE> res(l);
    for(int k=1; k<SIZE; ++k)
        res[k] += res[k-1];
    return res;
}

    /// product of the vector's elements
template <class V, int SIZE, class D1, class D2>
inline
typename NumericTraits<V>::Promote
prod(TinyArrayBase<V, SIZE, D1, D2> const & l)
{
    typename NumericTraits<V>::Promote res = l[0];
    for(int k=1; k<SIZE; ++k)
        res *= l[k];
    return res;
}

    /// cumulative product of the vector's elements
template <class V, int SIZE, class D1, class D2>
inline
TinyVector<typename NumericTraits<V>::Promote, SIZE>
cumprod(TinyArrayBase<V, SIZE, D1, D2> const & l)
{
    TinyVector<typename NumericTraits<V>::Promote, SIZE> res(l);
    for(int k=1; k<SIZE; ++k)
        res[k] *= res[k-1];
    return res;
}

using std::min;

    /// element-wise minimum
template <class V1, int SIZE, class D1, class D2, class V2, class D3, class D4>
inline
TinyVector<typename PromoteTraits<V1, V2>::Promote, SIZE>
min(TinyArrayBase<V1, SIZE, D1, D2> const & l,
    TinyArrayBase<V2, SIZE, D3, D4> const & r)
{
    typedef typename detail::LoopType<SIZE>::type ltype;
    TinyVector<typename PromoteTraits<V1, V2>::Promote, SIZE> res(l);
    ltype::min(res.begin(), r.begin());
    return res;
}

// we also have to overload min for like-typed argument to prevent match of std::min()
template <class V1, int SIZE, class D1, class D2>
inline
TinyVector<V1, SIZE>
min(TinyArrayBase<V1, SIZE, D1, D2> const & l,
    TinyArrayBase<V1, SIZE, D1, D2> const & r)
{
    typedef typename detail::LoopType<SIZE>::type ltype;
    TinyVector<V1, SIZE> res(l);
    ltype::min(res.begin(), r.begin());
    return res;
}

template <class V1, int SIZE>
inline
TinyVector<V1, SIZE>
min(TinyVector<V1, SIZE> const & l,
    TinyVector<V1, SIZE> const & r)
{
    typedef typename detail::LoopType<SIZE>::type ltype;
    TinyVector<V1, SIZE> res(l);
    ltype::min(res.begin(), r.begin());
    return res;
}

    /// minimum element
template <class V, int SIZE, class D1, class D2>
inline
V const &
min(TinyArrayBase<V, SIZE, D1, D2> const & l)
{
    return l.minimum();
}

using std::max;

    /// element-wise maximum
template <class V1, int SIZE, class D1, class D2, class V2, class D3, class D4>
inline
TinyVector<typename PromoteTraits<V1, V2>::Promote, SIZE>
max(TinyArrayBase<V1, SIZE, D1, D2> const & l,
    TinyArrayBase<V2, SIZE, D3, D4> const & r)
{
    typedef typename detail::LoopType<SIZE>::type ltype;
    TinyVector<typename PromoteTraits<V1, V2>::Promote, SIZE> res(l);
    ltype::max(res.begin(), r.begin());
    return res;
}

// we also have to overload max for like-typed argument to prevent match of std::max()
template <class V1, int SIZE, class D1, class D2>
inline
TinyVector<V1, SIZE>
max(TinyArrayBase<V1, SIZE, D1, D2> const & l,
    TinyArrayBase<V1, SIZE, D1, D2> const & r)
{
    typedef typename detail::LoopType<SIZE>::type ltype;
    TinyVector<V1, SIZE> res(l);
    ltype::max(res.begin(), r.begin());
    return res;
}

template <class V1, int SIZE>
inline
TinyVector<V1, SIZE>
max(TinyVector<V1, SIZE> const & l,
    TinyVector<V1, SIZE> const & r)
{
    typedef typename detail::LoopType<SIZE>::type ltype;
    TinyVector<V1, SIZE> res(l);
    ltype::max(res.begin(), r.begin());
    return res;
}

    /// maximum element
template <class V, int SIZE, class D1, class D2>
inline
V const &
max(TinyArrayBase<V, SIZE, D1, D2> const & l)
{
    return l.maximum();
}

    /// squared norm
template <class V1, int SIZE, class D1, class D2>
inline
typename TinyArrayBase<V1, SIZE, D1, D2>::SquaredNormType
squaredNorm(TinyArrayBase<V1, SIZE, D1, D2> const & t)
{
    return t.squaredMagnitude();
}

    /// squared norm
template <class V, int SIZE>
inline
typename TinyVector<V, SIZE>::SquaredNormType
squaredNorm(TinyVector<V, SIZE> const & t)
{
    return t.squaredMagnitude();
}

using std::reverse;

    /// reversed copy
template <class V, int SIZE>
inline
TinyVector<V, SIZE>
reverse(TinyVector<V, SIZE> const & t)
{
    return TinyVector<V, SIZE>(t.begin(), TinyVector<V, SIZE>::ReverseCopy);
}

    /** \brief transposed copy
    
        Elements are arranged such that <tt>res[k] = t[permutation[k]]</tt>.
    */
template <class V, int SIZE, class T>
inline
TinyVector<V, SIZE>
transpose(TinyVector<V, SIZE> const & t, TinyVector<T, SIZE> const & permutation)
{
    TinyVector<V, SIZE> res(SkipInitialization);
    for(int k=0; k<SIZE; ++k)
    {
        MULI_ASSERT_INSIDE(permutation[k]);
        res[k] = t[permutation[k]];
    }
    return res;
}

    /** \brief Clip negative values.
    
        All elements smaller than 0 are set to zero.
    */
template<class V,int SIZE>
inline
TinyVector<V, SIZE> clipLower(TinyVector<V, SIZE> const & t){
    return clipLower(t, V(0));
}

    /** \brief Clip values below a threshold.
    
        All elements smaller than \a val are set to \a val.
    */
template<class V,int SIZE>
inline
TinyVector<V, SIZE> clipLower(TinyVector<V, SIZE> const & t,const V val){
    TinyVector<V, SIZE> res(SkipInitialization);
    for(int k=0; k<SIZE; ++k){
        res[k]=t[k]< val ? val :  t[k];
    }
    return res;
}

    /** \brief Clip values above a threshold.
    
        All elements bigger than \a val are set to \a val.
    */
template<class V,int SIZE>
inline
TinyVector<V, SIZE> clipUpper(TinyVector<V, SIZE> const & t,const V val){
    TinyVector<V, SIZE> res(SkipInitialization);
    for(int k=0; k<SIZE; ++k){
        res[k]=t[k]> val ? val :  t[k];
    }
    return res;
}

    /** \brief Clip values to an interval.
    
        All elements less than \a valLower are set to \a valLower, all elements
        bigger than \a valUpper are set to \a valUpper.
    */
template<class V,int SIZE>
inline
TinyVector<V, SIZE> clip(TinyVector<V, SIZE> const & t,const V valLower, const V valUpper){
    TinyVector<V, SIZE> res(SkipInitialization);
    for(int k=0; k<SIZE; ++k){
        res[k] =  (t[k] < valLower)
                       ? valLower 
                       : (t[k] > valUpper)
                             ? valUpper
                             : t[k];
    }
    return res;
}

    /** \brief Clip values to a vector of intervals.
    
        All elements less than \a valLower are set to \a valLower, all elements
        bigger than \a valUpper are set to \a valUpper.
    */
template<class V,int SIZE>
inline
TinyVector<V, SIZE> clip(TinyVector<V, SIZE> const & t,
                         TinyVector<V, SIZE> const & valLower, 
                         TinyVector<V, SIZE> const & valUpper)
{
    TinyVector<V, SIZE> res(SkipInitialization);
    for(int k=0; k<SIZE; ++k){
        res[k] =  (t[k] < valLower[k])
                       ? valLower[k] 
                       : (t[k] > valUpper[k])
                             ? valUpper[k] 
                             : t[k];
    }
    return res;
}

template<class V,int SIZE>
inline
bool isZero(TinyVector<V, SIZE> const & t){
    for(int k=0; k<SIZE; ++k){
        if(t[k]!=static_cast<V>(0))
            return false;
    }
    return true;
}

template<class V,int SIZE>
inline typename NumericTraits<V>::RealPromote
mean(TinyVector<V, SIZE> const & t){
    const V sumVal = sum(t);
    return static_cast< typename NumericTraits<V>::RealPromote>(sumVal)/SIZE;
}


template<class V,int SIZE>
inline typename NumericTraits<V>::RealPromote
sizeDividedSquaredNorm(TinyVector<V, SIZE> const & t){
    return squaredNorm(t)/SIZE;
}

template<class V,int SIZE>
inline typename NumericTraits<V>::RealPromote
sizeDividedNorm(TinyVector<V, SIZE> const & t){
    return norm(t)/SIZE;
}



//@}

#endif // #if 0

// mask cl.exe shortcomings [end]
#if defined(_MSC_VER)
#pragma warning( pop )
#endif

} // namespace muli
#undef MULI_ASSERT_INSIDE


#endif // MULI_TINYARRAY_HXX
