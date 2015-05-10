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
#include "config.hxx"
#include "numeric_traits.hxx"
// #include "error.hxx"
// #include "mathutil.hxx"

#ifdef MULI_CHECK_BOUNDS
#define MULI_ASSERT_INSIDE(diff) \
  muli_precondition(diff >= 0, "Index out of bounds");\
  muli_precondition(diff < TinySize<N...>::value, "Index out of bounds");
#else
#define MULI_ASSERT_INSIDE(diff)
#endif

namespace muli {

// mask cl.exe shortcomings [begin]
#if defined(_MSC_VER)
#pragma warning( push )
#pragma warning( disable : 4503 )
#endif

using std::sqrt;

using ArrayIndex = std::ptrdiff_t;

enum SkipInitialization { DontInit };
enum ReverseCopyTag { ReverseCopy };

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

template <ArrayIndex ... N>
struct TinySize
{
    static const ArrayIndex value = TinyShapeImpl<0, N...>::total_size;
    static const ArrayIndex ndim  = sizeof...(N);
};

template <class VALUETYPE, ArrayIndex ... N>
class TinyArray;

template <class VALUETYPE, ArrayIndex ... N>
class TinyArrayView;

/********************************************************/
/*                                                      */
/*                    TinyArrayBase                     */
/*                                                      */
/********************************************************/

/** \brief Base class for fixed size vectors and matrices.

    This class contains functionality shared by
    \ref TinyArray and \ref TinyArrayView, and enables these classes
    to be freely mixed within expressions. It is typically not used directly.

    <b>\#include</b> \<muli/tinyarray.hxx\><br>
    Namespace: muli
**/
template <class VALUETYPE, class DERIVED, ArrayIndex ... N>
class TinyArrayBase
{
  protected:
    using ShapeHelper = TinyShapeImpl<0, N...>;
    
    static const bool derived_is_view = std::is_same<DERIVED, TinyArrayView<VALUETYPE, N...> >::value;
    using data_array_type = typename std::conditional<derived_is_view, 
                                                VALUETYPE *, 
                                                VALUETYPE[ShapeHelper::total_size]>::type;
    
    template <int LEVEL, class V1, class ... V2>
    void initImpl(V1 v1, V2... v2)
    {
        data_[LEVEL] = v1;
        initImpl<LEVEL+1>(v2...);
    }
    
    template <int LEVEL, class V1>
    void initImpl(V1 v1)
    {
        data_[LEVEL] = v1;
    }
    
  public:
  
    template <class NEW_VALUETYPE>
    using AsType = TinyArray<NEW_VALUETYPE, N...>;
    
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
    
  protected:
  
    TinyArrayBase(SkipInitialization)
    {}
    
    // constructores to be used by TinyArray
    
    template <class OTHER, class OTHER_DERIVED>
    TinyArrayBase(TinyArrayBase<OTHER, OTHER_DERIVED, N...> const & other)
    {
        for(ArrayIndex i=0; i<static_size; ++i)
            data_[i] = other[i];
    }
    
    // constructor for zero or one argument
    // activate 'constexpr' in C++ 14
    /* constexpr */ TinyArrayBase(value_type v = value_type())
    {
        for(ArrayIndex i=0; i<static_size; ++i)
            data_[i] = v;
    }
    
    // constructor for two or arguments
    template <class ... V>
    constexpr TinyArrayBase(value_type v0, value_type v1, V ... v)
    : data_{v0, v1, v...}
    {
        static_assert(sizeof...(V)+2 == static_size, 
                      "TinyArrayBase(): wrong number of arguments.");
    }
    
    constexpr TinyArrayBase(value_type const (&v)[static_ndim])
    : data_{v}
    {}
    
    template <class U>
    explicit TinyArrayBase(U const * u)
    {
        for(ArrayIndex i=0; i<static_size; ++i)
            data_[i] = u[i];
    }
    
    template <class U>
    explicit TinyArrayBase(U const * u, ReverseCopyTag)
    {
        for(ArrayIndex i=0; i<static_size; ++i)
            data_[i] = u[static_size-1-i];
    }

  public:
    
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
    
    DERIVED & init(value_type v = value_type())
    {
        for(ArrayIndex i=0; i<static_size; ++i)
            data_[i] = v;
        return static_cast<DERIVED &>(*this);
    }
    
    template <class ... V>
    DERIVED & init(value_type v0, value_type v1, V... v)
    {
        static_assert(sizeof...(V)+2 == static_size, 
                      "TinyArrayBase::init(): wrong number of arguments.");
        initImpl<0>(v0, v1, v...);
        return static_cast<DERIVED &>(*this);
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
    
    RealPromote<value_type> norm() const
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
    
        /** Check that all elements of this vector are zero (or 'false' if T is bool).
        */
    bool allZero() const
    {
        for(ArrayIndex i=0; i<static_size; ++i)
            if(data_[i] != value_type())
                return false;
        return true;
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
    
        /** Get a view to the subarray with length <tt>(TO-FROM)</tt> starting at <tt>FROM</tt>.
            The bounds must fullfill <tt>0 <= FROM < TO <= static_size</tt>.
            Only available if <tt>static_ndim == 1</tt>.
        */
    template <ArrayIndex FROM, ArrayIndex TO>
    typename std::enable_if<static_ndim == 1,
                            TinyArrayView<value_type, TO-FROM> >::type
    subarray() const
    {
        static_assert(FROM >= 0 && FROM < TO && TO <= static_size, 
                      "TinyArrayBase::subarray(): range out of bounds.");
        return TinyArrayView<value_type, TO-FROM>(data_+FROM);
    }
    
    typename std::enable_if<sizeof...(N),
                            TinyArray<value_type, static_size-1> >::type
    dropIndex(ArrayIndex m) const
    {
        TinyArray<value_type, static_size-1> res(DontInit);
        for(ArrayIndex k=0; k<m; ++k)
            res[k] = data_[k];
        for(ArrayIndex k=m; k<static_size-1; ++k)
            res[k] = data_[k+1];
        return res;
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
    
    constexpr bool       empty() const { return static_size == 0; }
    constexpr ArrayIndex size()  const { return static_size; }
    constexpr ArrayIndex max_size()  const { return static_size; }
    constexpr index_type shape() const { return static_shape; }
    constexpr ArrayIndex ndim()  const { return static_ndim; }
    
    template <class OTHER, class OTHER_DERIVED>
    void swap(TinyArrayBase<OTHER, OTHER_DERIVED, N...> & other)
    {
        for(int k=0; k<static_size; ++k)
        {
            Promote<value_type, OTHER> t = data_[k];
            data_[k] = static_cast<value_type>(other[k]);
            other[k] = static_cast<OTHER>(t);
        }
    }
    
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
    o << "{";
    if(TinyArrayBase<T, DERIVED, N...>::static_size > 0)
        o << v[0];
    for(ArrayIndex i=1; i < TinyArrayBase<T, DERIVED, N...>::static_size; ++i)
        o << ", " << v[i];
    o << "}";
    return o;
}

/** \brief Class for fixed size arrays.
    \ingroup RangesAndPoints

    This class contains an array of the specified VALUETYPE with
    (possibly multi-dimensional) shape given by the sequence <tt>ArrayIndex ... N</tt>.
    The interface conforms to STL vector, except that there are no functions
    that change the size of a TinyArray.

    \ref TinyArrayOperators "Arithmetic operations"
    on TinyArrays are defined as component-wise applications of these
    operations.

    See also:<br>
    <UL style="list-style-image:url(documents/bullet.gif)">
        <LI> \ref muli::TinyArrayBase
        <LI> \ref muli::TinyArrayView
        <LI> \ref TinyArrayOperators
    </UL>

    <b>\#include</b> \<vigra/TinyArray.hxx\><br>
    Namespace: vigra
**/
template <class VALUETYPE, ArrayIndex ... N>
class TinyArray
: public TinyArrayBase<VALUETYPE, TinyArray<VALUETYPE, N...>, N...>
{
    using BaseType = TinyArrayBase<VALUETYPE, TinyArray<VALUETYPE, N...>, N...>;
    
  public:
  
    typedef typename BaseType::value_type value_type;
    static const ArrayIndex static_ndim = BaseType::static_ndim;
  
    constexpr TinyArray(TinyArray const &) = default;
    
    template <class OTHER, class DERIVED>
    TinyArray(TinyArrayBase<OTHER, DERIVED, N...> const & other)
    : BaseType(other)
    {}
    
    TinyArray(SkipInitialization)
    : BaseType(DontInit)
    {}
    
        /** Construction from lemon::Invalid.
        
            Initializes all vector elements with -1.
        */
    TinyArray(lemon::Invalid const &)
    : BaseType(-1)
    {}
    
    TinyArray(value_type v = value_type())
    : BaseType(v)
    {}
    
    template <class ... V>
    constexpr TinyArray(value_type v0, value_type v1, V... v)
    : BaseType(v0, v1, v...)
    {}
    
    constexpr TinyArray(value_type const (&v)[static_ndim])
    : BaseType(v)
    {}
    
    template <class U>
    explicit TinyArray(U const * u)
    : BaseType(u)
    {}
    
    template <class U>
    explicit TinyArray(U const * u, ReverseCopyTag)
    : BaseType(u, ReverseCopy)
    {}

    TinyArray & operator=(TinyArray const &) = default;

    TinyArray & operator=(value_type v)
    {
        BaseType::operator=(v);
        return *this;
    }
    
    TinyArray & operator=(value_type const (&v)[static_ndim])
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

    /// factory function for fixed-size unit matrix
template <class VALUETYPE, ArrayIndex N>
inline 
TinyArray<VALUETYPE, N, N>
tinyEye()
{
    TinyArray<VALUETYPE, N, N> res;
    for(ArrayIndex k=0; k<N; ++k)
        res(k,k) = 1;
    return res;
}

    /// factory function for the fixed-size k-th unit vector 
template <class VALUETYPE, ArrayIndex N>
inline 
TinyArray<VALUETYPE, N>
unitVector(ArrayIndex k)
{
    TinyArray<VALUETYPE, N> res;
    res(k) = 1;
    return res;
}

    /// factory function for fixed-size linear sequence starting at <tt>start</tt> with stepsize <tt>step</tt>
template <class VALUETYPE, ArrayIndex ... N>
inline 
TinyArray<VALUETYPE, N...>
linearSequence(VALUETYPE start = VALUETYPE(), VALUETYPE step = VALUETYPE(1))
{
    TinyArray<VALUETYPE, N...> res;
    for(ArrayIndex k=0; k < TinySize<N...>::value; ++k)
        res[k] = k;
    return res;
}

/** \brief Wrapper for fixed size arrays.

    This class wraps the mamory of an array of the specified VALUETYPE
    with (possibly multi-dimensional) shape given by <tt>ArrayIndex....N</tt>.
    Thus, the array can be accessed with an interface similar to
    that of std::vector (except that there are no functions
    that change the size of a TinyArrayView). The TinyArrayView
    does <em>not</em> assume ownership of the given memory.

    \ref TinyArrayOperators "Arithmetic operations"
    on TinyArrayViews are defined as component-wise applications of these
    operations. 

    <b>See also:</b>
    <ul>
        <li> \ref muli::TinyArrayBase
        <li> \ref muli::TinyArray
        <li> \ref TinyArrayOperators
    </ul>

    <b>\#include</b> \<muli/tinyarray.hxx\><br>
    Namespace: muli
**/
template <class VALUETYPE, ArrayIndex ... N>
class TinyArrayView
: public TinyArrayBase<VALUETYPE, TinyArrayView<VALUETYPE, N...>, N...>
{
    using BaseType = TinyArrayBase<VALUETYPE, TinyArrayView<VALUETYPE, N...>, N...>;
    
  public:
  
    typedef typename BaseType::value_type value_type;
    typedef typename BaseType::pointer pointer;
    typedef typename BaseType::const_pointer const_pointer;
    static const ArrayIndex static_size = BaseType::static_size;
    static const ArrayIndex static_ndim = BaseType::static_ndim;
  
    TinyArrayView()
    : BaseType(DontInit)
    {
        BaseType::data_ = nullptr;
    }

        /** Construct view for given data array
        */
    TinyArrayView(const_pointer data)
    : BaseType(DontInit)
    {
        BaseType::data_ = const_cast<pointer>(data);
    }

        /** Copy constructor (shallow copy).
        */
    TinyArrayView(TinyArrayView const & other)
    : BaseType(DontInit)
    {
        BaseType::data_ = const_cast<pointer>(other.data());
    }

        /** Construct view from other TinyArray.
        */
    template <class OTHER_DERIVED>
    TinyArrayView(TinyArrayBase<value_type, OTHER_DERIVED, N...> const & other)
    : BaseType(DontInit)
    {
        BaseType::data_ = const_cast<pointer>(other.data());
    }

        /** Copy the data (not the pointer) of the rhs.
        */
    TinyArrayView & operator=(TinyArrayView const & r)
    {
        for(int k=0; k<static_size; ++k)
            BaseType::data_[k] = r[k];
        return *this;
    }

        /** Copy the data of the rhs with cast.
        */
    template <class U, class OTHER_DERIVED>
    TinyArrayView & operator=(TinyArrayBase<U, OTHER_DERIVED, N...> const & r)
    {
        for(int k=0; k<static_size; ++k)
            BaseType::data_[k] = static_cast<value_type>(r[k]);
        return *this;
    }
};

/********************************************************/
/*                                                      */
/*                TinyArray Comparison                  */
/*                                                      */
/********************************************************/

/** \addtogroup TinyArrayOperators Functions for TinyArray

    \brief Implement basic arithmetic and equality for TinyArray.

    These functions fulinit the requirements of a Linear Space (vector space).
    Return types are determined according to \ref TinyArrayTraits.

    <b>\#include</b> \<muli/TinyArray.hxx\><br>
    Namespace: muli
*/
//@{
    /// element-wise equal
template <class V1, class D1, class V2, class D2, ArrayIndex ... N>
inline bool
operator==(TinyArrayBase<V1, D1, N...> const & l,
           TinyArrayBase<V2, D2, N...> const & r)
{
    for(ArrayIndex k=0; k < TinySize<N...>::value; ++k)
        if(l[k] != r[k])
            return false;
    return true;
}

    /// element-wise not equal
template <class V1, class D1, class V2, class D2, ArrayIndex ... N>
inline bool
operator!=(TinyArrayBase<V1, D1, N...> const & l,
           TinyArrayBase<V2, D2, N...> const & r)
{
    for(ArrayIndex k=0; k < TinySize<N...>::value; ++k)
        if(l[k] != r[k])
            return true;
    return false;
}

    /// lexicographical comparison
template <class V1, class D1, class V2, class D2, ArrayIndex ... N>
inline bool
operator<(TinyArrayBase<V1, D1, N...> const & l,
          TinyArrayBase<V2, D2, N...> const & r)
{
    for(ArrayIndex k=0; k < TinySize<N...>::value; ++k)
    {
        if(l[k] < r[k])
            return true;
        if(r[k] < l[k])
            return false;
    }
    return false;
}

    /// check if all elements are zero
template <class V, class D, ArrayIndex ... N>
inline bool
isZero(TinyArrayBase<V, D, N...> const & t)
{
    return t.allZero();
}

    /// pointwise less-than
template <class V1, class D1, class V2, class D2, ArrayIndex ... N>
inline bool
allLess(TinyArrayBase<V1, D1, N...> const & l,
        TinyArrayBase<V2, D2, N...> const & r)
{
    for(int k=0; k < TinySize<N...>::value; ++k)
        if (l[k] >= r[k])
            return false;
    return true;
}

    /// pointwise greater-than
template <class V1, class D1, class V2, class D2, ArrayIndex ... N>
inline bool
allGreater(TinyArrayBase<V1, D1, N...> const & l,
           TinyArrayBase<V2, D2, N...> const & r)
{
    for(int k=0; k < TinySize<N...>::value; ++k)
        if(l[k] <= r[k])
            return false;
    return true;
}

    /// pointwise less-equal
template <class V1, class D1, class V2, class D2, ArrayIndex ... N>
inline bool
allLessEqual(TinyArrayBase<V1, D1, N...> const & l,
             TinyArrayBase<V2, D2, N...> const & r)
{
    for(int k=0; k < TinySize<N...>::value; ++k)
        if (l[k] > r[k])
            return false;
    return true;
}

    /// pointwise greater-equal
template <class V1, class D1, class V2, class D2, ArrayIndex ... N>
inline bool
allGreaterEqual(TinyArrayBase<V1, D1, N...> const & l,
                TinyArrayBase<V2, D2, N...> const & r)
{
    for(int k=0; k < TinySize<N...>::value; ++k)
        if (l[k] < r[k])
            return false;
    return true;
}

template <class V1, class D1, class V2, class D2, ArrayIndex ... N>
inline bool
closeAtTolerance(TinyArrayBase<V1, D1, N...> const & l,
                 TinyArrayBase<V2, D2, N...> const & r, 
                 Promote<V1, V2> epsilon = 2.0*NumericTraits<Promote<V1, V2> >::epsilon())
{
    for(ArrayIndex k=0; k < TinySize<N...>::value; ++k)
        if(!closeAtTolerance(l[k], r[k], epsilon))
            return false;
    return true;
}

/********************************************************/
/*                                                      */
/*                 TinyArray-Arithmetic                 */
/*                                                      */
/********************************************************/

/** \addtogroup TinyArrayOperators
 */
//@{

    /// element-wise addition
template <class V1, class D1, class V2, class D2, ArrayIndex ... N>
inline
TinyArray<Promote<V1, V2>, N...>
operator+(TinyArrayBase<V1, D1, N...> const & l,
          TinyArrayBase<V2, D2, N...> const & r)
{
    return TinyArray<Promote<V1, V2>, N...>(l) += r;
}

    /// element-wise scalar addition
template <class V1, class D1, class V2, ArrayIndex ... N>
inline
TinyArray<PromoteIf<std::is_arithmetic<V2>::value, V1, V2>, N...>
operator+(TinyArrayBase<V1, D1, N...> const & l,
          V2 r)
{
    return TinyArray<Promote<V1, V2>, N...>(l) += r;
}

    /// element-wise left scalar addition
template <class V1, class V2, class D2, ArrayIndex ... N>
inline
TinyArray<PromoteIf<std::is_arithmetic<V1>::value, V1, V2>, N...>
operator+(V1 l,
          TinyArrayBase<V2, D2, N...> const & r)
{
    return TinyArray<Promote<V1, V2>, N...>(l) += r;
}

    /// element-wise subtraction
template <class V1, class D1, class V2, class D2, ArrayIndex ... N>
inline
TinyArray<Promote<V1, V2>, N...>
operator-(TinyArrayBase<V1, D1, N...> const & l,
          TinyArrayBase<V2, D2, N...> const & r)
{
    return TinyArray<Promote<V1, V2>, N...>(l) -= r;
}

    /// element-wise scalar subtraction
template <class V1, class D1, class V2, ArrayIndex ... N>
inline
TinyArray<PromoteIf<std::is_arithmetic<V2>::value, V1, V2>, N...>
operator-(TinyArrayBase<V1, D1, N...> const & l,
          V2 r)
{
    return TinyArray<Promote<V1, V2>, N...>(l) -= r;
}

    /// element-wise left scalar subtraction
template <class V1, class V2, class D2, ArrayIndex ... N>
inline
TinyArray<PromoteIf<std::is_arithmetic<V1>::value, V1, V2>, N...>
operator-(V1 l,
          TinyArrayBase<V2, D2, N...> const & r)
{
    return TinyArray<Promote<V1, V2>, N...>(l) -= r;
}

    /// Unary negation
template <class V, class D, ArrayIndex ... N>
inline
TinyArray<V, N...>
operator-(TinyArrayBase<V, D, N...> const & v)
{
    TinyArray<V, N...> res(DontInit);
    for(ArrayIndex k=0; k < TinySize<N...>::value; ++k)
        res[k] = -v[k];
    return res;
}

    /// element-wise multiplication
template <class V1, class D1, class V2, class D2, ArrayIndex ... N>
inline
TinyArray<Promote<V1, V2>, N...>
operator*(TinyArrayBase<V1, D1, N...> const & l,
          TinyArrayBase<V2, D2, N...> const & r)
{
    return TinyArray<Promote<V1, V2>, N...>(l) *= r;
}

    /// element-wise scalar multiplication
template <class V1, class D1, class V2, ArrayIndex ... N>
inline
TinyArray<PromoteIf<std::is_arithmetic<V2>::value, V1, V2>, N...>
operator*(TinyArrayBase<V1, D1, N...> const & l,
          V2 r)
{
    return TinyArray<Promote<V1, V2>, N...>(l) *= r;
}

    /// element-wise left scalar multiplication
template <class V1, class V2, class D2, ArrayIndex ... N>
inline
TinyArray<PromoteIf<std::is_arithmetic<V1>::value, V1, V2>, N...>
operator*(V1 l,
          TinyArrayBase<V2, D2, N...> const & r)
{
    return TinyArray<Promote<V1, V2>, N...>(l) *= r;
}

    /// element-wise division
template <class V1, class D1, class V2, class D2, ArrayIndex ... N>
inline
TinyArray<Promote<V1, V2>, N...>
operator/(TinyArrayBase<V1, D1, N...> const & l,
          TinyArrayBase<V2, D2, N...> const & r)
{
    return TinyArray<Promote<V1, V2>, N...>(l) /= r;
}

    /// element-wise scalar division
template <class V1, class D1, class V2, ArrayIndex ... N>
inline
TinyArray<PromoteIf<std::is_arithmetic<V2>::value, V1, V2>, N...>
operator/(TinyArrayBase<V1, D1, N...> const & l,
          V2 r)
{
    return TinyArray<Promote<V1, V2>, N...>(l) /= r;
}

    /// element-wise left scalar division
template <class V1, class V2, class D2, ArrayIndex ... N>
inline
TinyArray<PromoteIf<std::is_arithmetic<V1>::value, V1, V2>, N...>
operator/(V1 l,
          TinyArrayBase<V2, D2, N...> const & r)
{
    return TinyArray<Promote<V1, V2>, N...>(l) /= r;
}

    /// element-wise scalar division without type promotion
template <class V1, class D1, class V2, ArrayIndex ... N>
inline
typename std::enable_if<std::is_arithmetic<V2>::value,
                        TinyArray<V1, N...> >::type
div(TinyArrayBase<V1, D1, N...> const & l, V2 r)
{
    return TinyArray<V1, N...>(l) /= r;
}

    /// element-wise modulo
template <class V1, class D1, class V2, class D2, ArrayIndex ... N>
inline
TinyArray<Promote<V1, V2>, N...>
operator%(TinyArrayBase<V1, D1, N...> const & l,
          TinyArrayBase<V2, D2, N...> const & r)
{
    return TinyArray<Promote<V1, V2>, N...>(l) %= r;
}

    /// element-wise scalar modulo
template <class V1, class D1, class V2, ArrayIndex ... N>
inline
TinyArray<PromoteIf<std::is_arithmetic<V2>::value, V1, V2>, N...>
operator%(TinyArrayBase<V1, D1, N...> const & l,
          V2 r)
{
    return TinyArray<Promote<V1, V2>, N...>(l) %= r;
}

    /// element-wise left scalar modulo
template <class V1, class V2, class D2, ArrayIndex ... N>
inline
TinyArray<PromoteIf<std::is_arithmetic<V1>::value, V1, V2>, N...>
operator%(V1 l,
          TinyArrayBase<V2, D2, N...> const & r)
{
    return TinyArray<Promote<V1, V2>, N...>(l) %= r;
}

using std::abs;

    /// element-wise absolute value
template <class V, class D, ArrayIndex ... N>
inline
TinyArray<V, N...>
abs(TinyArrayBase<V, D, N...> const & v)
{
    TinyArray<V, N...> res(DontInit);
    for(ArrayIndex k=0; k < TinySize<N...>::value; ++k)
        res[k] = abs(v[k]);
    return res;
}

using std::ceil;

    /** Apply ceil() function to each vector component.
    */
template <class V, class D, ArrayIndex ... N>
inline
TinyArray<V, N...>
ceil(TinyArrayBase<V, D, N...> const & v)
{
    TinyArray<V, N...> res(DontInit);
    for(ArrayIndex k=0; k < TinySize<N...>::value; ++k)
        res[k] = ceil(v[k]);
    return res;
}

using std::floor;

    /** Apply floor() function to each vector component.
    */
template <class V, class D, ArrayIndex ... N>
inline
TinyArray<V, N...>
floor(TinyArrayBase<V, D, N...> const & v)
{
    TinyArray<V, N...> res(DontInit);
    for(ArrayIndex k=0; k < TinySize<N...>::value; ++k)
        res[k] = floor(v[k]);
    return res;
}

using std::round;

    /** Apply round() function to each vector component.
    */
template <class V, class D, ArrayIndex ... N>
inline
TinyArray<V, N...>
round(TinyArrayBase<V, D, N...> const & v)
{
    TinyArray<V, N...> res(DontInit);
    for(ArrayIndex k=0; k < TinySize<N...>::value; ++k)
        res[k] = round(v[k]);
    return res;
}

    // /** Apply roundi() function to each vector component, i.e. return an integer vector.
    // */
// template <class V, class D, ArrayIndex ... N>
// inline
// TinyArray<ArrayIndex, SIZE>
// roundi(TinyArrayBase<V, D, N...> const & v)
// {
    // TinyArray<V, N...> res(DontInit);
    // for(ArrayIndex k=0; k < TinySize<N...>::value; ++k)
        // res[k] = roundi(v[k]);
    // return res;
// }

using std::sqrt;

    /** Apply sqrt() function to each vector component.
    */
template <class V, class D, ArrayIndex ... N>
inline
TinyArray<RealPromote<V>, N...>
sqrt(TinyArrayBase<V, D, N...> const & v)
{
    TinyArray<RealPromote<V>, N...> res(DontInit);
    for(ArrayIndex k=0; k < TinySize<N...>::value; ++k)
        res[k] = sqrt(v[k]);
    return res;
}

using std::pow;

    /** Apply pow() function to each vector component.
    */
template <class V, class D, class E, ArrayIndex ... N>
inline
TinyArray<Promote<V, E>, N...>
pow(TinyArrayBase<V, D, N...> const & v, E exponent)
{
    TinyArray<Promote<V, E>, N...> res(DontInit);
    for(ArrayIndex k=0; k < TinySize<N...>::value; ++k)
        res[k] = pow(v[k], exponent);
    return res;
}

    /// cross product
template <class V1, class D1, class V2, class D2>
inline
TinyArray<Promote<V1, V2>, 3>
cross(TinyArrayBase<V1, D1, 3> const & r1,
      TinyArrayBase<V2, D2, 3> const & r2)
{
    typedef TinyArray<Promote<V1, V2>, 3> Res;
    return  Res(r1[1]*r2[2] - r1[2]*r2[1],
                r1[2]*r2[0] - r1[0]*r2[2],
                r1[0]*r2[1] - r1[1]*r2[0]);
}

    /// dot product
template <class V1, class D1, class V2, class D2, ArrayIndex ... N>
inline
Promote<V1, V2>
dot(TinyArrayBase<V1, D1, N...> const & l,
    TinyArrayBase<V2, D2, N...> const & r)
{
    Promote<V1, V2> res = l[0] * r[0];
    for(ArrayIndex k=1; k < TinySize<N...>::value; ++k)
        res += l[k] * r[k];
    return res;
}

    /// sum of the vector's elements
template <class V, class D, ArrayIndex ... N>
inline
Promote<V>
sum(TinyArrayBase<V, D, N...> const & l)
{
    Promote<V> res = l[0];
    for(int k=1; k < TinySize<N...>::value; ++k)
        res += l[k];
    return res;
}

    /// mean of the vector's elements
template <class V, class D, ArrayIndex ... N>
inline RealPromote<V>
mean(TinyArrayBase<V, D, N...> const & t)
{
    const V sumVal = sum(t);
    return static_cast<RealPromote<V> >(sumVal) / TinySize<N...>::value;
}

    /// cumulative sum of the vector's elements
template <class V, class D, ArrayIndex ... N>
inline
TinyArray<Promote<V>, N...>
cumsum(TinyArrayBase<V, D, N...> const & l)
{
    TinyArray<Promote<V>, N...> res(l);
    for(int k=1; k < TinySize<N...>::value; ++k)
        res[k] += res[k-1];
    return res;
}

    /// product of the vector's elements
template <class V, class D, ArrayIndex ... N>
inline
Promote<V>
prod(TinyArrayBase<V, D, N...> const & l)
{
    Promote<V> res = l[0];
    for(int k=1; k < TinySize<N...>::value; ++k)
        res *= l[k];
    return res;
}

    /// cumulative product of the vector's elements
template <class V, class D, ArrayIndex ... N>
inline
TinyArray<Promote<V>, N...>
cumprod(TinyArrayBase<V, D, N...> const & l)
{
    TinyArray<Promote<V>, N...> res(l);
    for(int k=1; k < TinySize<N...>::value; ++k)
        res[k] *= res[k-1];
    return res;
}

using std::min;

    /// element-wise minimum
template <class V1, class D1, class V2, class D2, ArrayIndex ... N>
inline
TinyArray<Promote<V1, V2>, N...>
min(TinyArrayBase<V1, D1, N...> const & l,
    TinyArrayBase<V2, D2, N...> const & r)
{
    TinyArray<Promote<V1, V2>, N...> res(DontInit);
    for(int k=0; k < TinySize<N...>::value; ++k)
        res[k] = min<Promote<V1, V2> >(l[k], r[k]);
    return res;
}

// we also have to overload min for like-typed argument to prevent match of std::min()
template <class V, class D, ArrayIndex ...N>
inline
TinyArray<V, N...>
min(TinyArrayBase<V, D, N...> const & l,
    TinyArrayBase<V, D, N...> const & r)
{
    TinyArray<V, N...> res(DontInit);
    for(int k=0; k < TinySize<N...>::value; ++k)
        res[k] = min(l[k], r[k]);
    return res;
}

    /// minimum element
template <class V, class D, ArrayIndex ... N>
inline
V const &
min(TinyArrayBase<V, D, N...> const & l)
{
    return l.minimum();
}

using std::max;

    /// element-wise maximum
template <class V1, class D1, class V2, class D2, ArrayIndex ... N>
inline
TinyArray<Promote<V1, V2>, N...>
max(TinyArrayBase<V1, D1, N...> const & l,
    TinyArrayBase<V2, D2, N...> const & r)
{
    TinyArray<Promote<V1, V2>, N...> res(DontInit);
    for(int k=0; k < TinySize<N...>::value; ++k)
        res[k] = max<Promote<V1, V2> >(l[k], r[k]);
    return res;
}

// we also have to overload max for like-typed argument to prevent match of std::max()
template <class V, class D, ArrayIndex ...N>
inline
TinyArray<V, N...>
max(TinyArrayBase<V, D, N...> const & l,
    TinyArrayBase<V, D, N...> const & r)
{
    TinyArray<V, N...> res(DontInit);
    for(int k=0; k < TinySize<N...>::value; ++k)
        res[k] = max(l[k], r[k]);
    return res;
}


    /// maximum element
template <class V, class D, ArrayIndex ... N>
inline
V const &
max(TinyArrayBase<V, D, N...> const & l)
{
    return l.maximum();
}

    /// squared norm
template <class V, class D, ArrayIndex ... N>
inline
Promote<V>
squaredNorm(TinyArrayBase<V, D, N...> const & t)
{
    return t.squaredNorm();
}

template <class V, class D, ArrayIndex ... N>
inline RealPromote<V>
sizeDividedSquaredNorm(TinyArrayBase<V, D, N...> const & t)
{
    return static_cast<RealPromote<V> >(t.squaredNorm()) / TinySize<N...>::value;
}

template <class V, class D, ArrayIndex ... N>
inline RealPromote<V>
sizeDividedNorm(TinyArrayBase<V, D, N...> const & t)
{
    return static_cast<RealPromote<V> >(t.norm()) / TinySize<N...>::value;
}

using std::reverse;

    /// reversed copy
template <class V, class D, ArrayIndex ... N>
inline
TinyArray<V, N...>
reverse(TinyArrayBase<V, D, N...> const & t)
{
    return TinyArray<V, N...>(t.begin(), ReverseCopy);
}

    /** \brief transposed copy
    
        Elements are arranged such that <tt>res[k] = t[permutation[k]]</tt>.
    */
template <class V1, class D1, class V2, class D2, ArrayIndex ...N>
inline
TinyArray<V1, N...>
transpose(TinyArrayBase<V1, D1, N...> const & v, 
          TinyArrayBase<V2, D2, N...> const & permutation)
{
    TinyArray<V1, N...> res(DontInit);
    for(int k=0; k < TinySize<N...>::value; ++k)
    {
        MULI_ASSERT_INSIDE(permutation[k]);
        res[k] = v[permutation[k]];
    }
    return res;
}

    /** \brief Clip negative values.
    
        All elements smaller than 0 are set to zero.
    */
template <class V, class D, ArrayIndex ... N>
inline
TinyArray<V, N...> 
clipLower(TinyArrayBase<V, D, N...> const & t)
{
    return clipLower(t, V());
}

    /** \brief Clip values below a threshold.
    
        All elements smaller than \a val are set to \a val.
    */
template <class V, class D, ArrayIndex ... N>
inline
TinyArray<V, N...> 
clipLower(TinyArrayBase<V, D, N...> const & t, const V val)
{
    TinyArray<V, N...> res(DontInit);
    for(int k=0; k < TinySize<N...>::value; ++k)
    {
        res[k] = t[k] < val ? val :  t[k];
    }
    return res;
}

    /** \brief Clip values above a threshold.
    
        All elements bigger than \a val are set to \a val.
    */
template <class V, class D, ArrayIndex ... N>
inline
TinyArray<V, N...> 
clipUpper(TinyArrayBase<V, D, N...> const & t, const V val)
{
    TinyArray<V, N...> res(DontInit);
    for(int k=0; k < TinySize<N...>::value; ++k)
    {
        res[k] = t[k] > val ? val :  t[k];
    }
    return res;
}

    /** \brief Clip values to an interval.
    
        All elements less than \a valLower are set to \a valLower, all elements
        bigger than \a valUpper are set to \a valUpper.
    */
template <class V, class D, ArrayIndex ... N>
inline
TinyArray<V, N...> 
clip(TinyArrayBase<V, D, N...> const & t,
     const V valLower, const V valUpper)
{
    TinyArray<V, N...> res(DontInit);
    for(int k=0; k < TinySize<N...>::value; ++k)
    {
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
template <class V, class D1, class D2, class D3, ArrayIndex ... N>
inline
TinyArray<V, N...> 
clip(TinyArrayBase<V, D1, N...> const & t,
     TinyArrayBase<V, D2, N...> const & valLower, 
     TinyArrayBase<V, D3, N...> const & valUpper)
{
    TinyArray<V, N...> res(DontInit);
    for(int k=0; k < TinySize<N...>::value; ++k)
    {
        res[k] =  (t[k] < valLower[k])
                       ? valLower[k] 
                       : (t[k] > valUpper[k])
                             ? valUpper[k] 
                             : t[k];
    }
    return res;
}

//@}

// mask cl.exe shortcomings [end]
#if defined(_MSC_VER)
#pragma warning( pop )
#endif

} // namespace muli

#undef MULI_ASSERT_INSIDE


#endif // MULI_TINYARRAY_HXX
