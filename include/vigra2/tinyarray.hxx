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

#ifndef VIGRA_TINYARRAY_HXX
#define VIGRA_TINYARRAY_HXX

#include "config.hxx"
#include "numeric_traits.hxx"
#include "error.hxx"
#include "math.hxx"
#include <cmath>    // abs(double)
#include <cstdlib>  // abs(int)
#include <iosfwd>   // ostream
#include <algorithm>
#include <iterator>

#if defined(VIGRA_WITH_LEMON)
    #include <lemon/core.h>
#endif

#ifdef VIGRA_CHECK_BOUNDS
    #define VIGRA_ASSERT_INSIDE(diff) \
      vigra_precondition(diff >= 0, "Index out of bounds");\
      vigra_precondition(diff < TinySize<N...>::value, "Index out of bounds");
#else
    #define VIGRA_ASSERT_INSIDE(diff)
#endif

namespace vigra {

namespace lemon {

#if defined(VIGRA_WITH_LEMON)

using Invalid = ::lemon::Invalid;

#else

struct Invalid 
{
  public:
    bool operator==(Invalid) { return true;  }
    bool operator!=(Invalid) { return false; }
    bool operator< (Invalid) { return false; }
};

#endif

const Invalid INVALID = Invalid();

} // namespace lemon

// mask cl.exe shortcomings [begin]
#if defined(_MSC_VER)
#pragma warning( push )
#pragma warning( disable : 4503 )
#endif

using std::sqrt;

enum SkipInitialization { DontInit };
enum ReverseCopyTag { ReverseCopy };

// FIXME: document this
template <ArrayIndex LEVEL, int ... N>
struct TinyShapeImpl;

template <ArrayIndex LEVEL, int N, int ... REST>
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

template <ArrayIndex LEVEL, int N>
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

template <int ... N>
struct TinySize
{
    static const ArrayIndex value = TinyShapeImpl<0, N...>::total_size;
    static const ArrayIndex ndim  = sizeof...(N);
};

template <class VALUETYPE, int ... N>
class TinyArray;

template <class VALUETYPE, int ... N>
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

    <b>\#include</b> \<vigra/tinyarray.hxx\><br>
    Namespace: vigra
**/
template <class VALUETYPE, class DERIVED, int ... N>
class TinyArrayBase
{
  protected:
    using ShapeHelper = TinyShapeImpl<0, N...>;
    
    static const bool derived_is_array = std::is_same<DERIVED, TinyArray<VALUETYPE, N...> >::value;
    using data_array_type = typename std::conditional<derived_is_array, 
                                                VALUETYPE[ShapeHelper::total_size],
                                                VALUETYPE *>::type;
    
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
    
    static const ArrayIndex static_ndim  = sizeof...(N);
    static const ArrayIndex static_size  = ShapeHelper::total_size;
    static constexpr index_type static_shape = index_type(N...);
    
    // constructors
        
    constexpr TinyArrayBase(TinyArrayBase const &) = default;
    
  protected:
  
    TinyArrayBase(SkipInitialization)
    {
        static_assert(static_size > 0, 
                      "TinyArrayBase(): array must have non-zero size.");
    }
    
    // constructors to be used by TinyArray
    
    template <class OTHER, class OTHER_DERIVED>
    TinyArrayBase(TinyArrayBase<OTHER, OTHER_DERIVED, N...> const & other)
    {
        static_assert(static_size > 0, 
                      "TinyArrayBase(): array must have non-zero size.");
        for(int i=0; i<static_size; ++i)
            data_[i] = detail::RequiresExplicitCast<value_type>::cast(other[i]);
    }
    
    // constructor for zero or one argument
    // activate 'constexpr' in C++ 14
    /* constexpr */ TinyArrayBase(value_type v = value_type())
    {
        static_assert(static_size > 0, 
                      "TinyArrayBase(): array must have non-zero size.");
        for(int i=0; i<static_size; ++i)
            data_[i] = v;
    }
    
    // constructor for two or arguments
    template <class ... V>
    constexpr TinyArrayBase(value_type v0, value_type v1, V ... v)
    : data_{v0, v1, v...}
    {
        static_assert(static_size > 0, 
                      "TinyArrayBase(): array must have non-zero size.");
        static_assert(sizeof...(V)+2 == static_size, 
                      "TinyArrayBase(): wrong number of arguments.");
    }
    
    constexpr TinyArrayBase(value_type const (&v)[static_ndim])
    : data_{v}
    {
        static_assert(static_size > 0, 
                      "TinyArrayBase(): array must have non-zero size.");
    }
    
    template <class U>
    explicit TinyArrayBase(U const * u)
    {
        static_assert(static_size > 0, 
                      "TinyArrayBase(): array must have non-zero size.");
        for(int i=0; i<static_size; ++i)
            data_[i] = detail::RequiresExplicitCast<value_type>::cast(u[i]);
    }
    
    template <class U>
    TinyArrayBase(U const * u, ReverseCopyTag)
    {
        static_assert(static_size > 0, 
                      "TinyArrayBase(): array must have non-zero size.");
        for(int i=0; i<static_size; ++i)
            data_[i] = detail::RequiresExplicitCast<value_type>::cast(u[static_size-1-i]);
    }
    
        // for compatibility with TinyArrayBase<..., runtime_size>
    template <class U>
    TinyArrayBase(U const * u, U const * /* end */, ReverseCopyTag)
    : TinyArrayBase(u, ReverseCopy)
    {}

  public:
    
    // assignment
    
    TinyArrayBase & operator=(TinyArrayBase const &) = default;

    TinyArrayBase & operator=(value_type v)
    {
        for(int i=0; i<static_size; ++i)
            data_[i] = v;
        return *this;
    }
    
    TinyArrayBase & operator=(value_type const (&v)[static_ndim])
    {
        for(int i=0; i<static_size; ++i)
            data_[i] = v[i];
        return *this;
    }
    
    template <class OTHER, class OTHER_DERIVED>
    TinyArrayBase & operator=(TinyArrayBase<OTHER, OTHER_DERIVED, N...> const & other)
    {
        for(int i=0; i<static_size; ++i)
            data_[i] = detail::RequiresExplicitCast<value_type>::cast(other[i]);
        return *this;
    }
    
    DERIVED & init(value_type v = value_type())
    {
        for(int i=0; i<static_size; ++i)
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
    
    template <class Iterator>
    DERIVED & init(Iterator first, Iterator end)
    {
        int range = std::distance(first, end);
        if(static_size < range)
            range = static_size;
        for(int i=0; i<range; ++i, ++first)
            data_[i] = detail::RequiresExplicitCast<value_type>::cast(*first);
        return static_cast<DERIVED &>(*this);
    }
    
    // vector functions
    
    SquaredNormType<value_type> squaredNormImpl() const
    {
        SquaredNormType<value_type> result = SquaredNormType<value_type>();
        for(int i=0; i<static_size; ++i)
            result += squaredNorm(data_[i]);
        return result;
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
    template <int FROM, int TO>
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
        for(int k=0; k<m; ++k)
            res[k] = data_[k];
        for(int k=m; k<static_size-1; ++k)
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
    
    reverse_iterator rbegin() { return reverse_iterator(data_ + static_size); }
    reverse_iterator rend()   { return reverse_iterator(data_); }
    const_reverse_iterator rbegin() const { return const_reverse_iterator(data_ + static_size); }
    const_reverse_iterator rend()   const { return const_reverse_iterator(data_); }
    const_reverse_iterator crbegin() const { return const_reverse_iterator(data_ + static_size); }
    const_reverse_iterator crend()   const { return const_reverse_iterator(data_); }
    
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
    
    TinyArrayBase & reverse()
    {
        ArrayIndex i=0, j=size()-1;
        while(i < j)
             std::swap(data_[i++], data_[j--]);
        return *this;
    }

    template <class OTHER, class OTHER_DERIVED>
    void swap(TinyArrayBase<OTHER, OTHER_DERIVED, N...> & other)
    {
        for(int k=0; k<static_size; ++k)
        {
            PromoteType<value_type, OTHER> t = data_[k];
            data_[k] = static_cast<value_type>(other[k]);
            other[k] = static_cast<OTHER>(t);
        }
    }
    
        /// factory function for fixed-size unit matrix
    template <int SIZE>
    static inline 
    TinyArray<value_type, SIZE, SIZE>
    eye()
    {
        TinyArray<value_type, SIZE, SIZE> res;
        for(int k=0; k<SIZE; ++k)
            res(k,k) = 1;
        return res;
    }

        /// factory function for the fixed-size k-th unit vector 
    template <int SIZE>
    static inline 
    TinyArray<value_type, SIZE>
    unitVector(ArrayIndex k)
    {
        TinyArray<value_type, SIZE> res;
        res(k) = 1;
        return res;
    }

        /// factory function for fixed-size linear sequence starting at <tt>start</tt> with stepsize <tt>step</tt>
    static inline 
    TinyArray<value_type, N...>
    linearSequence(value_type start = value_type(), value_type step = value_type(1))
    {
        TinyArray<value_type, N...> res;
        for(int k=0; k < static_size; ++k, start += step)
            res[k] = start;
        return res;
    }
    
  protected:
    data_array_type data_;
};

template <class T, class DERIVED, int ... N>
constexpr 
typename TinyArrayBase<T, DERIVED, N...>::index_type 
TinyArrayBase<T, DERIVED, N...>::static_shape;

/********************************************************/
/*                                                      */
/*                TinyArrayBase output                  */
/*                                                      */
/********************************************************/

template <class T, class DERIVED, int ... N>
std::ostream & operator<<(std::ostream & o, TinyArrayBase<T, DERIVED, N...> const & v)
{
    o << "{";
    if(v.size() > 0)
        o << v[0];
    for(int i=1; i < v.size(); ++i)
        o << ", " << v[i];
    o << "}";
    return o;
}

template <class T, class DERIVED, int N1, int N2>
std::ostream & operator<<(std::ostream & o, TinyArrayBase<T, DERIVED, N1, N2> const & v)
{
    o << "{";
    for(int i=0; i<N1; ++i)
    {
        if(i > 0)
            o << ",\n ";
        o << v(i,0);
        for(int j=1; j<N2; ++j)
        {
            o << ", " << v(i, j);
        }
    }
    o << "}";
    return o;
}

/********************************************************/
/*                                                      */
/*         TinyArrayBase<..., runtime_size>             */
/*                                                      */
/********************************************************/

/** \brief Specialization of TinyArrayBase for dynamic arrays.

    This class contains functionality shared by
    \ref TinyArray and \ref TinyArrayView, and enables these classes
    to be freely mixed within expressions. It is typically not used directly.

    <b>\#include</b> \<vigra/tinyarray.hxx\><br>
    Namespace: vigra
**/
template <class VALUETYPE, class DERIVED>
class TinyArrayBase<VALUETYPE, DERIVED, runtime_size>
{
  public:
  
    template <class NEW_VALUETYPE>
    using AsType = TinyArray<NEW_VALUETYPE, runtime_size>;
    
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
    
    template <class V>
    DERIVED & init(std::initializer_list<V> l)
    {
        return init(l.begin(), l.end());
    }
        
    // vector functions
    
    SquaredNormType<value_type> squaredNormImpl() const
    {
        SquaredNormType<value_type> result = SquaredNormType<value_type>();
        for(int i=0; i<size_; ++i)
            result += squaredNorm(data_[i]);
        return result;
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
    
        /** Get a view to the subarray with length <tt>(TO-FROM)</tt> starting at <tt>FROM</tt>.
            The bounds must fullfill <tt>0 <= FROM < TO <= size_</tt>.
        */
    TinyArrayBase
    subarray(ArrayIndex FROM, ArrayIndex TO) const
    {
        vigra_precondition(FROM >= 0 && FROM < TO && TO <= size_, 
                      "TinyArrayBase::subarray(): range out of bounds.");
        return TinyArrayBase(TO-FROM, data_+FROM);
    }
    
    TinyArray<value_type, runtime_size>
    dropIndex(ArrayIndex m) const
    {
        TinyArray<value_type, runtime_size> res(size_-1, DontInit);
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
    TinyArray<value_type, runtime_size>
    unitVector(ArrayIndex size, ArrayIndex k)
    {
        TinyArray<value_type, runtime_size> res(size);
        res[k] = 1;
        return res;
    }

        /// factory function for a linear sequence from <tt>start</tt> to <tt>stop</tt> 
        /// (exclusive) with stepsize <tt>step</tt>
    static inline 
    TinyArray<value_type, runtime_size>
    range(value_type start, 
          value_type stop, 
          value_type step = value_type(1))
    {
        TinyArray<value_type, runtime_size> res(stop - start, DontInit);
        for(int k=0; k < res.size(); ++k, start += step)
            res[k] = start;
        return res;
    }

        /// factory function for a linear sequence from 0 to <tt>stop</tt> 
        /// (exclusive) with stepsize 1
    static inline 
    TinyArray<value_type, runtime_size>
    range(value_type stop)
    {
        TinyArray<value_type, runtime_size> res(stop, DontInit);
        value_type start = value_type();
        for(int k=0; k < res.size(); ++k, ++start)
            res[k] = start;
        return res;
    }
    
  protected:
    ArrayIndex size_;
    pointer data_;
};

/********************************************************/
/*                                                      */
/*                       TinyArray                      */
/*                                                      */
/********************************************************/


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
        <LI> \ref vigra::TinyArrayBase
        <LI> \ref vigra::TinyArrayView
        <LI> \ref TinyArrayOperators
    </UL>

    <b>\#include</b> \<vigra/TinyArray.hxx\><br>
    Namespace: vigra
**/
template <class VALUETYPE, int ... N>
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
    
    // for compatibility with TinyArray<VALUETYPE, runtime_size>
    TinyArray(ArrayIndex /* size */, SkipInitialization)
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

        // for compatibility with TinyArray<..., runtime_size>
    template <class U>
    TinyArray(U const * u, U const * /* end */, ReverseCopyTag)
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

/********************************************************/
/*                                                      */
/*            TinyArray<..., runtime_size>              */
/*                                                      */
/********************************************************/

/** \brief Specialization of TinyArray for dynamic arrays.
    \ingroup RangesAndPoints

    This class contains an array of the specified VALUETYPE with
    size specified at runtim.
    The interface conforms to STL vector, except that there are no functions
    that change the size of a TinyArray.

    \ref TinyArrayOperators "Arithmetic operations"
    on TinyArrays are defined as component-wise applications of these
    operations.

    See also:<br>
    <UL style="list-style-image:url(documents/bullet.gif)">
        <LI> \ref vigra::TinyArrayBase
        <LI> \ref vigra::TinyArrayView
        <LI> \ref TinyArrayOperators
    </UL>

    <b>\#include</b> \<vigra/TinyArray.hxx\><br>
    Namespace: vigra
**/
template <class VALUETYPE>
class TinyArray<VALUETYPE, runtime_size>
: public TinyArrayBase<VALUETYPE, TinyArray<VALUETYPE, runtime_size>, runtime_size>
{
    using BaseType = TinyArrayBase<VALUETYPE, TinyArray<VALUETYPE, runtime_size>, runtime_size>;
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
    
    template <class U>
    TinyArray(U const * i, U const * end, ReverseCopyTag)
    : TinyArray(std::distance(i, end), DontInit)
    {
         for(int i=0; i<this->size_; ++i, --end)
            new(this->data_+i) value_type(detail::RequiresExplicitCast<value_type>::cast(*(end-1)));
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

/********************************************************/
/*                                                      */
/*                    TinyArrayView                     */
/*                                                      */
/********************************************************/

/** \brief Wrapper for fixed size arrays.

    This class wraps the memory of an array of the specified VALUETYPE
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
        <li> \ref vigra::TinyArrayBase
        <li> \ref vigra::TinyArray
        <li> \ref vigra::TinySymmetricView
        <li> \ref TinyArrayOperators
    </ul>

    <b>\#include</b> \<vigra/tinyarray.hxx\><br>
    Namespace: vigra
**/
template <class VALUETYPE, int ... N>
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

        /** Reset to the other array's pointer.
        */
    template <class OTHER_DERIVED>
    void reset(TinyArrayBase<value_type, OTHER_DERIVED, N...> const & other)
    {
        BaseType::data_ = const_cast<pointer>(other.data());
    }

        /** Copy the data (not the pointer) of the rhs.
        */
    TinyArrayView & operator=(TinyArrayView const & r)
    {
        for(int k=0; k<static_size; ++k)
            BaseType::data_[k] = detail::RequiresExplicitCast<value_type>::cast(r[k]);
        return *this;
    }

        /** Copy the data of the rhs with cast.
        */
    template <class U, class OTHER_DERIVED>
    TinyArrayView & operator=(TinyArrayBase<U, OTHER_DERIVED, N...> const & r)
    {
        for(int k=0; k<static_size; ++k)
            BaseType::data_[k] = detail::RequiresExplicitCast<value_type>::cast(r[k]);
        return *this;
    }
};

/********************************************************/
/*                                                      */
/*               TinyArraySymmetricView                 */
/*                                                      */
/********************************************************/

/** \brief Wrapper for fixed size arrays.

    This class wraps the memory of an 1D array of the specified VALUETYPE
    with size <tt>N*(N+1)/2</tt> and interprets this array as a symmetric
    matrix. Specifically, the data are interpreted as the row-wise
    representation of the upper triangular part of the symmetric matrix.
    All index access operations are overloaded such that the view appears 
    as if it were a full matrix. The TinySymmetricView
    does <em>not</em> assume ownership of the given memory.

    \ref TinyArrayOperators "Arithmetic operations"
    on TinySymmetricView are defined as component-wise applications of these
    operations. 

    <b>See also:</b>
    <ul>
        <li> \ref vigra::TinyArrayBase
        <li> \ref vigra::TinyArray
        <li> \ref vigra::TinyArrayView
        <li> \ref TinyArrayOperators
    </ul>

    <b>\#include</b> \<vigra/tinyarray.hxx\><br>
    Namespace: vigra
**/
template <class VALUETYPE, int N>
class TinySymmetricView
: public TinyArrayBase<VALUETYPE, TinySymmetricView<VALUETYPE, N>, N*(N+1)/2>
{
    using BaseType = TinyArrayBase<VALUETYPE, TinySymmetricView<VALUETYPE, N>, N*(N+1)/2>;
    
  public:
  
    typedef typename BaseType::value_type value_type;
    typedef typename BaseType::pointer pointer;
    typedef typename BaseType::const_pointer const_pointer;
    typedef typename BaseType::reference reference;
    typedef typename BaseType::const_reference const_reference;
    using index_type             = TinyArray<ArrayIndex, 2>;

    static const ArrayIndex static_size = BaseType::static_size;
    static const ArrayIndex static_ndim = 2;
    static constexpr index_type static_shape = index_type(N, N);
  
    TinySymmetricView()
    : BaseType(DontInit)
    {
        BaseType::data_ = nullptr;
    }

        /** Construct view for given data array
        */
    TinySymmetricView(const_pointer data)
    : BaseType(DontInit)
    {
        BaseType::data_ = const_cast<pointer>(data);
    }

        /** Copy constructor (shallow copy).
        */
    TinySymmetricView(TinySymmetricView const & other)
    : BaseType(DontInit)
    {
        BaseType::data_ = const_cast<pointer>(other.data());
    }

        /** Construct view from other TinyArray.
        */
    template <class OTHER_DERIVED>
    TinySymmetricView(TinyArrayBase<value_type, OTHER_DERIVED, N*(N+1)/2> const & other)
    : BaseType(DontInit)
    {
        BaseType::data_ = const_cast<pointer>(other.data());
    }

        /** Reset to the other array's pointer.
        */
    template <class OTHER_DERIVED>
    void reset(TinyArrayBase<value_type, OTHER_DERIVED, N*(N+1)/2> const & other)
    {
        BaseType::data_ = const_cast<pointer>(other.data());
    }

        /** Copy the data (not the pointer) of the rhs.
        */
    TinySymmetricView & operator=(TinySymmetricView const & r)
    {
        for(int k=0; k<static_size; ++k)
            BaseType::data_[k] = detail::RequiresExplicitCast<value_type>::cast(r[k]);
        return *this;
    }

        /** Copy the data of the rhs with cast.
        */
    template <class U, class OTHER_DERIVED>
    TinySymmetricView & operator=(TinyArrayBase<U, OTHER_DERIVED, N*(N+1)/2> const & r)
    {
        for(int k=0; k<static_size; ++k)
            BaseType::data_[k] = detail::RequiresExplicitCast<value_type>::cast(r[k]);
        return *this;
    }
    
    // index access
    
    reference operator[](ArrayIndex i)
    {
        return BaseType::operator[](i);
    }
    
    constexpr const_reference operator[](ArrayIndex i) const
    {
        return BaseType::operator[](i);
    }
    
    reference at(ArrayIndex i)
    {
        return BaseType::at(i);
    }
    
    const_reference at(ArrayIndex i) const
    {
        return BaseType::at(i);
    }
    
    reference operator[](ArrayIndex const (&i)[2])
    {
        return this->operator()(i[0], i[1]);
    }
    
    constexpr const_reference operator[](ArrayIndex const (&i)[2]) const
    {
        return this->operator()(i[0], i[1]);
    }
    
    reference at(ArrayIndex const (&i)[static_ndim])
    {
        return this->at(i[0], i[1]);
    }
    
    const_reference at(ArrayIndex const (&i)[static_ndim]) const
    {
        return this->at(i[0], i[1]);
    }
    
    reference operator[](index_type const & i)
    {
        return this->operator()(i[0], i[1]);
    }
    
    constexpr const_reference operator[](index_type const & i) const
    {
        return this->operator()(i[0], i[1]);
    }
    
    reference at(index_type const & i)
    {
        return this->at(i[0], i[1]);
    }
    
    const_reference at(index_type const & i) const
    {
        return this->at(i[0], i[1]);
    }
    
    reference operator()(ArrayIndex i, ArrayIndex j)
    {
        if(i > j)
            return BaseType::data_[i + N*j - (j*(j+1) >> 1)];
        else
            return BaseType::data_[N*i + j - (i*(i+1) >> 1)];
    }
    
    constexpr const_reference operator()(ArrayIndex i, ArrayIndex j) const
    {
        if(i > j)
            return BaseType::data_[i + (j*((2*N-1) - j) >> 1)];
        else
            return BaseType::data_[j + (i*((2*N-1) - i) >> 1)];
    }
    
    reference at(ArrayIndex i, ArrayIndex j)
    {
        ArrayIndex k = (i > j)
                           ? i + (j*((2*N-1) - j) >> 1)
                           : j + (i*((2*N-1) - i) >> 1);
        if(k < 0 || k >= static_size)
            throw std::out_of_range("TinySymmetricView::at()");
        return BaseType::data_[k];
    }
    
    const_reference at(ArrayIndex i, ArrayIndex j) const
    {
        ArrayIndex k = (i > j)
                           ? i + N*j - (j*(j+1) >> 1)
                           : N*i + j - (i*(i+1) >> 1);
        if(k < 0 || k >= static_size)
            throw std::out_of_range("TinySymmetricView::at()");
        return BaseType::data_[k];
    }
    
    constexpr index_type shape() const { return static_shape; }
    constexpr ArrayIndex ndim () const { return static_ndim; }
};

template <class T, int N>
constexpr 
typename TinySymmetricView<T, N>::index_type 
TinySymmetricView<T, N>::static_shape;

/********************************************************/
/*                                                      */
/*            TinyArraySymmetricView output             */
/*                                                      */
/********************************************************/

template <class T, int N>
std::ostream & operator<<(std::ostream & o, TinySymmetricView<T, N> const & v)
{
    o << "{";
    for(int i=0; i<N; ++i)
    {
        if(i > 0)
            o << ",\n ";
        o << v(i,0);
        for(int j=1; j<N; ++j)
        {
            o << ", " << v(i, j);
        }
    }
    o << "}";
    return o;
}

namespace detail {

template <int N0, int ... N>
struct TinyArrayIsStatic
{
    static const int first = N0;
    static const int ndim = sizeof...(N)+1;
    static const bool value = ndim > 1 || N0 != runtime_size;
};

} // namespace detail

#define vigra_tinyarray_check_size(PREDICATE, MESSAGE) \
  if(detail::TinyArrayIsStatic<N...>::value) {} else \
    if((PREDICATE)) {} else \
      vigra::throw_contract_error("Precondition violation!", MESSAGE, __FILE__, __LINE__)


/********************************************************/
/*                                                      */
/*                TinyArray Comparison                  */
/*                                                      */
/********************************************************/

/** \addtogroup TinyArrayOperators Functions for TinyArray

    \brief Implement basic arithmetic and equality for TinyArray.

    These functions fulfill the requirements of a Linear Space (vector space).
    Return types are determined according to \ref PromoteType or \ref RealPromoteType.

    <b>\#include</b> \<vigra/TinyArray.hxx\><br>
    Namespace: vigra
*/
//@{
    /// element-wise equal
template <class V1, class D1, class V2, class D2, int ... N>
inline bool
operator==(TinyArrayBase<V1, D1, N...> const & l,
           TinyArrayBase<V2, D2, N...> const & r)
{
    vigra_tinyarray_check_size(l.size() == r.size(),
        "TinyArrayBase::operator==(): size mismatch.");
    for(int k=0; k < l.size(); ++k)
        if(l[k] != r[k])
            return false;
    return true;
}

    /// element-wise not equal
template <class V1, class D1, class V2, class D2, int ... N>
inline bool
operator!=(TinyArrayBase<V1, D1, N...> const & l,
           TinyArrayBase<V2, D2, N...> const & r)
{
    vigra_tinyarray_check_size(l.size() == r.size(),
        "TinyArrayBase::operator!=(): size mismatch.");
    for(int k=0; k < l.size(); ++k)
        if(l[k] != r[k])
            return true;
    return false;
}

    /// lexicographical comparison
template <class V1, class D1, class V2, class D2, int ... N>
inline bool
operator<(TinyArrayBase<V1, D1, N...> const & l,
          TinyArrayBase<V2, D2, N...> const & r)
{
    vigra_tinyarray_check_size(l.size() == r.size(),
        "TinyArrayBase::operator<(): size mismatch.");
    for(int k=0; k < l.size(); ++k)
    {
        if(l[k] < r[k])
            return true;
        if(r[k] < l[k])
            return false;
    }
    return false;
}

    /// check if all elements are non-zero (or 'true' if V is bool)
template <class V, class D, int ... N>
inline bool
all(TinyArrayBase<V, D, N...> const & t)
{
    for(int i=0; i<t.size(); ++i)
        if(t[i] == V())
            return false;
    return true;
}

    /// check if at least one element is non-zero (or 'true' if V is bool)
template <class V, class D, int ... N>
inline bool
any(TinyArrayBase<V, D, N...> const & t)
{
    for(int i=0; i<t.size(); ++i)
        if(t[i] != V())
            return true;
    return false;
}

    /// check if all elements are zero (or 'false' if V is bool)
template <class V, class D, int ... N>
inline bool
allZero(TinyArrayBase<V, D, N...> const & t)
{
    for(int i=0; i<t.size(); ++i)
        if(t[i] != V())
            return false;
    return true;
}

    /// pointwise less-than
template <class V1, class D1, class V2, class D2, int ... N>
inline bool
allLess(TinyArrayBase<V1, D1, N...> const & l,
        TinyArrayBase<V2, D2, N...> const & r)
{
    vigra_tinyarray_check_size(l.size() == r.size(),
        "TinyArrayBase::allLess(): size mismatch.");
    for(int k=0; k < l.size(); ++k)
        if (l[k] >= r[k])
            return false;
    return true;
}

    /// pointwise greater-than
template <class V1, class D1, class V2, class D2, int ... N>
inline bool
allGreater(TinyArrayBase<V1, D1, N...> const & l,
           TinyArrayBase<V2, D2, N...> const & r)
{
    vigra_tinyarray_check_size(l.size() == r.size(),
        "TinyArrayBase::allGreater(): size mismatch.");
    for(int k=0; k < l.size(); ++k)
        if(l[k] <= r[k])
            return false;
    return true;
}

    /// pointwise less-equal
template <class V1, class D1, class V2, class D2, int ... N>
inline bool
allLessEqual(TinyArrayBase<V1, D1, N...> const & l,
             TinyArrayBase<V2, D2, N...> const & r)
{
    vigra_tinyarray_check_size(l.size() == r.size(),
        "TinyArrayBase::allLessEqual(): size mismatch.");
    for(int k=0; k < l.size(); ++k)
        if (l[k] > r[k])
            return false;
    return true;
}

    /// pointwise greater-equal
template <class V1, class D1, class V2, class D2, int ... N>
inline bool
allGreaterEqual(TinyArrayBase<V1, D1, N...> const & l,
                TinyArrayBase<V2, D2, N...> const & r)
{
    vigra_tinyarray_check_size(l.size() == r.size(),
        "TinyArrayBase::allGreaterEqual(): size mismatch.");
    for(int k=0; k < l.size(); ++k)
        if (l[k] < r[k])
            return false;
    return true;
}

template <class V1, class D1, class V2, class D2, int ... N>
inline bool
closeAtTolerance(TinyArrayBase<V1, D1, N...> const & l,
                 TinyArrayBase<V2, D2, N...> const & r, 
                 PromoteType<V1, V2> epsilon = 2.0*NumericTraits<PromoteType<V1, V2> >::epsilon())
{
    vigra_tinyarray_check_size(l.size() == r.size(),
        "TinyArrayBase::closeAtTolerance(): size mismatch.");
    for(int k=0; k < l.size(); ++k)
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


// Declare arithmetic functions for documentation,
// the implementations are provided by a macro below.

    /// scalar add-assignment
template <class VALUE, class DERIVED, int ... N>
DERIVED & 
operator+=(TinyArrayBase<VALUE, DERIVED, N...> & l,
           VALUE r);

    /// element-wise add-assignment
template <class VALUE, class DERIVED, class OTHER, class OTHER_DERIVED, int ... N>
DERIVED & 
operator+=(TinyArrayBase<VALUE, DERIVED, N...> & l,
           TinyArrayBase<OTHER, OTHER_DERIVED, N...> const & r);
           
    /// element-wise addition
template <class V1, class D1, class V2, class D2, int ... N>
TinyArray<PromoteType<V1, V2>, N...>
operator+(TinyArrayBase<V1, D1, N...> const & l,
          TinyArrayBase<V2, D2, N...> const & r);

    /// element-wise scalar addition
template <class V1, class D1, class V2, int ... N>
TinyArray<PromoteTypeIf<std::is_arithmetic<V2>::value, V1, V2>, N...>
operator+(TinyArrayBase<V1, D1, N...> const & l,
          V2 r);

    /// element-wise left scalar addition
template <class V1, class V2, class D2, int ... N>
TinyArray<PromoteTypeIf<std::is_arithmetic<V1>::value, V1, V2>, N...>
operator+(V1 l,
          TinyArrayBase<V2, D2, N...> const & r);

    /// scalar subtract-assignment
template <class VALUE, class DERIVED, int ... N>
DERIVED & 
operator-=(TinyArrayBase<VALUE, DERIVED, N...> & l,
           VALUE r);

    /// element-wise subtract-assignment
template <class VALUE, class DERIVED, class OTHER, class OTHER_DERIVED, int ... N>
DERIVED & 
operator-=(TinyArrayBase<VALUE, DERIVED, N...> & l,
           TinyArrayBase<OTHER, OTHER_DERIVED, N...> const & r);
           
    /// element-wise subtraction
template <class V1, class D1, class V2, class D2, int ... N>
TinyArray<PromoteType<V1, V2>, N...>
operator-(TinyArrayBase<V1, D1, N...> const & l,
          TinyArrayBase<V2, D2, N...> const & r);

    /// element-wise scalar subtraction
template <class V1, class D1, class V2, int ... N>
TinyArray<PromoteTypeIf<std::is_arithmetic<V2>::value, V1, V2>, N...>
operator-(TinyArrayBase<V1, D1, N...> const & l,
          V2 r);

    /// element-wise left scalar subtraction
template <class V1, class V2, class D2, int ... N>
TinyArray<PromoteTypeIf<std::is_arithmetic<V1>::value, V1, V2>, N...>
operator-(V1 l,
          TinyArrayBase<V2, D2, N...> const & r);

    /// scalar multiply-assignment
template <class VALUE, class DERIVED, int ... N>
DERIVED & 
operator*=(TinyArrayBase<VALUE, DERIVED, N...> & l,
           VALUE r);

    /// element-wise multiply-assignment
template <class VALUE, class DERIVED, class OTHER, class OTHER_DERIVED, int ... N>
DERIVED & 
operator*=(TinyArrayBase<VALUE, DERIVED, N...> & l,
           TinyArrayBase<OTHER, OTHER_DERIVED, N...> const & r);
           
    /// element-wise multiplication
template <class V1, class D1, class V2, class D2, int ... N>
TinyArray<PromoteType<V1, V2>, N...>
operator*(TinyArrayBase<V1, D1, N...> const & l,
          TinyArrayBase<V2, D2, N...> const & r);

    /// element-wise scalar multiplication
template <class V1, class D1, class V2, int ... N>
TinyArray<PromoteTypeIf<std::is_arithmetic<V2>::value, V1, V2>, N...>
operator*(TinyArrayBase<V1, D1, N...> const & l,
          V2 r);

    /// element-wise left scalar multiplication
template <class V1, class V2, class D2, int ... N>
TinyArray<PromoteTypeIf<std::is_arithmetic<V1>::value, V1, V2>, N...>
operator*(V1 l,
          TinyArrayBase<V2, D2, N...> const & r);

    /// scalar divide-assignment
template <class VALUE, class DERIVED, int ... N>
DERIVED & 
operator/=(TinyArrayBase<VALUE, DERIVED, N...> & l,
           VALUE r);

    /// element-wise divide-assignment
template <class VALUE, class DERIVED, class OTHER, class OTHER_DERIVED, int ... N>
DERIVED & 
operator/=(TinyArrayBase<VALUE, DERIVED, N...> & l,
           TinyArrayBase<OTHER, OTHER_DERIVED, N...> const & r);
           
    /// element-wise division
template <class V1, class D1, class V2, class D2, int ... N>
TinyArray<PromoteType<V1, V2>, N...>
operator/(TinyArrayBase<V1, D1, N...> const & l,
          TinyArrayBase<V2, D2, N...> const & r);

    /// element-wise scalar division
template <class V1, class D1, class V2, int ... N>
TinyArray<PromoteTypeIf<std::is_arithmetic<V2>::value, V1, V2>, N...>
operator/(TinyArrayBase<V1, D1, N...> const & l,
          V2 r);

    /// element-wise left scalar division
template <class V1, class V2, class D2, int ... N>
TinyArray<PromoteTypeIf<std::is_arithmetic<V1>::value, V1, V2>, N...>
operator/(V1 l,
          TinyArrayBase<V2, D2, N...> const & r);

    /// scalar modulo-assignment
template <class VALUE, class DERIVED, int ... N>
DERIVED & 
operator%=(TinyArrayBase<VALUE, DERIVED, N...> & l,
           VALUE r);

    /// element-wise modulo-assignment
template <class VALUE, class DERIVED, class OTHER, class OTHER_DERIVED, int ... N>
DERIVED & 
operator%=(TinyArrayBase<VALUE, DERIVED, N...> & l,
           TinyArrayBase<OTHER, OTHER_DERIVED, N...> const & r);
           
    /// element-wise modulo
template <class V1, class D1, class V2, class D2, int ... N>
TinyArray<PromoteType<V1, V2>, N...>
operator%(TinyArrayBase<V1, D1, N...> const & l,
          TinyArrayBase<V2, D2, N...> const & r);

    /// element-wise scalar modulo
template <class V1, class D1, class V2, int ... N>
TinyArray<PromoteTypeIf<std::is_arithmetic<V2>::value, V1, V2>, N...>
operator%(TinyArrayBase<V1, D1, N...> const & l,
          V2 r);

    /// element-wise left scalar modulo
template <class V1, class V2, class D2, int ... N>
TinyArray<PromoteTypeIf<std::is_arithmetic<V1>::value, V1, V2>, N...>
operator%(V1 l,
          TinyArrayBase<V2, D2, N...> const & r);

#define VIGRA_TINYARRAY_OPERATORS(op) \
template <class VALUE, class DERIVED, int ... N> \
inline DERIVED &  \
operator op##=(TinyArrayBase<VALUE, DERIVED, N...> & l, \
               VALUE r) \
{ \
    for(int i=0; i<l.size(); ++i) \
        l[i] op##= r; \
    return static_cast<DERIVED &>(l); \
} \
 \
template <class VALUE, class DERIVED, class OTHER, class OTHER_DERIVED, int ... N> \
inline DERIVED &  \
operator op##=(TinyArrayBase<VALUE, DERIVED, N...> & l, \
               TinyArrayBase<OTHER, OTHER_DERIVED, N...> const & r) \
{ \
    vigra_tinyarray_check_size(l.size() == r.size(), \
        "TinyArrayBase::operator" #op "=(): size mismatch."); \
    for(int i=0; i<l.size(); ++i) \
        l[i] op##= r[i]; \
    return static_cast<DERIVED &>(l); \
} \
template <class V1, class D1, class V2, class D2, int ... N> \
inline \
TinyArray<PromoteType<V1, V2>, N...> \
operator op(TinyArrayBase<V1, D1, N...> const & l, \
             TinyArrayBase<V2, D2, N...> const & r) \
{ \
    TinyArray<PromoteType<V1, V2>, N...> res(l); \
    return res op##= r; \
} \
 \
template <class V1, class D1, class V2, int ... N> \
inline \
TinyArray<PromoteTypeIf<std::is_arithmetic<V2>::value, V1, V2>, N...> \
operator op(TinyArrayBase<V1, D1, N...> const & l, \
             V2 r) \
{ \
    TinyArray<PromoteType<V1, V2>, N...> res(l); \
    return res op##= r; \
} \
 \
template <class V1, class V2, class D2, int ... N> \
inline \
TinyArray<PromoteTypeIf<std::is_arithmetic<V1>::value, V1, V2>, N...> \
operator op(V1 l, \
            TinyArrayBase<V2, D2, N...> const & r) \
{ \
    TinyArray<PromoteType<V1, V2>, N...> res(l); \
    return res op##= r; \
} \
 \
template <class V1, class V2, class D2> \
inline \
TinyArray<PromoteTypeIf<std::is_arithmetic<V1>::value, V1, V2>, runtime_size> \
operator op(V1 l, \
            TinyArrayBase<V2, D2, runtime_size> const & r) \
{ \
    TinyArray<PromoteType<V1, V2>, runtime_size> res(r.size(), l); \
    return res op##= r; \
}

VIGRA_TINYARRAY_OPERATORS(+)
VIGRA_TINYARRAY_OPERATORS(-)
VIGRA_TINYARRAY_OPERATORS(*)
VIGRA_TINYARRAY_OPERATORS(/)
VIGRA_TINYARRAY_OPERATORS(%)

#undef VIGRA_TINYARRAY_OPERATORS

    /// Unary negation
template <class V, class D, int ... N>
inline
TinyArray<V, N...>
operator-(TinyArrayBase<V, D, N...> const & v)
{
    TinyArray<V, N...> res(v.size(), DontInit);
    for(int k=0; k < v.size(); ++k)
        res[k] = -v[k];
    return res;
}

using std::abs;

    /// element-wise absolute value
template <class V, class D, int ... N>
inline
TinyArray<V, N...>
abs(TinyArrayBase<V, D, N...> const & v)
{
    TinyArray<V, N...> res(v.size(), DontInit);
    for(int k=0; k < v.size(); ++k)
        res[k] = abs(v[k]);
    return res;
}

using std::ceil;

    /** Apply ceil() function to each vector component.
    */
template <class V, class D, int ... N>
inline
TinyArray<V, N...>
ceil(TinyArrayBase<V, D, N...> const & v)
{
    TinyArray<V, N...> res(v.size(), DontInit);
    for(int k=0; k < v.size(); ++k)
        res[k] = ceil(v[k]);
    return res;
}

using std::floor;

    /** Apply floor() function to each vector component.
    */
template <class V, class D, int ... N>
inline
TinyArray<V, N...>
floor(TinyArrayBase<V, D, N...> const & v)
{
    TinyArray<V, N...> res(v.size(), DontInit);
    for(int k=0; k < v.size(); ++k)
        res[k] = floor(v[k]);
    return res;
}

using std::round;

    /** Apply round() function to each vector component.
    */
template <class V, class D, int ... N>
inline
TinyArray<V, N...>
round(TinyArrayBase<V, D, N...> const & v)
{
    TinyArray<V, N...> res(v.size(), DontInit);
    for(int k=0; k < v.size(); ++k)
        res[k] = round(v[k]);
    return res;
}

    /** Apply roundi() function to each vector component, i.e. return an integer vector.
    */
template <class V, class D, int ... N>
inline
TinyArray<ArrayIndex, N...>
roundi(TinyArrayBase<V, D, N...> const & v)
{
    TinyArray<ArrayIndex, N...> res(v.size(), DontInit);
    for(int k=0; k < v.size(); ++k)
        res[k] = roundi(v[k]);
    return res;
}

using std::sqrt;

    /** Apply sqrt() function to each vector component.
    */
template <class V, class D, int ... N>
inline
TinyArray<RealPromoteType<V>, N...>
sqrt(TinyArrayBase<V, D, N...> const & v)
{
    TinyArray<RealPromoteType<V>, N...> res(v.size(), DontInit);
    for(int k=0; k < v.size(); ++k)
        res[k] = sqrt(v[k]);
    return res;
}

using std::pow;

    /** Apply pow() function to each vector component.
    */
template <class V, class D, class E, int ... N>
inline
TinyArray<PromoteType<V, E>, N...>
pow(TinyArrayBase<V, D, N...> const & v, E exponent)
{
    TinyArray<PromoteType<V, E>, N...> res(v.size(), DontInit);
    for(int k=0; k < v.size(); ++k)
        res[k] = pow(v[k], exponent);
    return res;
}

    /// cross product
template <class V1, class D1, class V2, class D2>
inline
TinyArray<PromoteType<V1, V2>, 3>
cross(TinyArrayBase<V1, D1, 3> const & r1,
      TinyArrayBase<V2, D2, 3> const & r2)
{
    typedef TinyArray<PromoteType<V1, V2>, 3> Res;
    return  Res(r1[1]*r2[2] - r1[2]*r2[1],
                r1[2]*r2[0] - r1[0]*r2[2],
                r1[0]*r2[1] - r1[1]*r2[0]);
}

template <class V1, class D1, class V2, class D2>
inline
TinyArray<PromoteType<V1, V2>, runtime_size>
cross(TinyArrayBase<V1, D1, runtime_size> const & r1,
      TinyArrayBase<V2, D2, runtime_size> const & r2)
{
    vigra_precondition(r1.size() == 3 && r2.size() == 3,
        "TinyArrayBase::cross(): cross product requires size() == 3.");
    typedef TinyArray<PromoteType<V1, V2>, runtime_size> Res;
    return  Res({r1[1]*r2[2] - r1[2]*r2[1],
                 r1[2]*r2[0] - r1[0]*r2[2],
                 r1[0]*r2[1] - r1[1]*r2[0]});
}

    /// dot product of two vectors
template <class V1, class D1, class V2, class D2, int N>
inline
PromoteType<V1, V2>
dot(TinyArrayBase<V1, D1, N> const & l,
    TinyArrayBase<V2, D2, N> const & r)
{
    if(N == runtime_size)
        vigra_precondition(l.size() == r.size(),
            "TinyArrayBase::dot(): size mismatch.");
    PromoteType<V1, V2> res = l[0] * r[0];
    for(int k=1; k < l.size(); ++k)
        res += l[k] * r[k];
    return res;
}

    /// sum of the vector's elements
template <class V, class D, int ... N>
inline
PromoteType<V>
sum(TinyArrayBase<V, D, N...> const & l)
{
    PromoteType<V> res = l[0];
    for(int k=1; k < l.size(); ++k)
        res += l[k];
    return res;
}

    /// mean of the vector's elements
template <class V, class D, int ... N>
inline RealPromoteType<V>
mean(TinyArrayBase<V, D, N...> const & t)
{
    const V sumVal = sum(t);
    return static_cast<RealPromoteType<V> >(sumVal) / t.size();
}

    /// cumulative sum of the vector's elements
template <class V, class D, int ... N>
inline
TinyArray<PromoteType<V>, N...>
cumsum(TinyArrayBase<V, D, N...> const & l)
{
    TinyArray<PromoteType<V>, N...> res(l);
    for(int k=1; k < l.size(); ++k)
        res[k] += res[k-1];
    return res;
}

    /// product of the vector's elements
template <class V, class D, int ... N>
inline
PromoteType<V>
prod(TinyArrayBase<V, D, N...> const & l)
{
    PromoteType<V> res = l[0];
    for(int k=1; k < l.size(); ++k)
        res *= l[k];
    return res;
}

    /// cumulative product of the vector's elements
template <class V, class D, int ... N>
inline
TinyArray<PromoteType<V>, N...>
cumprod(TinyArrayBase<V, D, N...> const & l)
{
    TinyArray<PromoteType<V>, N...> res(l);
    for(int k=1; k < l.size(); ++k)
        res[k] *= res[k-1];
    return res;
}

    /// \brief compute the C-order default stride of a given shape.
    /// Example: {200, 100, 50}  =>  {5000, 50, 1}
template <class V, class D, int N>
inline
TinyArray<PromoteType<V>, N>
shapeToStride(TinyArrayBase<V, D, N> const & shape)
{
    TinyArray<PromoteType<V>, N> res(shape.size(), DontInit);
    res[shape.size()-1] = 1;
    for(int k=shape.size()-2; k >= 0; --k)
        res[k] = res[k+1] * shape[k+1];
    return res;
}

using std::min;

template <class V1, class D1, class V2, class D2, int ... N>
inline
TinyArray<PromoteType<V1, V2>, N...>
minImpl(TinyArrayBase<V1, D1, N...> const & l,
        TinyArrayBase<V2, D2, N...> const & r)
{
    vigra_tinyarray_check_size(l.size() == r.size(),
        "TinyArrayBase::min(): size mismatch.");
    TinyArray<PromoteType<V1, V2>, N...> res(l.size(), DontInit);
    for(int k=0; k < l.size(); ++k)
        res[k] =  min<PromoteType<V1, V2> >(l[k], r[k]);
    return res;
}

    /// element-wise minimum
template <class V1, class D1, class V2, class D2, int ... N>
inline
TinyArray<PromoteType<V1, V2>, N...>
min(TinyArrayBase<V1, D1, N...> const & l,
    TinyArrayBase<V2, D2, N...> const & r)
{
    return minImpl(l, r);
}

// we also have to overload min for like-typed argument to prevent match of std::min()
template <class V, class D, int ...N>
inline TinyArray<V, N...>
min(TinyArrayBase<V, D, N...> const & l,
    TinyArrayBase<V, D, N...> const & r)
{
    return minImpl(l, r);
}

template <class V, int ...N>
inline TinyArray<V, N...>
min(TinyArray<V, N...> const & l,
    TinyArray<V, N...> const & r)
{
    return minImpl(l, r);
}

template <class V, int ...N>
inline TinyArray<V, N...>
min(TinyArrayView<V, N...> const & l,
    TinyArrayView<V, N...> const & r)
{
    return minImpl(l, r);
}

template <class V, int N>
inline TinyArray<V, N*(N+1)/2>
min(TinySymmetricView<V, N> const & l,
    TinySymmetricView<V, N> const & r)
{
    return minImpl(l, r);
}

    /// minimal element
template <class V, class D, int ... N>
inline
V const &
min(TinyArrayBase<V, D, N...> const & l)
{
    int m = 0;
    for(int i=1; i<l.size(); ++i)
        if(l[i] < l[m])
            m = i;
    return l[m];
}

using std::max;

template <class V1, class D1, class V2, class D2, int ... N>
inline
TinyArray<PromoteType<V1, V2>, N...>
maxImpl(TinyArrayBase<V1, D1, N...> const & l,
        TinyArrayBase<V2, D2, N...> const & r)
{
    vigra_tinyarray_check_size(l.size() == r.size(),
        "TinyArrayBase::max(): size mismatch.");
    TinyArray<PromoteType<V1, V2>, N...> res(l.size(), DontInit);
    for(int k=0; k < l.size(); ++k)
        res[k] =  max<PromoteType<V1, V2> >(l[k], r[k]);
    return res;
}

    /// element-wise maximum
template <class V1, class D1, class V2, class D2, int ... N>
inline
TinyArray<PromoteType<V1, V2>, N...>
max(TinyArrayBase<V1, D1, N...> const & l,
    TinyArrayBase<V2, D2, N...> const & r)
{
    return maxImpl(l, r);
}

// we also have to overload max for like-typed argument to prevent match of std::max()
template <class V, class D, int ...N>
inline TinyArray<V, N...>
max(TinyArrayBase<V, D, N...> const & l,
    TinyArrayBase<V, D, N...> const & r)
{
    return maxImpl(l, r);
}

template <class V, int ...N>
inline TinyArray<V, N...>
max(TinyArray<V, N...> const & l,
    TinyArray<V, N...> const & r)
{
    return maxImpl(l, r);
}

template <class V, int ...N>
inline TinyArray<V, N...>
max(TinyArrayView<V, N...> const & l,
    TinyArrayView<V, N...> const & r)
{
    return maxImpl(l, r);
}

template <class V, int N>
inline TinyArray<V, N*(N+1)/2>
max(TinySymmetricView<V, N> const & l,
    TinySymmetricView<V, N> const & r)
{
    return maxImpl(l, r);
}

    /// maximal element
template <class V, class D, int ... N>
inline
V const &
max(TinyArrayBase<V, D, N...> const & l)
{
    int m = 0;
    for(int i=1; i<l.size(); ++i)
        if(l[m] < l[i])
            m = i;
    return l[m];
}

    /// squared norm
template <class V, class D, int ... N>
inline
SquaredNormType<V>
squaredNorm(TinyArrayBase<V, D, N...> const & t)
{
    return t.squaredNormImpl();
}

template <class V, class D, int ... N>
inline 
NormType<V>
sizeDividedSquaredNorm(TinyArrayBase<V, D, N...> const & t)
{
    return NormType<V>(squaredNorm(t)) / t.size();
}

template <class V, class D, int ... N>
inline 
NormType<V>
sizeDividedNorm(TinyArrayBase<V, D, N...> const & t)
{
    return NormType<V>(norm(t)) / t.size();
}

    /// reversed copy
template <class V, class D, int ... N>
inline
TinyArray<V, N...>
reversed(TinyArrayBase<V, D, N...> const & t)
{
    return TinyArray<V, N...>(t.begin(), t.end(), ReverseCopy);
}

    /** \brief transposed copy
    
        Elements are arranged such that <tt>res[k] = t[permutation[k]]</tt>.
    */
template <class V1, class D1, class V2, class D2, int N>
inline
TinyArray<V1, N>
transpose(TinyArrayBase<V1, D1, N> const & v, 
          TinyArrayBase<V2, D2, N> const & permutation)
{
    TinyArray<V1, N> res(DontInit);
    for(int k=0; k < N; ++k)
    {
        VIGRA_ASSERT_INSIDE(permutation[k]);
        res[k] = v[permutation[k]];
    }
    return res;
}

template <class V1, class D1, int N>
inline
TinyArray<V1, N>
transpose(TinyArrayBase<V1, D1, N> const & v)
{
    return reversed(v);
}

template <class V1, class D1, int N1, int N2>
inline
TinyArray<V1, N2, N1>
transpose(TinyArrayBase<V1, D1, N1, N2> const & v)
{
    TinyArray<V1, N2, N1> res(DontInit);
    for(int i=0; i < N1; ++i)
    {
        for(int j=0; j < N2; ++j)
        {
            res(j,i) = v(i,j);
        }
    }
    return res;
}

template <class V, int N>
inline
TinySymmetricView<V, N>
transpose(TinySymmetricView<V, N> const & v)
{
    return v;
}

    /** \brief Clip negative values.
    
        All elements smaller than 0 are set to zero.
    */
template <class V, class D, int ... N>
inline
TinyArray<V, N...> 
clipLower(TinyArrayBase<V, D, N...> const & t)
{
    return clipLower(t, V());
}

    /** \brief Clip values below a threshold.
    
        All elements smaller than \a val are set to \a val.
    */
template <class V, class D, int ... N>
inline
TinyArray<V, N...> 
clipLower(TinyArrayBase<V, D, N...> const & t, const V val)
{
    TinyArray<V, N...> res(t.size(), DontInit);
    for(int k=0; k < t.size(); ++k)
    {
        res[k] = t[k] < val ? val :  t[k];
    }
    return res;
}

    /** \brief Clip values above a threshold.
    
        All elements bigger than \a val are set to \a val.
    */
template <class V, class D, int ... N>
inline
TinyArray<V, N...> 
clipUpper(TinyArrayBase<V, D, N...> const & t, const V val)
{
    TinyArray<V, N...> res(t.size(), DontInit);
    for(int k=0; k < t.size(); ++k)
    {
        res[k] = t[k] > val ? val :  t[k];
    }
    return res;
}

    /** \brief Clip values to an interval.
    
        All elements less than \a valLower are set to \a valLower, all elements
        bigger than \a valUpper are set to \a valUpper.
    */
template <class V, class D, int ... N>
inline
TinyArray<V, N...> 
clip(TinyArrayBase<V, D, N...> const & t,
     const V valLower, const V valUpper)
{
    TinyArray<V, N...> res(t.size(), DontInit);
    for(int k=0; k < t.size(); ++k)
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
template <class V, class D1, class D2, class D3, int ... N>
inline
TinyArray<V, N...> 
clip(TinyArrayBase<V, D1, N...> const & t,
     TinyArrayBase<V, D2, N...> const & valLower, 
     TinyArrayBase<V, D3, N...> const & valUpper)
{
    vigra_tinyarray_check_size(t.size() == valLower.size() && t.size() == valUpper.size(),
        "TinyArrayBase::clip(): size mismatch.");
    TinyArray<V, N...> res(t.size(), DontInit);
    for(int k=0; k < t.size(); ++k)
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

////////////////////////////////////////////////////////////
// PromoteTraits specializations

template <class T, class D, int ...N>
struct PromoteTraits<TinyArrayBase<T, D, N...>, TinyArrayBase<T, D, N...> >
{
    typedef TinyArrayBase<PromoteType<T>, D, N...>       Promote;
    typedef decltype(sqrt(TinyArrayBase<T, D, N...>()))  RealPromote;
};

template <class T, int ...N>
struct PromoteTraits<TinyArray<T, N...>, TinyArray<T, N...> >
{
    typedef TinyArray<PromoteType<T>, N...>       Promote;
    typedef decltype(sqrt(TinyArray<T, N...>()))  RealPromote;
};

template <class T, int ...N>
struct PromoteTraits<TinyArrayView<T, N...>, TinyArrayView<T, N...> >
{
    typedef TinyArrayView<PromoteType<T>, N...>       Promote;
    typedef decltype(sqrt(TinyArrayView<T, N...>()))  RealPromote;
};

template <class T, int N>
struct PromoteTraits<TinySymmetricView<T, N>, TinySymmetricView<T, N> >
{
    typedef TinySymmetricView<PromoteType<T>, N>       Promote;
    typedef decltype(sqrt(TinySymmetricView<T, N>()))  RealPromote;
};

////////////////////////////////////////////////////////////
// NormTraits specializations

template <class T, class D, int ...N>
struct NormTraits<TinyArrayBase<T, D, N...> >
{
    typedef vigra::SquaredNormType<T>         SquaredNormType;
    typedef RealPromoteType<SquaredNormType>  NormType;
};

template <class T, int ...N>
struct NormTraits<TinyArray<T, N...> >
{
    typedef vigra::SquaredNormType<T>         SquaredNormType;
    typedef RealPromoteType<SquaredNormType>  NormType;
};

template <class T, int ...N>
struct NormTraits<TinyArrayView<T, N...> >
{
    typedef vigra::SquaredNormType<T>         SquaredNormType;
    typedef RealPromoteType<SquaredNormType>  NormType;
};

template <class T, int N>
struct NormTraits<TinySymmetricView<T, N> >
{
    typedef vigra::SquaredNormType<T>         SquaredNormType;
    typedef RealPromoteType<SquaredNormType>  NormType;
};

// mask cl.exe shortcomings [end]
#if defined(_MSC_VER)
#pragma warning( pop )
#endif

} // namespace vigra

namespace vigra1 {

template <class T, int SIZE>
using TinyVector = vigra::TinyArray<T, SIZE>;

template <class T, int SIZE>
using TinyVectorView = vigra::TinyArray<T, SIZE>;

} // namespace vigra1

#undef VIGRA_ASSERT_INSIDE


#endif // VIGRA_TINYARRAY_HXX
