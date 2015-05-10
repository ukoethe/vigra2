/************************************************************************/
/*                                                                      */
/*               Copyright 2014-2015 by Ullrich Koethe                  */
/*                                                                      */
/*    This file is part of the MULI computer vision library.            */
/*    The MULI Website is                                               */
/*        http://ukoethe.github.io/muli                                 */
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

#ifndef MULI_NUMERIC_TRAITS_HXX
#define MULI_NUMERIC_TRAITS_HXX

#include <type_traits>
#include <cmath>    // abs(double)
#include <cstdlib>  // abs(int)
#include <complex>
#include "config.hxx"

namespace muli {

namespace numeric_traits_detail {

using std::sqrt;
// using muli::sqrt;

}

template <class T1, class T2 = T1>
using Promote = decltype(T1()+T2());

template <bool Cond, class T1, class T2 = T1>
using PromoteIf = typename std::enable_if<Cond, Promote<T1, T2> >::type;

template <class T1, class T2 = T1>
using RealPromote = decltype(numeric_traits_detail::sqrt(Promote<T1, T2>()));

struct Error_NumericTraits_not_specialized_for_this_case { };
struct Error_NumericTraits_char_is_not_a_numeric_type__use_signed_char_or_unsigned_char { };

template<class A>
struct NumericTraits
{
    typedef Error_NumericTraits_not_specialized_for_this_case type;
    typedef Error_NumericTraits_not_specialized_for_this_case promote_type;
    typedef Error_NumericTraits_not_specialized_for_this_case unsigned_type;
    typedef Error_NumericTraits_not_specialized_for_this_case real_promote_type;
    typedef Error_NumericTraits_not_specialized_for_this_case complex_promote_type;
    typedef Error_NumericTraits_not_specialized_for_this_case value_type;
};

template<>
struct NumericTraits<char>
{
    typedef Error_NumericTraits_char_is_not_a_numeric_type__use_signed_char_or_unsigned_char type;
    typedef Error_NumericTraits_char_is_not_a_numeric_type__use_signed_char_or_unsigned_char promote_type;
    typedef Error_NumericTraits_char_is_not_a_numeric_type__use_signed_char_or_unsigned_char unsigned_type;
    typedef Error_NumericTraits_char_is_not_a_numeric_type__use_signed_char_or_unsigned_char real_promote_type;
    typedef Error_NumericTraits_char_is_not_a_numeric_type__use_signed_char_or_unsigned_char complex_promote_type;
    typedef Error_NumericTraits_char_is_not_a_numeric_type__use_signed_char_or_unsigned_char value_type;
};

template<>
struct NumericTraits<bool>
{
    typedef bool type;
    typedef int promote_type;
    typedef unsigned int unsigned_type;
    typedef double real_promote_type;
    typedef std::complex<real_promote_type> complex_promote_type;
    typedef type value_type;

    static constexpr bool zero() noexcept { return false; }
    static constexpr bool one() noexcept { return true; }
    static constexpr bool nonZero() noexcept { return true; }
    static constexpr type epsilon() noexcept { return true; }
    static constexpr type smallestPositive() noexcept { return true; }
    static constexpr bool min() noexcept { return false; }
    static constexpr bool max() noexcept { return true; }
    
    static const bool minConst = false;
    static const bool maxConst = true;
    
    static promote_type toPromote(bool v) { return v ? 1 : 0; }
    static real_promote_type toRealPromote(bool v) { return v ? 1.0 : 0.0; }
    static bool fromPromote(promote_type v) { 
        return (v == 0) ? false : true; 
    }
    static bool fromRealPromote(real_promote_type v) {
        return (v == 0.0) ? false : true; 
    }
};

template<class T>
struct SignedNumericTraits
{
    typedef T             type;
    typedef type          value_type;
    typedef Promote<type> promote_type;
    typedef typename std::make_unsigned<promote_type>::type unsigned_type;
    typedef RealPromote<type> real_promote_type;
    typedef std::complex<real_promote_type> complex_promote_type;
    
    static constexpr type zero() noexcept    { return 0; }
    static constexpr type one() noexcept     { return 1; }
    static constexpr type nonZero() noexcept { return 1; }
    static constexpr type epsilon() noexcept { return 1; }
    static constexpr type smallestPositive() noexcept { return 1; }
    static constexpr type min() noexcept     { return std::numeric_limits<T>::min(); }
    static constexpr type max() noexcept     { return std::numeric_limits<T>::max(); }
    
    static const type minConst = min();
    static const type maxConst = max();
    
    static promote_type      toPromote(type v)     { return v; }
    static real_promote_type toRealPromote(type v) { return v; }
    
    static type fromPromote(promote_type v) 
    {
        return v <= static_cast<promote_type>(min())
                   ? min()
                   : v >= static_cast<promote_type>(max()) 
                          ? max()
                          : static_cast<type>(v);
    }
    static type fromRealPromote(real_promote_type v) 
    {
        return v <= static_cast<real_promote_type>(min())
                   ? min()
                   : v >= static_cast<real_promote_type>(max())
                          ? max()
                          : static_cast<type>(std::round(v));
    }
};

template<>
struct NumericTraits<signed char> : public SignedNumericTraits<signed char> {};
template<>
struct NumericTraits<signed short> : public SignedNumericTraits<signed short> {};
template<>
struct NumericTraits<signed int> : public SignedNumericTraits<signed int> {};
template<>
struct NumericTraits<signed long> : public SignedNumericTraits<signed long> {};
template<>
struct NumericTraits<signed long long> : public SignedNumericTraits<signed long long> {};

template<class T>
struct UnsignedNumericTraits
{
    typedef T             type;
    typedef type          value_type;
    typedef Promote<type> promote_type;
    typedef typename std::make_unsigned<promote_type>::type unsigned_type;
    typedef RealPromote<type> real_promote_type;
    typedef std::complex<real_promote_type> complex_promote_type;
    
    static constexpr type zero() noexcept    { return 0; }
    static constexpr type one() noexcept     { return 1; }
    static constexpr type nonZero() noexcept { return 1; }
    static constexpr type epsilon() noexcept { return 1; }
    static constexpr type smallestPositive() noexcept { return 1; }
    static constexpr type min() noexcept     { return std::numeric_limits<T>::min(); }
    static constexpr type max() noexcept     { return std::numeric_limits<T>::max(); }
    
    static const type minConst = min();
    static const type maxConst = max();
    
    static promote_type      toPromote(type v)     { return v; }
    static real_promote_type toRealPromote(type v) { return v; }
    
    static type fromPromote(promote_type v) 
    {
        return v <= static_cast<promote_type>(zero())
                   ? zero()
                   : v >= static_cast<promote_type>(max()) 
                          ? max()
                          : static_cast<type>(v);
    }
    static type fromRealPromote(real_promote_type v) 
    {
        return v <= static_cast<real_promote_type>(zero())
                   ? zero()
                   : v >= static_cast<real_promote_type>(max())
                          ? max()
                          : static_cast<type>(std::round(v));
    }
};

template<>
struct NumericTraits<unsigned char> : public UnsignedNumericTraits<unsigned char> {};
template<>
struct NumericTraits<unsigned short> : public UnsignedNumericTraits<unsigned short> {};
template<>
struct NumericTraits<unsigned int> : public UnsignedNumericTraits<unsigned int> {};
template<>
struct NumericTraits<unsigned long> : public UnsignedNumericTraits<unsigned long> {};
template<>
struct NumericTraits<unsigned long long> : public UnsignedNumericTraits<unsigned long long> {};

template<class T>
struct FloatNumericTraits
{
    typedef T    type;
    typedef type value_type;    
    typedef type promote_type;
    typedef type unsigned_type;
    typedef type real_promote_type;
    typedef std::complex<real_promote_type> complex_promote_type;
    
    static constexpr type zero() noexcept { return 0.0; }
    static constexpr type one() noexcept { return 1.0; }
    static constexpr type nonZero() noexcept { return 1.0; }
    static constexpr type epsilon() noexcept { return std::numeric_limits<type>::epsilon(); }
    static constexpr type smallestPositive() noexcept { return std::numeric_limits<type>::min(); }
    static constexpr type min() noexcept { return std::numeric_limits<type>::lowest(); }
    static constexpr type max() noexcept { return std::numeric_limits<type>::max(); }
    
    static promote_type toPromote(type v) { return v; }
    static real_promote_type toRealPromote(type v) { return v; }
    static type fromPromote(promote_type v) { return v; }
    static type fromRealPromote(real_promote_type v) 
    {
        return v <= static_cast<real_promote_type>(min())
                   ? min()
                   : v >= static_cast<real_promote_type>(max())
                          ? max()
                          : static_cast<type>(v);
    }
};

template<>
struct NumericTraits<float> : public FloatNumericTraits<float> {};
template<>
struct NumericTraits<double> : public FloatNumericTraits<double> {};
template<>
struct NumericTraits<long double> : public FloatNumericTraits<long double> {};


template<class T>
struct NumericTraits<std::complex<T> >
{
    typedef std::complex<T> type;
    typedef std::complex<typename NumericTraits<T>::promote_type> promote_type;
    typedef std::complex<typename NumericTraits<T>::unsigned_type> unsigned_type;
    typedef std::complex<typename NumericTraits<T>::real_promote_type> real_promote_type;
    typedef std::complex<real_promote_type> complex_promote_type;
    typedef T value_type;

    static type zero() { return type(0.0); }
    static type one() { return type(1.0); }
    static type nonZero() { return one(); }
    static type epsilon() { return type(NumericTraits<T>::epsilon()); }
    static type smallestPositive() { return type(NumericTraits<T>::smallestPositive()); }

    static promote_type toPromote(type const & v) { return v; }
    static type fromPromote(promote_type const & v) { return v; }
    static type fromRealPromote(real_promote_type v) { return type(v); }
};

namespace numeric_traits_detail {

// both f1 and f2 are unsigned here
template<class FPT>
inline
FPT safeFloatDivision( FPT f1, FPT f2 )
{
    return  f2 < NumericTraits<FPT>::one() && f1 > f2 * NumericTraits<FPT>::max()
                ? NumericTraits<FPT>::max() 
                : (f2 > NumericTraits<FPT>::one() && f1 < f2 * NumericTraits<FPT>::smallestPositive()) || 
                   f1 == NumericTraits<FPT>::zero()
                     ? NumericTraits<FPT>::zero() 
                     : f1/f2;
}

} // namespace numeric_traits_detail
    
    /** \brief Tolerance based floating-point equality.

        Check whether two floating point numbers are equal within the given tolerance.
        This is useful because floating point numbers that should be equal in theory are
        rarely exactly equal in practice. If the tolerance \a epsilon is not given,
        twice the machine epsilon is used.

        <b>\#include</b> \<vigra/mathutil.hxx\><br>
        Namespace: vigra
    */
template <class T1, class T2>
typename std::enable_if<std::is_floating_point<Promote<T1, T2> >::value,
                        bool>::type
closeAtTolerance(T1 l, T2 r, Promote<T1, T2> epsilon)
{
    using std::abs;
    typedef Promote<T1, T2> T;
    if(l == 0.0)
        return abs(r) <= epsilon;
    if(r == 0.0)
        return abs(l) <= epsilon;
    T diff = abs( l - r );
    T d1   = numeric_traits_detail::safeFloatDivision<T>( diff, abs( r ) );
    T d2   = numeric_traits_detail::safeFloatDivision<T>( diff, abs( l ) );

    return (d1 <= epsilon && d2 <= epsilon);
}

template <class T1, class T2>
inline
typename std::enable_if<std::is_floating_point<Promote<T1, T2> >::value,
                        bool>::type 
closeAtTolerance(T1 l, T2 r)
{
    typedef Promote<T1, T2> T;
    return closeAtTolerance(l, r, T(2.0) * NumericTraits<T>::epsilon());
}
    
    /** \brief Tolerance based floating-point less-or-equal.

        Check whether two floating point numbers are less or equal within the given tolerance.
        That is, \a l can actually be greater than \a r within the given \a epsilon.
        This is useful because floating point numbers that should be equal in theory are
        rarely exactly equal in practice. If the tolerance \a epsilon is not given,
        twice the machine epsilon is used.

        <b>\#include</b> \<vigra/mathutil.hxx\><br>
        Namespace: vigra
    */
template <class T1, class T2>
inline
typename std::enable_if<std::is_floating_point<Promote<T1, T2> >::value,
                        bool>::type 
lessEqualAtTolerance(T1 l, T2 r, Promote<T1, T2> epsilon)
{
    return l < r || closeAtTolerance(l, r, epsilon);
}

template <class T1, class T2>
inline
typename std::enable_if<std::is_floating_point<Promote<T1, T2> >::value,
                        bool>::type 
lessEqualAtTolerance(T1 l, T2 r)
{
    typedef Promote<T1, T2> T;
    return lessEqualAtTolerance(l, r, T(2.0) * NumericTraits<T>::epsilon());
}
    
    /** \brief Tolerance based floating-point greater-or-equal.

        Check whether two floating point numbers are greater or equal within the given tolerance.
        That is, \a l can actually be less than \a r within the given \a epsilon.
        This is useful because floating point numbers that should be equal in theory are
        rarely exactly equal in practice. If the tolerance \a epsilon is not given,
        twice the machine epsilon is used.

        <b>\#include</b> \<vigra/mathutil.hxx\><br>
        Namespace: vigra
    */
template <class T1, class T2>
inline
typename std::enable_if<std::is_floating_point<Promote<T1, T2> >::value,
                        bool>::type 
greaterEqualAtTolerance(T1 l, T2 r, Promote<T1, T2> epsilon)
{
    return r < l || closeAtTolerance(l, r, epsilon);
}

template <class T1, class T2>
inline bool greaterEqualAtTolerance(T1 l, T2 r)
{
    typedef Promote<T1, T2> T;
    return greaterEqualAtTolerance(l, r, T(2.0) * NumericTraits<T>::epsilon());
}

} // namespace muli

#endif // MULI_NUMERIC_TRAITS_HXX
