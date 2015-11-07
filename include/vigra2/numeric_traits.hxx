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

#ifndef VIGRA2_NUMERIC_TRAITS_HXX_HXX
#define VIGRA2_NUMERIC_TRAITS_HXX_HXX

#include <type_traits>
#include <cmath>    // abs(double)
#include <cstdlib>  // abs(int)
#include <complex>
#include "config.hxx"

namespace vigra {

/////////////////////////////////////////////////////////
// PromoteTraits

template <class T1, class T2=T1>
struct PromoteTraits
{
    typedef decltype(*(T1*)0 + *(T2*)0)                   Promote;
    typedef typename PromoteTraits<Promote>::RealPromote  RealPromote;
};

template <class T>
struct PromoteTraits<T, T>
{
    typedef decltype(*(T*)0 + *(T*)0)  Promote;
    typedef decltype(sqrt(*(T*)0))     RealPromote;
};

template <class T>
struct PromoteTraits<T *, T *>
{
    typedef T *  Promote;
    typedef T *  RealPromote;
};

template <class T1, class T2=T1>
using PromoteType = typename PromoteTraits<T1, T2>::Promote;

template <bool Cond, class T1, class T2 = T1>
using PromoteTypeIf = typename std::enable_if<Cond, PromoteType<T1, T2> >::type;

template <class T1, class T2=T1>
using RealPromoteType = typename PromoteTraits<T1, T2>::RealPromote;

///////////////////////////////////////////////////////////////
// NumericTraits

template<class T>
struct NumericTraits
{
    typedef T                                          Type;
    typedef typename PromoteTraits<T>::Promote         Promote;
    typedef typename PromoteTraits<T>::RealPromote     RealPromote;
    typedef std::complex<RealPromote>                  ComplexPromote;
    typedef T                                          ValueType;
};

///////////////////////////////////////////////////////////////
// NormTraits

template<class T>
struct NormTraits
{
    typedef PromoteType<T>                    SquaredNormType;
    typedef RealPromoteType<SquaredNormType>  NormType;
};

template <class T>
using SquaredNormType = typename NormTraits<T>::SquaredNormType;

template <class T>
using NormType = typename NormTraits<T>::NormType;

////////////////////////////////////////////////////////////////
// PromoteTraits specializations

template <>
struct PromoteTraits<float, float>
{
    typedef float Promote;
    typedef float RealPromote;
};

template <>
struct PromoteTraits<long double, long double>
{
    typedef long double Promote;
    typedef long double RealPromote;
};


//////////////////////////////////////////////////////
// NumericTraits specializations

struct Error_NumericTraits_char_is_not_a_numeric_type__use_signed_char_or_unsigned_char { };

template<>
struct NumericTraits<char>
{
    typedef Error_NumericTraits_char_is_not_a_numeric_type__use_signed_char_or_unsigned_char Type;
    typedef Error_NumericTraits_char_is_not_a_numeric_type__use_signed_char_or_unsigned_char Promote;
    typedef Error_NumericTraits_char_is_not_a_numeric_type__use_signed_char_or_unsigned_char UnsignedPromote;
    typedef Error_NumericTraits_char_is_not_a_numeric_type__use_signed_char_or_unsigned_char RealPromote;
    typedef Error_NumericTraits_char_is_not_a_numeric_type__use_signed_char_or_unsigned_char ComplexPromote;
    typedef Error_NumericTraits_char_is_not_a_numeric_type__use_signed_char_or_unsigned_char ValueType;
};

template<>
struct NumericTraits<bool>
{
    typedef bool Type;
    typedef int Promote;
    typedef unsigned int UnsignedPromote;
    typedef double RealPromote;
    typedef std::complex<RealPromote> ComplexPromote;
    typedef Type ValueType;

    static constexpr bool zero() noexcept { return false; }
    static constexpr bool one() noexcept { return true; }
    static constexpr bool nonZero() noexcept { return true; }
    static constexpr Type epsilon() noexcept { return true; }
    static constexpr Type smallestPositive() noexcept { return true; }
    static constexpr bool min() noexcept { return false; }
    static constexpr bool max() noexcept { return true; }

    static const bool minConst = false;
    static const bool maxConst = true;

    static Promote toPromote(bool v) { return v ? 1 : 0; }
    static RealPromote toRealPromote(bool v) { return v ? 1.0 : 0.0; }
    static bool fromPromote(Promote v) {
        return (v == 0) ? false : true;
    }
    static bool fromRealPromote(RealPromote v) {
        return (v == 0.0) ? false : true;
    }
};

template<class T>
struct SignedNumericTraits
{
    typedef T             Type;
    typedef Type          ValueType;
    typedef PromoteType<Type> Promote;
    typedef typename std::make_unsigned<Promote>::type UnsignedPromote;
    typedef RealPromoteType<Type> RealPromote;
    typedef std::complex<RealPromote> ComplexPromote;

    static constexpr Type zero() noexcept    { return 0; }
    static constexpr Type one() noexcept     { return 1; }
    static constexpr Type nonZero() noexcept { return 1; }
    static constexpr Type epsilon() noexcept { return 1; }
    static constexpr Type smallestPositive() noexcept { return 1; }
    static constexpr Type min() noexcept     { return std::numeric_limits<T>::min(); }
    static constexpr Type max() noexcept     { return std::numeric_limits<T>::max(); }

    static const Type minConst = min();
    static const Type maxConst = max();

    static Promote      toPromote(Type v)     { return v; }
    static RealPromote  toRealPromote(Type v) { return v; }

    static Type fromPromote(Promote v)
    {
        return v <= static_cast<Promote>(min())
                   ? min()
                   : v >= static_cast<Promote>(max())
                          ? max()
                          : static_cast<Type>(v);
    }
    static Type fromRealPromote(RealPromote v)
    {
        return v <= static_cast<RealPromote>(min())
                   ? min()
                   : v >= static_cast<RealPromote>(max())
                          ? max()
                          : static_cast<Type>(std::round(v));
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
    typedef T             Type;
    typedef Type          ValueType;
    typedef PromoteType<Type> Promote;
    typedef typename std::make_unsigned<Promote>::type UnsignedPromote;
    typedef RealPromoteType<Type> RealPromote;
    typedef std::complex<RealPromote> ComplexPromote;

    static constexpr Type zero() noexcept    { return 0; }
    static constexpr Type one() noexcept     { return 1; }
    static constexpr Type nonZero() noexcept { return 1; }
    static constexpr Type epsilon() noexcept { return 1; }
    static constexpr Type smallestPositive() noexcept { return 1; }
    static constexpr Type min() noexcept     { return std::numeric_limits<T>::min(); }
    static constexpr Type max() noexcept     { return std::numeric_limits<T>::max(); }

    static const Type minConst = min();
    static const Type maxConst = max();

    static Promote      toPromote(Type v)     { return v; }
    static RealPromote  toRealPromote(Type v) { return v; }

    static Type fromPromote(Promote v)
    {
        return v <= static_cast<Promote>(zero())
                   ? zero()
                   : v >= static_cast<Promote>(max())
                          ? max()
                          : static_cast<Type>(v);
    }
    static Type fromRealPromote(RealPromote v)
    {
        return v <= static_cast<RealPromote>(zero())
                   ? zero()
                   : v >= static_cast<RealPromote>(max())
                          ? max()
                          : static_cast<Type>(std::round(v));
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
    typedef T    Type;
    typedef Type ValueType;
    typedef Type Promote;
    typedef Type UnsignedPromote;
    typedef Type RealPromote;
    typedef std::complex<RealPromote> ComplexPromote;

    static constexpr Type zero() noexcept { return 0.0; }
    static constexpr Type one() noexcept { return 1.0; }
    static constexpr Type nonZero() noexcept { return 1.0; }
    static constexpr Type epsilon() noexcept { return std::numeric_limits<Type>::epsilon(); }
    static constexpr Type smallestPositive() noexcept { return std::numeric_limits<Type>::min(); }
    static constexpr Type min() noexcept { return std::numeric_limits<Type>::lowest(); }
    static constexpr Type max() noexcept { return std::numeric_limits<Type>::max(); }

    static Promote      toPromote(Type v) { return v; }
    static RealPromote  toRealPromote(Type v) { return v; }
    static Type         fromPromote(Promote v) { return v; }
    static Type         fromRealPromote(RealPromote v)
    {
        return v <= static_cast<RealPromote>(min())
                   ? min()
                   : v >= static_cast<RealPromote>(max())
                          ? max()
                          : static_cast<Type>(v);
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
    typedef std::complex<T> Type;
    typedef std::complex<typename NumericTraits<T>::Promote> Promote;
    typedef std::complex<typename NumericTraits<T>::UnsignedPromote> UnsignedPromote;
    typedef std::complex<typename NumericTraits<T>::RealPromote> RealPromote;
    typedef std::complex<RealPromote> ComplexPromote;
    typedef T ValueType;

    static Type zero() { return Type(0.0); }
    static Type one() { return Type(1.0); }
    static Type nonZero() { return one(); }
    static Type epsilon() { return Type(NumericTraits<T>::epsilon()); }
    static Type smallestPositive() { return Type(NumericTraits<T>::smallestPositive()); }

    static Promote toPromote(Type const & v) { return v; }
    static Type    fromPromote(Promote const & v) { return v; }
    static Type    fromRealPromote(RealPromote v) { return Type(v); }
};

////////////////////////////////////////////////
// NormTraits specializations

template <class T>
struct FundamentalNormTraits
{
    typedef PromoteType<T>  SquaredNormType;
    typedef T               NormType;
};

template<>
struct NormTraits<signed char> : public FundamentalNormTraits<signed char> {};
template<>
struct NormTraits<signed short> : public FundamentalNormTraits<signed short> {};
template<>
struct NormTraits<signed int> : public FundamentalNormTraits<signed int> {};
template<>
struct NormTraits<signed long> : public FundamentalNormTraits<signed long> {};
template<>
struct NormTraits<signed long long> : public FundamentalNormTraits<signed long long> {};
template<>
struct NormTraits<unsigned char> : public FundamentalNormTraits<unsigned char> {};
template<>
struct NormTraits<unsigned short> : public FundamentalNormTraits<unsigned short> {};
template<>
struct NormTraits<unsigned int> : public FundamentalNormTraits<unsigned int> {};
template<>
struct NormTraits<unsigned long> : public FundamentalNormTraits<unsigned long> {};
template<>
struct NormTraits<unsigned long long> : public FundamentalNormTraits<unsigned long long> {};
template<>
struct NormTraits<float> : public FundamentalNormTraits<float> {};
template<>
struct NormTraits<double> : public FundamentalNormTraits<double> {};
template<>
struct NormTraits<long double> : public FundamentalNormTraits<long double> {};

///////////////////////////////////////////////////////////////
// RequiresExplicitCast

namespace detail {

template <class T>
struct RequiresExplicitCast {
    template <class U>
    static U const & cast(U const & v)
        { return v; }
};

#if !defined(_MSC_VER) || _MSC_VER >= 1300
#  define VIGRA_SPECIALIZED_CAST(type) \
    template <> \
    struct RequiresExplicitCast<type> { \
        static type cast(float v) \
            { return NumericTraits<type>::fromRealPromote(v); } \
        static type cast(double v) \
            { return NumericTraits<type>::fromRealPromote(v); } \
        static type cast(type v) \
            { return v; } \
        template <class U> \
        static type cast(U v) \
            { return static_cast<type>(v); } \
 \
    };
#else
#  define VIGRA_SPECIALIZED_CAST(type) \
    template <> \
    struct RequiresExplicitCast<type> { \
        static type cast(float v) \
            { return NumericTraits<type>::fromRealPromote(v); } \
        static type cast(double v) \
            { return NumericTraits<type>::fromRealPromote(v); } \
        static type cast(signed char v) \
            { return v; } \
        static type cast(unsigned char v) \
            { return v; } \
        static type cast(short v) \
            { return v; } \
        static type cast(unsigned short v) \
            { return v; } \
        static type cast(int v) \
            { return v; } \
        static type cast(unsigned int v) \
            { return v; } \
        static type cast(long v) \
            { return v; } \
        static type cast(unsigned long v) \
            { return v; } \
    };
#endif


VIGRA_SPECIALIZED_CAST(signed char)
VIGRA_SPECIALIZED_CAST(unsigned char)
VIGRA_SPECIALIZED_CAST(short)
VIGRA_SPECIALIZED_CAST(unsigned short)
VIGRA_SPECIALIZED_CAST(int)
VIGRA_SPECIALIZED_CAST(unsigned int)
VIGRA_SPECIALIZED_CAST(long)
VIGRA_SPECIALIZED_CAST(unsigned long)

template <>
struct RequiresExplicitCast<bool> {
    template <class U>
    static bool cast(U v)
    { return v == NumericTraits<U>::zero()
                ? false
                : true; }
};

template <>
struct RequiresExplicitCast<float> {
    static float cast(int v)
        { return (float)v; }

    static float cast(unsigned int v)
        { return (float)v; }

    static float cast(long v)
        { return (float)v; }

    static float cast(unsigned long v)
        { return (float)v; }

    static float cast(long long v)
        { return (float)v; }

    static float cast(unsigned long long v)
        { return (float)v; }

    static float cast(double v)
        { return (float)v; }

    static float cast(long double v)
        { return (float)v; }

    template <class U>
    static U cast(U v)
        { return v; }
};

template <>
struct RequiresExplicitCast<double> {
    static double cast(long long v)
        { return (double)v; }

    static double cast(unsigned long long v)
        { return (double)v; }

    template <class U>
    static U cast(U v)
        { return v; }
};

#undef VIGRA_SPECIALIZED_CAST

} // namespace detail

// comparison with tolerance

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
typename std::enable_if<std::is_floating_point<PromoteType<T1, T2> >::value,
                        bool>::type
closeAtTolerance(T1 l, T2 r, PromoteType<T1, T2> epsilon)
{
    using std::abs;
    typedef PromoteType<T1, T2> T;
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
typename std::enable_if<std::is_floating_point<PromoteType<T1, T2> >::value,
                        bool>::type
closeAtTolerance(T1 l, T2 r)
{
    typedef PromoteType<T1, T2> T;
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
typename std::enable_if<std::is_floating_point<PromoteType<T1, T2> >::value,
                        bool>::type
lessEqualAtTolerance(T1 l, T2 r, PromoteType<T1, T2> epsilon)
{
    return l < r || closeAtTolerance(l, r, epsilon);
}

template <class T1, class T2>
inline
typename std::enable_if<std::is_floating_point<PromoteType<T1, T2> >::value,
                        bool>::type
lessEqualAtTolerance(T1 l, T2 r)
{
    typedef PromoteType<T1, T2> T;
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
typename std::enable_if<std::is_floating_point<PromoteType<T1, T2> >::value,
                        bool>::type
greaterEqualAtTolerance(T1 l, T2 r, PromoteType<T1, T2> epsilon)
{
    return r < l || closeAtTolerance(l, r, epsilon);
}

template <class T1, class T2>
inline bool greaterEqualAtTolerance(T1 l, T2 r)
{
    typedef PromoteType<T1, T2> T;
    return greaterEqualAtTolerance(l, r, T(2.0) * NumericTraits<T>::epsilon());
}

} // namespace vigra

#endif // VIGRA2_NUMERIC_TRAITS_HXX_HXX
