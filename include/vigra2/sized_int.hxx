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

#ifndef VIGRA2_SIZED_INT_HXX
#define VIGRA2_SIZED_INT_HXX

#include <cstdint>

namespace vigra {

using std::int8_t;
using std::int16_t;
using std::int32_t;
using std::int64_t;
using std::uint8_t;
using std::uint16_t;
using std::uint32_t;
using std::uint64_t;
using std::intmax_t;
using std::uintmax_t;

} // namespace vigra

#if 0 // the old method, maybe still needed for MSVC

#include "metaprogramming.hxx"
#include <limits>

#if   SHRT_MAX  == 0x7FL
# define VIGRA_BITSOF_SHORT 8
#elif SHRT_MAX  == 0x7FFFL
# define VIGRA_BITSOF_SHORT 16
#elif SHRT_MAX  == 0x7FFFFFFFL
# define VIGRA_BITSOF_SHORT 32
#elif SHRT_MAX  > 0xFFFFFFFFL
# define VIGRA_BITSOF_SHORT 64
#else
# define VIGRA_BITSOF_SHORT -1
#endif

#if   INT_MAX  == 0x7FL
# define VIGRA_BITSOF_INT 8
#elif INT_MAX  == 0x7FFFL
# define VIGRA_BITSOF_INT 16
#elif INT_MAX  == 0x7FFFFFFFL
# define VIGRA_BITSOF_INT 32
#elif INT_MAX  > 0xFFFFFFFFL
# define VIGRA_BITSOF_INT 64
#else
# define VIGRA_BITSOF_INT -1
#endif

#if   LONG_MAX  == 0x7FL
# define VIGRA_BITSOF_LONG 8
#elif LONG_MAX  == 0x7FFFL
# define VIGRA_BITSOF_LONG 16
#elif LONG_MAX  == 0x7FFFFFFFL
# define VIGRA_BITSOF_LONG 32
#elif LONG_MAX  > 0xFFFFFFFFL
# define VIGRA_BITSOF_LONG 64
#else
# define VIGRA_BITSOF_LONG -1
#endif

#if   LLONG_MAX  == 0x7FL
# define VIGRA_BITSOF_LONG_LONG 8
#elif LLONG_MAX  == 0x7FFFL
# define VIGRA_BITSOF_LONG_LONG 16
#elif LLONG_MAX  == 0x7FFFFFFFL
# define VIGRA_BITSOF_LONG_LONG 32
#elif LLONG_MAX  > 0xFFFFFFFFL
# define VIGRA_BITSOF_LONG_LONG 64
#else
# define VIGRA_BITSOF_LONG_LONG -1
#endif

namespace vigra {

class Int_type_not_supported_on_this_platform {};

#ifndef NO_PARTIAL_TEMPLATE_SPECIALIZATION

namespace detail {

template<class T, class NEXT>
struct IntTypeList
{
    enum { size = sizeof(T)*8 };
    typedef T type;
    typedef NEXT next;
};

template<int SIZE, class LIST>
struct SelectIntegerType
{
    typedef typename
       IfBool<(SIZE == LIST::size),
           typename LIST::type,
           typename SelectIntegerType<SIZE, typename LIST::next>::type >::type
       type;
};

template<int SIZE>
struct SelectIntegerType<SIZE, Int_type_not_supported_on_this_platform>
{
    typedef Int_type_not_supported_on_this_platform type;
};

template<class LIST>
struct SelectBiggestIntegerType
{
    enum { cursize = static_cast<int>(LIST::size),
           nextsize = static_cast<int>(SelectBiggestIntegerType<typename LIST::next>::size),
           size = (cursize < nextsize) ? nextsize : cursize };
    typedef typename
       IfBool<(cursize < nextsize),
           typename SelectBiggestIntegerType<typename LIST::next>::type,
           typename LIST::type>::type
       type;
};

template<>
struct SelectBiggestIntegerType<Int_type_not_supported_on_this_platform>
{
    enum { size = 0 };
    typedef Int_type_not_supported_on_this_platform type;
};

typedef IntTypeList<signed char,
        IntTypeList<signed short,
        IntTypeList<signed int,
        IntTypeList<signed long,
        IntTypeList<signed long long,
        Int_type_not_supported_on_this_platform > > > > > SignedIntTypes;
typedef IntTypeList<unsigned char,
        IntTypeList<unsigned short,
        IntTypeList<unsigned int,
        IntTypeList<unsigned long,
        IntTypeList<unsigned long long,
        Int_type_not_supported_on_this_platform > > > > > UnsignedIntTypes;

} // namespace detail

/** \addtogroup FixedSizeInt Fixed Size Integer Types

    Since the C++ standard does only specify minimal sizes for the built-in
    integer types, one cannot rely on them to have a specific size. But
    pixel types with a specific size are often required in image processing,
    especially when reading or writing binary files. The VIGRA typedefs
    are guaranteed to have exactly the correct size. If the system
    does not provide a suitable type, the typedef will evaluate to
    <tt>Int_type_not_supported_on_this_platform</tt>.
*/
//@{

    /// 8-bit signed int
typedef detail::SelectIntegerType<8,  detail::SignedIntTypes>::type int8_t;
    /// 16-bit signed int
typedef detail::SelectIntegerType<16, detail::SignedIntTypes>::type int16_t;
    /// 32-bit signed int
typedef detail::SelectIntegerType<32, detail::SignedIntTypes>::type int32_t;
    /// 64-bit signed int
typedef detail::SelectIntegerType<64, detail::SignedIntTypes>::type int64_t;
    /// 8-bit unsigned int
typedef detail::SelectIntegerType<8,  detail::UnsignedIntTypes>::type uint8_t;
    /// 16-bit unsigned int
typedef detail::SelectIntegerType<16, detail::UnsignedIntTypes>::type uint16_t;
    /// 32-bit unsigned int
typedef detail::SelectIntegerType<32, detail::UnsignedIntTypes>::type uint32_t;
    /// 64-bit unsigned int
typedef detail::SelectIntegerType<64, detail::UnsignedIntTypes>::type uint64_t;

    /// the biggest signed integer type of the system
typedef detail::SelectBiggestIntegerType<detail::SignedIntTypes>::type   IntBiggest;
    /// the biggest unsigned integer type of the system
typedef detail::SelectBiggestIntegerType<detail::UnsignedIntTypes>::type UIntBiggest;

//@}

#else // NO_PARTIAL_TEMPLATE_SPECIALIZATION

typedef signed char    int8_t;
typedef signed short   int16_t;
typedef signed int     int32_t;
typedef Int_type_not_supported_on_this_platform int64_t;
typedef unsigned char  uint8_t;
typedef unsigned short uint16_t;
typedef unsigned int   uint32_t;
typedef Int_type_not_supported_on_this_platform uint64_t;

typedef int32_t  IntBiggest;
typedef uint32_t UIntBiggest;

#endif // NO_PARTIAL_TEMPLATE_SPECIALIZATION

} // namespace vigra

#endif // #if 0

#endif /* VIGRA2_SIZED_INT_HXX */
