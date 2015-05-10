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


#ifndef MULI_CONFIG_HXX
#define MULI_CONFIG_HXX

#include "config_version.hxx"
#include <stdexcept>

///////////////////////////////////////////////////////////
//                                                       //
//                     VisualC++                         //
//                                                       //
///////////////////////////////////////////////////////////

#ifdef _MSC_VER
    // make sure that we use muli/windows.h so that incompatibilities are fixed
    #include "windows.h"

    #if(_MSC_VER < 1100)    // before VisualC++ 5.0
        #error "Need VisualC++ 5.0, Service Pack 2, or later"
    #endif // _MSC_VER < 1100

    #if (_MSC_VER < 1300)
        #define NO_TYPENAME         // no 'typename' keyword
        #define TEMPLATE_COPY_CONSTRUCTOR_BUG
        #define NO_STL_MEMBER_TEMPLATES
        #define NO_INLINE_STATIC_CONST_DEFINITION
        #define CMATH_NOT_IN_STD
        #define NO_COVARIANT_RETURN_TYPES

        #ifdef MULI_NO_STD_MINMAX  // activate if necessary
        namespace std {

        template<class T>
        const T& min(const T& x, const T& y)
        {
            return (y < x)
                ? y
                : x;
        }

        template<class T>
        const T& max(const T& x, const T& y)
        {
            return (x < y)
                ? y
                : x;
        }
        }
        #endif // MULI_NO_STD_MINMAX
    #endif // (_MSC_VER < 1300)

    #if _MSC_VER < 1310
        #pragma warning( disable : 4786 4250 4244 4305)

        #define NO_PARTIAL_TEMPLATE_SPECIALIZATION
        #define NO_OUT_OF_LINE_MEMBER_TEMPLATES
        #include <cmath>

        #ifdef _MSC_EXTENSIONS
        #ifndef CMATH_NOT_IN_STD
                namespace std {
        #endif // CMATH_NOT_IN_STD
                inline double abs(double v) { return fabs(v); }
                inline float  abs(float v)  { return fabs(v); }
        #ifndef CMATH_NOT_IN_STD
                }
        #endif // CMATH_NOT_IN_STD
        #endif // _MSC_EXTENSIONS
    #endif // _MSC_VER < 1310

    #if _MSC_VER < 1400
        #define MULI_NO_WORKING_STRINGSTREAM
    #endif
    
    #if _MSC_VER < 1600
        #define MULI_NO_UNIQUE_PTR
    #endif
    
    #define MULI_NEED_BIN_STREAMS
    
    #define MULI_NO_THREADSAFE_STATIC_INIT  // at least up to _MSC_VER <= 1600, probably higher
    
    // usage: 
    //   static int * p = MULI_SAFE_STATIC(p, new int(42));
    //
    #define MULI_SAFE_STATIC(p, v) \
    0; while(p == 0) ::muli::detail::safeStaticInit(&p, v)
    
    namespace muli { namespace detail {
    template <class T>
    inline void safeStaticInit(T ** p, T * v)
    {
        if (InterlockedCompareExchangePointer((PVOID *)p, v, 0) != 0)
            delete v;
    }
    }} // namespace muli::detail
    
    #ifndef MULI_ENABLE_ANNOYING_WARNINGS
        #pragma warning ( disable: 4244 4267) // implicit integer conversion warnings
    #endif

    #ifdef MULI_DLL
        #define MULI_EXPORT __declspec(dllexport)
    #elif defined(MULI_STATIC_LIB)
        #define MULI_EXPORT
    #else
        #define MULI_EXPORT __declspec(dllimport)
    #endif
#endif // _MSC_VER

///////////////////////////////////////////////////////////
//                                                       //
//                           gcc                         //
//                                                       //
///////////////////////////////////////////////////////////

#if defined(__GNUC__) && !defined(__clang__)
    #if  __GNUC__ < 2 || ((__GNUC__ == 2) && (__GNUC_MINOR__ <= 8))
        #error "Need at least g++ 2.95"
    #endif
    #if __GNUC__ < 3
        #define MULI_NO_WORKING_STRINGSTREAM
    #endif
    #define HAS_HASH_CONTAINERS
    
    // these warnings produce too many false positives to be useful
    #pragma GCC diagnostic ignored "-Wshadow"  
    
    #if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
        #if defined(__APPLE__)
            #define MULI_NO_UNIQUE_PTR
        #endif
    #else
        // C++98 mode.  No native unique_ptr support
        #define MULI_NO_UNIQUE_PTR
        #define MULI_SHARED_PTR_IN_TR1
    #endif
#endif  // __GNUC__

///////////////////////////////////////////////////////////
//                                                       //
//                         clang                         //
//                                                       //
///////////////////////////////////////////////////////////
#if defined(__clang__)
    // In Apple builds of clang, __clang_major__ and __clang_minor__
    // have totally different values than in other builds of clang.
    #if defined(__apple_build_version__)
        // (For Apple builds of clang, __clang_major__ tracks the XCode version.)
        // For Apple builds, C++11 only works well with libc++, not stdlibc++
        #define MULI_NO_UNIQUE_PTR
        #if __cplusplus >= 201103L
            // Must have at least XCode 4 and use libc++ to use std::shared_ptr, etc.
            // Otherwise, use tr1.
            #if !((__clang_major__ >= 4) && defined(_LIBCPP_VERSION))
                #define MULI_SHARED_PTR_IN_TR1
            #endif
        #else
            // C++98 mode.  No native unique_ptr support
            #define MULI_NO_UNIQUE_PTR
            #if !defined(_LIBCPP_VERSION)
                #define MULI_SHARED_PTR_IN_TR1
            #endif
        #endif
    #else
        // This is a conservative constraint:
        // Full C++11 support was achieved in clang-3.3,
        // but most support was available in 3.1 and 3.2
        #if __cplusplus >= 201103L
            #if (__clang_major__ < 3) || ((__clang_major__ == 3) && (__clang_minor__ < 3))
                #define MULI_SHARED_PTR_IN_TR1
                #define MULI_NO_UNIQUE_PTR
            #endif
        #else
            // C++98 mode.  No native shared_ptr/unique_ptr support
            #define MULI_NO_UNIQUE_PTR
            #define MULI_SHARED_PTR_IN_TR1
        #endif
    #endif
#endif // __clang__

///////////////////////////////////////////////////////////
//                                                       //
//                         MingW                         //
//                                                       //
///////////////////////////////////////////////////////////

#if defined(__MINGW32__)
    #define MULI_NEED_BIN_STREAMS

    #ifdef MULI_DLL
        #define MULI_EXPORT __declspec(dllexport)
    #elif defined(MULI_STATIC_LIB)
        #define MULI_EXPORT
    #else
        #define MULI_EXPORT __declspec(dllimport)
    #endif
#endif  // __MINGW32__

///////////////////////////////////////////////////////////
//                                                       //
//                      SGI C++ 7.2                      //
//                                                       //
///////////////////////////////////////////////////////////

#if defined(__sgi) && !defined(__GNUC__)
    #if _COMPILER_VERSION < 720
        #error "Need SGI C++ 7.2 or later"
    #endif
    #if (_COMPILER_VERSION  == 720) || (_COMPILER_VERSION  == 721)
        #define SPECIAL_STDEXCEPTION_DEFINITION_NEEDED

        namespace muli {
            typedef std::exception StdException; // must be above next #define !!
        }
        #define std
        #define NO_NAMESPACE_STD
    #endif // _COMPILER_VERSION
    #define HAS_HASH_CONTAINERS
#endif // __sgi

///////////////////////////////////////////////////////////
//                                                       //
//                      Sun C++ ???                      //
//                                                       //
///////////////////////////////////////////////////////////

#if defined(__sun) && !defined(__GNUC__)
    #define MULI_HAS_ERF
#endif // __sun

///////////////////////////////////////////////////////////
//                                                       //
//                        general                        //
//                                                       //
///////////////////////////////////////////////////////////

#ifdef CMATH_NOT_IN_STD
    #define MULI_CSTD
#else
    #define MULI_CSTD std
#endif

#ifdef NO_TYPENAME
    #define typename
#endif

#ifdef NO_EXPLICIT
    #define explicit
#endif

#ifndef MULI_EXPORT
    #define MULI_EXPORT
#endif

#ifdef MULI_NO_UNIQUE_PTR
#  define MULI_UNIQUE_PTR  std::auto_ptr
#else
#  ifdef _GLIBCXX_INCLUDE_AS_TR1
#    define MULI_UNIQUE_PTR  std::tr1::unique_ptr
#  else
#    define MULI_UNIQUE_PTR  std::unique_ptr
#  endif
#endif

#ifdef MULI_SHARED_PTR_IN_TR1
#  define MULI_SHARED_PTR  std::tr1::shared_ptr
#else
#  define MULI_SHARED_PTR  std::shared_ptr
#endif

#ifndef MULI_NO_THREADSAFE_STATIC_INIT    
    // usage: 
    //   static int * p = MULI_SAFE_STATIC(p, new int(42));
    //
    #define MULI_SAFE_STATIC(p, v) v
#endif

namespace muli {

#ifndef SPECIAL_STDEXCEPTION_DEFINITION_NEEDED
     typedef std::exception StdException;
#endif

} // namespace muli

#ifdef DOXYGEN
#  define doxygen_overloaded_function(fun) fun(...);
#else
#  define doxygen_overloaded_function(fun)
#endif


#endif // MULI_CONFIG_HXX