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

#include <typeinfo>
#include <iostream>
#include <string>
#include <vigra2/unittest.hxx>
#include <vigra2/numeric_traits.hxx>

using namespace vigra;

struct NumericTraitsTest
{
    NumericTraitsTest()
    {
    }
    
    void testPromote()
    {
        should((std::is_same<unsigned int, PromoteType<unsigned int, short> >::value));
        should((std::is_same<int, PromoteType<unsigned char> >::value));
        should((std::is_same<float, PromoteType<float> >::value));
        
        should((std::is_same<double, RealPromoteType<unsigned char> >::value));
        should((std::is_same<long double, RealPromoteType<long double> >::value));
        should((std::is_same<float, RealPromoteType<float> >::value));
        should((std::is_same<float, RealPromoteType<float, int> >::value));
        
        should((std::is_same<int, SquaredNormType<int> >::value));
        should((std::is_same<unsigned int, SquaredNormType<unsigned int> >::value));
        should((std::is_same<int, SquaredNormType<unsigned char> >::value));
        should((std::is_same<float, SquaredNormType<float> >::value));
        should((std::is_same<long double, SquaredNormType<long double> >::value));
        
        should((std::is_same<int, NormType<int> >::value));
        should((std::is_same<unsigned int, NormType<unsigned int> >::value));
        should((std::is_same<unsigned char, NormType<unsigned char> >::value));
        should((std::is_same<float, NormType<float> >::value));
        should((std::is_same<long double, NormType<long double> >::value));
    }
    
    void testNumericTraits()
    {
        {
            typedef NumericTraits<signed short> T;
            shouldEqual(0, T::zero());
            shouldEqual(1, T::one());
            shouldEqual(1, T::nonZero());
            shouldEqual(1, T::smallestPositive());
            shouldEqual(1, T::epsilon());
            shouldEqual(SHRT_MIN, T::min());
            shouldEqual(SHRT_MAX, T::max());
            
            should((std::is_same<int, T::Promote>::value));
            should((std::is_same<unsigned int, T::UnsignedPromote>::value));
            should((std::is_same<double, T::RealPromote>::value));
            should((std::is_same<std::complex<double>, T::ComplexPromote>::value));
            
            shouldEqual(T::fromPromote(INT_MIN), T::min());
            shouldEqual(T::fromPromote(INT_MAX), T::max());
            shouldEqual(T::fromPromote(1), 1);
            shouldEqual(T::fromRealPromote(-FLT_MAX), T::min());
            shouldEqual(T::fromRealPromote(FLT_MAX), T::max());
            shouldEqual(T::fromRealPromote(1.25), 1);
        }
        {
            typedef NumericTraits<unsigned short> T;
            shouldEqual(0, T::zero());
            shouldEqual(1, T::one());
            shouldEqual(1, T::nonZero());
            shouldEqual(1, T::smallestPositive());
            shouldEqual(1, T::epsilon());
            shouldEqual(0, T::min());
            shouldEqual(USHRT_MAX, T::max());
            
            should((std::is_same<int, T::Promote>::value));
            should((std::is_same<unsigned int, T::UnsignedPromote>::value));
            should((std::is_same<double, T::RealPromote>::value));
            should((std::is_same<std::complex<double>, T::ComplexPromote>::value));
            
            shouldEqual(T::fromPromote(INT_MIN), T::min());
            shouldEqual(T::fromPromote(INT_MAX), T::max());
            shouldEqual(T::fromPromote(1), 1);
            shouldEqual(T::fromRealPromote(-FLT_MAX), T::min());
            shouldEqual(T::fromRealPromote(FLT_MAX), T::max());
            shouldEqual(T::fromRealPromote(1.25), 1);
        }
        {
            typedef NumericTraits<float> T;
            shouldEqual(0.0f, T::zero());
            shouldEqual(1.0f, T::one());
            shouldEqual(1.0f, T::nonZero());
            shouldEqual(FLT_MIN, T::smallestPositive());
            shouldEqual(FLT_EPSILON, T::epsilon());
            shouldEqual(-FLT_MAX, T::min());
            shouldEqual(FLT_MAX, T::max());
            
            should((std::is_same<float, T::Promote>::value));
            should((std::is_same<float, T::UnsignedPromote>::value));
            should((std::is_same<float, T::RealPromote>::value));
            should((std::is_same<std::complex<float>, T::ComplexPromote>::value));
            
            shouldEqual(T::fromPromote(-FLT_MAX), T::min());
            shouldEqual(T::fromPromote(FLT_MAX), T::max());
            shouldEqual(T::fromPromote(1.25f), 1.25f);
            shouldEqual(T::fromRealPromote(-DBL_MAX), T::min());
            shouldEqual(T::fromRealPromote(DBL_MAX), T::max());
            shouldEqual(T::fromRealPromote(1.25), 1.25f);
        }
    }
    
    void testCloseAtTolerance()
    {
        typedef NumericTraits<float> T;
        
        should(closeAtTolerance(T::zero(), T::epsilon()));
        should(!closeAtTolerance(T::zero(), T::one()));
        should(closeAtTolerance(T::zero(), T::one(), 2.0f));
        
        should(lessEqualAtTolerance(T::zero(), T::epsilon()));
        should(lessEqualAtTolerance(T::epsilon(), T::zero()));
        should(lessEqualAtTolerance(T::zero(), T::one()));
        should(!lessEqualAtTolerance(T::one(), T::zero()));
        should(lessEqualAtTolerance(T::one(), T::zero(), 2.0f));
        
        should(greaterEqualAtTolerance(T::zero(), T::epsilon()));
        should(greaterEqualAtTolerance(T::epsilon(), T::zero()));
        should(!greaterEqualAtTolerance(T::zero(), T::one()));
        should(greaterEqualAtTolerance(T::zero(), T::one(), 2.0f));
        should(greaterEqualAtTolerance(T::one(), T::zero()));

        double a = 0.0, b = NumericTraits<double>::epsilon(), c = 1000.0, d = 1000.1;

        using namespace vigra;
        should(closeAtTolerance(a, b));
        should(closeAtTolerance(c, c + b));
        should(!closeAtTolerance(c, d));
        should(closeAtTolerance(c, d, 0.01));
        should(closeAtTolerance(-a, -b));
        should(closeAtTolerance(-c, -c + b));
        should(!closeAtTolerance(-c, -d));
        should(closeAtTolerance(-c, -d, 0.01));
        should(!closeAtTolerance(c, -c));
        should(!closeAtTolerance(a, c));
        should(!closeAtTolerance(c, d));
        should(!closeAtTolerance(-a, -c));
        should(!closeAtTolerance(-c, -d));

        should(lessEqualAtTolerance(a, c));
        should(lessEqualAtTolerance(c, c));
        should(lessEqualAtTolerance(c, c-b));
        should(lessEqualAtTolerance(c, d));
        should(!lessEqualAtTolerance(d, c));
        should(lessEqualAtTolerance(d, c, 0.01));

        should(greaterEqualAtTolerance(c, a));
        should(greaterEqualAtTolerance(c, c));
        should(greaterEqualAtTolerance(c-b, c));
        should(greaterEqualAtTolerance(d, c));
        should(!greaterEqualAtTolerance(c, d));
        should(greaterEqualAtTolerance(c, d, 0.01));
    }
};

struct NumericTraitsTestSuite
: public test_suite
{
    NumericTraitsTestSuite()
    : test_suite("NumericTraitsTest")
    {
        add( testCase(&NumericTraitsTest::testPromote));
        add( testCase(&NumericTraitsTest::testNumericTraits));
        add( testCase(&NumericTraitsTest::testCloseAtTolerance));
    }
};

int main(int argc, char ** argv)
{
    NumericTraitsTestSuite test;

    int failed = test.run(testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}
