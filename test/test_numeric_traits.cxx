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

#include <typeinfo>
#include <iostream>
#include <string>
#include <muli/unittest.hxx>
#include <muli/numeric_traits.hxx>

using namespace muli;

struct NumericTraitsTest
{
    NumericTraitsTest()
    {
    }
    
    void testPromote()
    {
        should((std::is_same<unsigned int, muli::Promote<unsigned int, short> >::value));
        should((std::is_same<int, muli::Promote<unsigned char> >::value));
        should((std::is_same<double, muli::RealPromote<unsigned char> >::value));
        should((std::is_same<float, muli::RealPromote<float> >::value));
        should((std::is_same<float, muli::RealPromote<float, int> >::value));
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
            
            should((std::is_same<int, T::promote_type>::value));
            should((std::is_same<unsigned int, T::unsigned_type>::value));
            should((std::is_same<double, T::real_promote_type>::value));
            should((std::is_same<std::complex<double>, T::complex_promote_type>::value));
            
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
            
            should((std::is_same<int, T::promote_type>::value));
            should((std::is_same<unsigned int, T::unsigned_type>::value));
            should((std::is_same<double, T::real_promote_type>::value));
            should((std::is_same<std::complex<double>, T::complex_promote_type>::value));
            
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
            
            should((std::is_same<float, T::promote_type>::value));
            should((std::is_same<float, T::unsigned_type>::value));
            should((std::is_same<float, T::real_promote_type>::value));
            should((std::is_same<std::complex<float>, T::complex_promote_type>::value));
            
            shouldEqual(T::fromPromote(-FLT_MAX), T::min());
            shouldEqual(T::fromPromote(FLT_MAX), T::max());
            shouldEqual(T::fromPromote(1.25f), 1.25f);
            shouldEqual(T::fromRealPromote(-DBL_MAX), T::min());
            shouldEqual(T::fromRealPromote(DBL_MAX), T::max());
            shouldEqual(T::fromRealPromote(1.25), 1.25f);
        }
    }
};

struct NumericTraitsTestSuite
: public muli::test_suite
{
    NumericTraitsTestSuite()
    : muli::test_suite("NumericTraitsTest")
    {
        add( testCase(&NumericTraitsTest::testPromote));
        add( testCase(&NumericTraitsTest::testNumericTraits));
    }
};

int main(int argc, char ** argv)
{
    NumericTraitsTestSuite test;

    int failed = test.run(muli::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}
