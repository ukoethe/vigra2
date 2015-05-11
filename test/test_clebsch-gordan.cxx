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
#include <vigra2/clebsch-gordan.hxx>

using namespace vigra;

struct ClebschGordanTest
{
    ClebschGordanTest()
    {
    }
    
    void test()
    {
        shouldEqualTolerance(clebschGordan(0.5, 0.5, 0.5, 0.5, 1.0, 1.0), std::sqrt(1.0), 1e-15);
        shouldEqualTolerance(clebschGordan(0.5, 0.5, 0.5, -0.5, 1.0, 0.0), std::sqrt(0.5), 1e-15);
        shouldEqualTolerance(clebschGordan(0.5, -0.5, 0.5, 0.5, 1.0, 0.0), std::sqrt(0.5), 1e-15);
        shouldEqualTolerance(clebschGordan(0.5, 0.5, 0.5, -0.5, 0.0, 0.0), std::sqrt(0.5), 1e-15);
        shouldEqualTolerance(clebschGordan(0.5, -0.5, 0.5, 0.5, 0.0, 0.0), -std::sqrt(0.5), 1e-15);

        shouldEqualTolerance(clebschGordan(2.0, 2.0, 0.5, 0.5, 2.5, 2.5), std::sqrt(1.0), 1e-15);
        shouldEqualTolerance(clebschGordan(2.0, 2.0, 0.5, -0.5, 2.5, 1.5), std::sqrt(0.2), 1e-15);
        shouldEqualTolerance(clebschGordan(2.0, 1.0, 0.5, 0.5, 2.5, 1.5), std::sqrt(0.8), 1e-15);
        shouldEqualTolerance(clebschGordan(2.0, 2.0, 0.5, -0.5, 1.5, 1.5), std::sqrt(0.8), 1e-15);
        shouldEqualTolerance(clebschGordan(2.0, 1.0, 0.5, 0.5, 1.5, 1.5), -std::sqrt(0.2), 1e-15);
    }
};

struct ClebschGordanTestSuite
: public vigra::test_suite
{
    ClebschGordanTestSuite()
    : vigra::test_suite("ClebschGordanTest")
    {
        add( testCase(&ClebschGordanTest::test));
    }
};

int main(int argc, char ** argv)
{
    ClebschGordanTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}
