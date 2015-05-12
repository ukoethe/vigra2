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
#include <vigra2/tiny_linalg.hxx>

using namespace vigra;

struct TinyLinalgTest
{
    TinyLinalgTest()
    {
    }
    
    void testOperations()
    {
        TinyArray<int, 2, 3> a = { 0, 1, 2, 3, 4, 5 };
        TinyArray<int, 3> b {4, 3, 2};
        
        // matrix-matrix
        {
            TinyArray<int, 2, 2> dotref = { 5, 14, 14, 50 };
            shouldEqual(dotref, dot(a, transpose(a)));
        }
        
        // vector-matrix
        {
            TinyArray<int, 2> dotref { 7, 34 };
            shouldEqual(dotref, dot(a, b));
            shouldEqual(dotref, dot(b, transpose(a)));
        }
        
        // vector - symmetric matrix
        {
            TinySymmetricView<int, 3> s(a.data());
            TinyArray<int, 3> dotref { 7, 21, 30 };
            shouldEqual(dotref, dot(b, s));
            shouldEqual(dotref, dot(s, b));
            shouldEqual(dotref, dot(b, transpose(s)));
        }
        
        // matrix - symmetric matrix
        {
            TinySymmetricView<int, 3> s(a.data());
            TinyArray<int, 2, 3> dotref { 5, 11, 14, 14, 35, 47 };
            shouldEqual(dotref, dot(a, s));
            shouldEqual(transpose(dotref), dot(s, transpose(a)));
        }
        
        // symmetric matrix - symmetric matrix
        {
            TinySymmetricView<int, 3> s(a.data());
            TinyArray<int, 6> dotref { 5, 11, 14, 26, 34, 45 };
            shouldEqual(dotref, dot(s, s));
        }
     }
    
    void testEigenvalues()
    {
        double data2[] = { 2.0, 1.0, 2.0 };
        TinySymmetricView<double, 2> s2 = data2;
        
        TinyArray<double, 2> ev2{3.0, 1.0};
        should(max(abs(ev2 - symmetricEigenvalues(s2))) < 1e-15);
        
        double data3[] = { 2.0, 0.0, 1.0, 2.0, 0.0, 2.0 };
        TinySymmetricView<double, 3> s3 = data3;
        
        TinyArray<double, 3> ev3{3.0, 2.0, 1.0};
        should(max(abs(ev3 - symmetricEigenvalues(s3))) < 1e-15);
    }
    
    void testException()
    {
        try
        {
            failTest("no exception thrown");
        }
        catch(std::runtime_error & e)
        {
            std::string expected("expected message");
            std::string message(e.what());
            should(0 == expected.compare(message.substr(0,expected.size())));
        }
    }
};

struct TinyLinalgTestSuite
: public vigra::test_suite
{
    TinyLinalgTestSuite()
    : vigra::test_suite("TinyLinalgTest")
    {
        add( testCase(&TinyLinalgTest::testOperations));
        add( testCase(&TinyLinalgTest::testEigenvalues));
    }
};

int main(int argc, char ** argv)
{
    TinyLinalgTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}
