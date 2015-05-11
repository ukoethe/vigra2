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
#include <vigra2/autodiff.hxx>
#include <vigra2/splines.hxx>

using namespace vigra;

struct AutodiffTest
{
    typedef autodiff::DualVector<double, 1> N1;
    typedef autodiff::DualVector<double, 2> N2;

    void testOStreamShifting()
    {
        std::ostringstream out;
        out << N1(1.0,2.0);
        out << "Testing.." << N1(42.0,23.0) << 3.141592653589793238 << std::endl;
    }

    void testOperators()
    {
        should(N2(3.0,4.0,5.0).value() == 3.0);
        should(N2(3.0,4.0,5.0).gradient() == N2::Gradient(4.0,5.0));
        should(N2(3.0,4.0,5.0) == N2(3.0, N2::Gradient(4.0,5.0)));

        should(N1(3.0,4.0) == N1(3.0,4.0));
        should(!(N1(3.0,4.0) == N1(2.0,4.0)));
        should(!(N1(3.0,4.0) == N1(3.0,2.0)));
        should(!(N1(3.0,4.0) != N1(3.0,4.0)));
        should(N1(3.0,4.0) != N1(2.0,4.0));
        should(N1(3.0,4.0) != N1(3.0,2.0));
        should(closeAtTolerance(N1(3.0,4.0), N1(3.0,4.0)));
        should(!closeAtTolerance(N1(3.0,4.0), N1(2.0,4.0)));
        should(!closeAtTolerance(N1(3.0,4.0), N1(3.0,2.0)));
        should(closeAtTolerance(N1(3.0,4.0), N1(2.0,4.0), 1.0));
        should(closeAtTolerance(N1(3.0,4.0), N1(3.0,2.0), 1.0));

        shouldEqual(N1(3.0,4.0), N1(3.0,4.0));
        shouldEqual(-N1(3.0,4.0), N1(-3.0,-4.0));

        TinyArray<N2, 2> v = autodiff::dualMatrix(TinyArray<double, 2>(2.0, 3.0));
        shouldEqual(v[0], N2(2.0, 1.0, 0.0));
        shouldEqual(v[1], N2(3.0, 0.0, 1.0));

        shouldEqual(N2(5.0,1.0,0.0) + N2(2.0,0.0,1.0), N2(7.0,1.0,1.0));
        shouldEqual(N2(5.0,1.0,0.0) - N2(2.0,0.0,1.0), N2(3.0,1.0,-1.0));
        shouldEqual(N2(5.0,1.0,0.0) * N2(2.0,0.0,1.0), N2(10.0,2.0,5.0));
        shouldEqual(N2(5.0,1.0,0.0) / N2(2.0,0.0,1.0), N2(2.5,0.5,-1.25));

        shouldEqual(5.0 + N2(2.0,0.0,1.0), N2(7.0,0.0,1.0));
        shouldEqual(5.0 - N2(2.0,0.0,1.0), N2(3.0,0.0,-1.0));
        shouldEqual(5.0 * N2(2.0,0.0,1.0), N2(10.0,0.0,5.0));
        shouldEqual(5.0 / N2(2.0,0.0,1.0), N2(2.5,0.0,-1.25));

        shouldEqual(N2(5.0,1.0,0.0) + 2.0, N2(7.0,1.0,0.0));
        shouldEqual(N2(5.0,1.0,0.0) - 2.0, N2(3.0,1.0,0.0));
        shouldEqual(N2(5.0,1.0,0.0) * 2.0, N2(10.0,2.0,0.0));
        shouldEqual(N2(5.0,1.0,0.0) / 2.0, N2(2.5,0.5,0.0));

        shouldEqual(abs(N1(3.0,4.0)), N1(3.0,4.0));
        shouldEqual(abs(N1(-3.0,4.0)), N1(3.0,-4.0));
        shouldEqual(abs(N1(-3.0,-4.0)), N1(3.0,4.0));
        shouldEqual(abs(N1(3.0,-4.0)), N1(3.0,-4.0));

        shouldEqual(max(N2(5.0,1.0,0.0), N2(2.0,0.0,1.0)), N2(5.0,1.0,0.0));
        shouldEqual(min(N2(5.0,1.0,0.0), N2(2.0,0.0,1.0)), N2(2.0,0.0,1.0));
        shouldEqual(max(5.0, N2(2.0,0.0,1.0)), N2(5.0,0.0,0.0));
        shouldEqual(min(5.0, N2(2.0,0.0,1.0)), N2(2.0,0.0,1.0));
        shouldEqual(max(N2(5.0,1.0,0.0), 2.0), N2(5.0,1.0,0.0));
        shouldEqual(min(N2(5.0,1.0,0.0), 2.0), N2(2.0,0.0,0.0));

        should(N2(1.0, 1.0, 0.0) < N2(2.0, 0.0, 1.0));
        should(!(N2(1.0, 1.0, 0.0) < N2(1.0, 0.0, 1.0)));
        should(N2(1.0, 1.0, 0.0) < 1.1);
        should(!(N2(1.0, 1.0, 0.0) < 1.0));
        should(0.8 < N2(2.0, 0.0, 1.0));
        should(!(2.0 < N2(2.0, 0.0, 1.0)));

        should(N2(1.0, 1.0, 0.0) <= N2(2.0, 0.0, 1.0));
        should(N2(1.0, 1.0, 0.0) <= N2(1.0, 0.0, 1.0));
        should(N2(1.0, 1.0, 0.0) <= 1.1);
        should(N2(1.0, 1.0, 0.0) <= 1.0);
        should(0.8 <= N2(2.0, 0.0, 1.0));
        should(2.0 <= N2(2.0, 0.0, 1.0));

        should(N2(2.0, 1.0, 0.0) > N2(1.0, 0.0, 1.0));
        should(!(N2(1.0, 1.0, 0.0) > N2(1.0, 0.0, 1.0)));
        should(N2(2.0, 1.0, 0.0) > 1.1);
        should(!(N2(1.0, 1.0, 0.0) > 1.0));
        should(2.8 > N2(2.0, 0.0, 1.0));
        should(!(2.0 > N2(2.0, 0.0, 1.0)));

        should(N2(2.0, 1.0, 0.0) >= N2(1.0, 0.0, 1.0));
        should(N2(1.0, 1.0, 0.0) >= N2(1.0, 0.0, 1.0));
        should(N2(2.0, 1.0, 0.0) >= 1.1);
        should(N2(1.0, 1.0, 0.0) >= 1.0);
        should(2.8 >= N2(2.0, 0.0, 1.0));
        should(2.0 >= N2(2.0, 0.0, 1.0));
    }

    void testFunctions()
    {
        // check numbers
        should(closeAtTolerance(log(N1(M_E, 1.0)), N1(1.0, 1.0 / M_E)));
        should(closeAtTolerance(exp(N1(0.0, 1.0)), N1(1.0, 1.0)));
        should(closeAtTolerance(sqrt(N1(4.0, 1.0)), N1(2.0, 0.25)));
        should(closeAtTolerance(sq(N1(3.0, 1.0)), N1(9.0, 6.0)));
        should(closeAtTolerance(sin(N1(M_PI_2, 1.0)), N1(1.0, 0.0)));
        should(closeAtTolerance(sin(N1(M_PI, 1.0)), N1(0.0, -1.0)));
        should(closeAtTolerance(sin_pi(N1(0.5, 1.0)), N1(1.0, 0.0)));
        should(closeAtTolerance(sin_pi(N1(1.0, 1.0)), N1(0.0, -M_PI)));
        should(closeAtTolerance(cos(N1(0.0, 1.0)), N1(1.0, 0.0)));
        should(closeAtTolerance(cos(N1(M_PI_2, 1.0)), N1(0.0, -1.0)));
        should(closeAtTolerance(cos_pi(N1(1.0, 1.0)), N1(-1.0, 0.0)));
        should(closeAtTolerance(cos_pi(N1(0.5, 1.0)), N1(0.0, -M_PI)));
        should(closeAtTolerance(asin(N1(0.5, 1.0)), N1(M_PI/6.0, 2.0/sqrt(3.0)), 1e-15));
        should(closeAtTolerance(acos(N1(0.5, 1.0)), N1(M_PI/3.0, -2.0/sqrt(3.0)), 1e-15));
        should(closeAtTolerance(tan(N1(M_PI/4.0, 1.0)), N1(1.0, 2.0)));
        should(closeAtTolerance(atan(N1(1.0, 1.0)), N1(0.25*M_PI, 0.5)));
        should(closeAtTolerance(sinh(N1(1.0, 1.0)), N1(sinh(1.0), cosh(1.0))));
        should(closeAtTolerance(cosh(N1(1.0, 1.0)), N1(cosh(1.0), sinh(1.0))));
        should(closeAtTolerance(tanh(N1(0.0, 1.0)), N1(0.0, 1.0)));
        should(closeAtTolerance(atan2(N2(-1.0, 1.0, 0.0), N2(-1.0, 0.0, 1.0)), N2(-0.75*M_PI, -0.5, 0.5)));
        should(closeAtTolerance(pow(N1(3.0, 1.0), 2.0), N1(9.0, 6.0)));
        should(closeAtTolerance(pow(3.0, N1(2.0, 1.0)), N1(9.0, 9.0*log(3.0))));
        should(closeAtTolerance(pow(N2(3.0, 1.0, 0.0), N2(2.0, 0.0, 1.0)), N2(9.0, 6.0, 9.0*log(3.0))));
        
        // check constraints
        N1 x(2.3, 1.0);
        N2 a(1.2,2.3,3.4), b(4.5,5.6,6.7);
        should(closeAtTolerance(sq(sqrt(a)), a));
        should(closeAtTolerance(exp(log(a)), a));
        should(closeAtTolerance(sq(sin(x)) + sq(cos(x)), N1(1.0, 0.0), 1e-13));
        should(closeAtTolerance(sin(2.0*a), 2.0*cos(a)*sin(a), 1e-13));
        should(closeAtTolerance(cos(2.0*a), sq(cos(a)) - sq(sin(a)), 1e-13));
        should(closeAtTolerance(sin(a) / cos(a), tan(a), 1e-13));
        should(closeAtTolerance(tan(atan(a)), a, 1e-13));
        should(closeAtTolerance(sq(cosh(x)) - sq(sinh(x)), N1(1.0, 0.0), 1e-13));
        should(closeAtTolerance(tanh(a+b), (tanh(a) + tanh(b)) / (1.0 + tanh(a) * tanh(b)), 1e-12));
        should(closeAtTolerance(atan2(b*sin(a), b*cos(a)), a, 1e-13));
        should(closeAtTolerance(pow(a, 1.0), a));
        should(closeAtTolerance(pow(pow(a, b), 1.0 / b), a, 1e-13));

        BSpline<0, double> s0;
        BSpline<1, double> s1;
        BSpline<2, double> s2;
        BSpline<3, double> s3;
        BSpline<4, double> s4;
        BSpline<5, double> s5;

        for(double x=-3.3; x < 3.5; x += 0.5)
        {
            N1 r = s0(N1(x, 0));
            should(closeAtTolerance(r.value(), s0(x), 1e-15));
            should(closeAtTolerance(r.gradient()[0], s0(x,1), 1e-13));

            r = s1(N1(x, 0));
            should(closeAtTolerance(r.value(), s1(x), 1e-15));
            should(closeAtTolerance(r.gradient()[0], s1(x,1), 1e-13));

            r = s2(N1(x, 0));
            should(closeAtTolerance(r.value(), s2(x), 1e-15));
            should(closeAtTolerance(r.gradient()[0], s2(x,1), 1e-13));

            r = s3(N1(x, 0));
            should(closeAtTolerance(r.value(), s3(x), 1e-15));
            should(closeAtTolerance(r.gradient()[0], s3(x,1), 1e-13));

            r = s4(N1(x, 0));
            should(closeAtTolerance(r.value(), s4(x), 1e-15));
            should(closeAtTolerance(r.gradient()[0], s4(x,1), 1e-13));

            r = s5(N1(x, 0));
            should(closeAtTolerance(r.value(), s5(x), 1e-15));
            should(closeAtTolerance(r.gradient()[0], s5(x,1), 1e-13));
        }
    }
};

struct AutodiffTestSuite
: public vigra::test_suite
{
    AutodiffTestSuite()
    : vigra::test_suite("AutodiffTest")
    {
        add( testCase(&AutodiffTest::testOStreamShifting));
        add( testCase(&AutodiffTest::testOperators));
        add( testCase(&AutodiffTest::testFunctions));
    }
};

int main(int argc, char ** argv)
{
    AutodiffTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}
