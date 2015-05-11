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
#include <vigra2/splines.hxx>

using namespace vigra;

template <int ORDER>
struct SplineTest
{
    typedef BSpline<ORDER, double> BS;
    typedef BSplineBase<ORDER, double> BSB;
    BS spline;
    BSB splineBase;

    void testValues()
    {
        double r = spline.radius();
        shouldEqual(r, splineBase.radius());

        for(int d = 0; d <= ORDER+1; ++d)
        {
            for(double x = -r-0.5; x <= r+0.5; x += 0.5)
                shouldEqualTolerance(spline(x, d), splineBase(x, d), 1e-15);
        }
    }

    // void testFixedPointValues()
    // {
        // double r = spline.radius();
        // shouldEqual(r, splineBase.radius());

        // for(double x = -r-0.5; x <= r+0.5; x += 0.5)
        // {
            // FixedPoint<11,20> fpx20(x);
            // FixedPoint<11,15> fpx15(x);
            // should(abs(fixed_point_cast<double>(spline(fpx20)) - spline(x)) < 1e-6);
            // should(abs(fixed_point_cast<double>(spline(fpx15)) - spline(x)) < 4e-5);

        // }
    // }

    void testPrefilterCoefficients()
    {
        int n = ORDER / 2;
        std::vector<double> const & ps = spline.prefilterCoefficients();
        std::vector<double> const & psb = splineBase.prefilterCoefficients();

        if(n == 0)
        {
            shouldEqual(ps.size(), 0u);
            shouldEqual(psb.size(), 0u);
        }
        else
        {
            std::vector<double> & psb1 =
                const_cast<std::vector<double> &>(psb);
            std::sort(psb1.begin(), psb1.end());

            for(int i = 0; i < n; ++i)
                shouldEqualTolerance(ps[i], psb[i], 1e-14);
        }
    }

    void testWeightMatrix()
    {
        int n = ORDER + 1;
        typename BS::WeightMatrix const & ws = BS::weights();
        typename BSB::WeightMatrix const & wsb = BSB::weights();

        for(int d = 0; d < n; ++d)
            for(int i = 0; i < n; ++i)
                shouldEqualTolerance(ws[d][i], wsb[d][i], 1e-14);
    }
};

struct SplineTestSuite
: public vigra::test_suite
{
    SplineTestSuite()
    : vigra::test_suite("SplineTest")
    {
        add( testCase(&SplineTest<0>::testValues));
        add( testCase(&SplineTest<1>::testValues));
        add( testCase(&SplineTest<2>::testValues));
        add( testCase(&SplineTest<3>::testValues));
        add( testCase(&SplineTest<4>::testValues));
        add( testCase(&SplineTest<5>::testValues));
        // add( testCase(&SplineTest<0>::testFixedPointValues));
        // add( testCase(&SplineTest<1>::testFixedPointValues));
        // add( testCase(&SplineTest<2>::testFixedPointValues));
        // add( testCase(&SplineTest<3>::testFixedPointValues));
        // add( testCase(&SplineTest<0>::testPrefilterCoefficients));
        add( testCase(&SplineTest<1>::testPrefilterCoefficients));
        add( testCase(&SplineTest<2>::testPrefilterCoefficients));
        add( testCase(&SplineTest<3>::testPrefilterCoefficients));
        add( testCase(&SplineTest<4>::testPrefilterCoefficients));
        add( testCase(&SplineTest<5>::testPrefilterCoefficients));
        add( testCase(&SplineTest<0>::testWeightMatrix));
        add( testCase(&SplineTest<1>::testWeightMatrix));
        add( testCase(&SplineTest<2>::testWeightMatrix));
        add( testCase(&SplineTest<3>::testWeightMatrix));
        add( testCase(&SplineTest<4>::testWeightMatrix));
        add( testCase(&SplineTest<5>::testWeightMatrix));
    }
};

int main(int argc, char ** argv)
{
    SplineTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}
