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
#include <vigra2/bessel.hxx>

using namespace vigra;

struct BesselTest
{
    BesselTest()
    {
    }
    
    void test()
    {
        // Reference values computed to 16 digits with Python.mpmath. 
        // Casual comparison showed no difference to Mathematica's results.
        double x[] = { 1.0, 4.0, 6.0 };
        double besseljnref[] = {
            7.6519768655796649e-01, -3.9714980986384740e-01, 1.5064525725099692e-01, 
            5.7672480775687329e-01, -3.2757913759146523e-01, -4.6828234823458334e-03, 
            4.8609126058589103e-01, -2.4287320996018547e-01, -1.1299172042407525e-01, 
            4.3017147387562193e-01, -1.6755558799533424e-01, -1.8093519033665686e-01, 
            3.9123236045864818e-01, -1.0535743487538894e-01, -2.1960268610200856e-01, 
            3.6208707488717234e-01, -5.5038855669513713e-02, -2.3828585178317879e-01, 
            3.3919660498317961e-01, -1.4458842084785106e-02, -2.4372476722886663e-01, 
            3.2058907797982628e-01, 1.8376032647858614e-02, -2.4057094958616052e-01, 
            3.0506707225300012e-01, 4.5095329080457235e-02, -2.3197310306707983e-01, 
            2.9185568526512001e-01, 6.6976198673670620e-02, -2.2004622511384700e-01, 
            2.8042823052537585e-01, 8.5006705446061023e-02, -2.0620569442259729e-01, 
            2.7041248255096445e-01, 9.9950477050301592e-02, -1.9139539469541733e-01, 
            2.6153687541034509e-01, 1.1240023492610679e-01, -1.7624117645477547e-01, 
            2.5359797330294920e-01, 1.2281915265293869e-01, -1.6115376768165826e-01, 
            2.4643993656993257e-01, 1.3157198580936999e-01, -1.4639794400255970e-01    
        };
        double besselynref[] = {
            8.8256964215676956e-02, -1.6940739325064992e-02, -2.8819468398157916e-01, 
            -1.0703243154093754e-01, 1.4786314339122683e-01, -3.0266723702418485e-01, 
            -1.6040039348492374e-01, 2.2985790254811306e-01, -2.6303660482037811e-01, 
            -1.8202211595348500e-01, 2.6808060304231507e-01, -2.0509487811877961e-01, 
            -1.9214228737369318e-01, 2.8294322431117191e-01, -1.4494951186809379e-01, 
            -1.9706088806443733e-01, 2.8511777841103764e-01, -8.9252841434580163e-02, 
            -1.9930679029227036e-01, 2.8035255955745608e-01, -4.0297251103395833e-02, 
            -2.0006390460040860e-01, 2.7184139484930947e-01, 1.5698795407253514e-03, 
            -1.9994686666043449e-01, 2.6140472921203017e-01, 3.6815736940746704e-02, 
            -1.9929926580524435e-01, 2.5009898312668521e-01, 6.6197858895869655e-02, 
            -1.9832403085028555e-01, 2.3854272714494473e-01, 9.0526604143921052e-02, 
            -1.9714613354518651e-01, 2.2709735924007149e-01, 1.1056356972736049e-01, 
            -1.9584504763522584e-01, 2.1597027298252575e-01, 1.2698414345087472e-01, 
            -1.9447256680104227e-01, 2.0527533641239212e-01, 1.4036965442780550e-01, 
            -1.9306306446008192e-01, 1.9506914688206353e-01, 1.5121244335755843e-01        
        };

        for(int n = 0; n < 15; ++n)
        {
            if(n == 0)
                shouldEqual(besselJ(n, 0.0), 1.0);
            else
                shouldEqual(besselJ(n, 0.0), 0.0);
            should(besselY(n, 0.0) == -std::numeric_limits<double>::infinity());

            for(int k=0; k<3; ++k)
            {
                double f = odd(n) ? -1.0 : 1.0;
                double eps = 1e-14;
                shouldEqualTolerance(besselJ(n, x[k]+n) - besseljnref[k+3*n], 0.0, eps);
                shouldEqualTolerance(besselJ(-n, x[k]+n) - f*besseljnref[k+3*n], 0.0, eps);
                shouldEqualTolerance(besselY(n, x[k]+n) - besselynref[k+3*n], 0.0, eps);
                shouldEqualTolerance(besselY(-n, x[k]+n) - f*besselynref[k+3*n], 0.0, eps);
            }
        }
    }
};

struct BesselTestSuite
: public vigra::test_suite
{
    BesselTestSuite()
    : vigra::test_suite("BesselTest")
    {
        add( testCase(&BesselTest::test));
    }
};

int main(int argc, char ** argv)
{
    BesselTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}
