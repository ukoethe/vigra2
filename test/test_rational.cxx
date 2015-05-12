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
#include <vigra2/rational.hxx>

using namespace vigra;

struct RationalTest
{
    typedef Rational<int> R;

    void testGcdLcm()
    {
        shouldEqual(gcd(24, 18), 6);
        shouldEqual(lcm(6, 4), 12);
        shouldEqual(gcd(18, 24), 6);
        shouldEqual(lcm(4, 6), 12);
    }

    void testOStreamShifting()
    {
        std::ostringstream out;
        out << R(1,2);
        out << "Testing.." << R(42,23) << 3.141592653589793238 << std::endl;
    }

    void testOperators()
    {
        shouldEqual(R(3,4), R(3,4));
        shouldEqual(-R(3,4), R(-3,4));

        shouldEqual(R(3,4) + R(12,6), R(11,4));
        shouldEqual(R(3,4) - R(12,6), R(-5,4));
        shouldEqual(R(3,4) * R(12,6), R(3,2));
        shouldEqual(R(3,4) / R(12,6), R(3,8));

        shouldEqual(abs(R(-3,4)), R(3,4));
        shouldEqual(norm(R(-3,4)), R(3,4));
        shouldEqual(squaredNorm(R(-3,4)), R(9,16));

        should(R(3,4) == R(9,12));
        should(R(3,4) != R(12,6));
        should(R(3,4) < R(12,6));
        should(R(19,4) > R(12,6));
        should(R(3,4) <= R(12,6));
        should(R(19,4) >= R(12,6));

        shouldEqual(R(3,4) + 2, R(11,4));
        shouldEqual(R(3,4) - 2, R(-5,4));
        shouldEqual(R(3,4) * 2, R(3,2));
        shouldEqual(R(3,4) / 2, R(3,8));
        should(!(R(3,4) == 2));
        should(R(3,4) != 2);
        should(R(3,4) < 2);
        should(R(19,4) > 2);
        should(R(3,4) <= 2);
        should(R(19,4) >= 2);

        shouldEqual(2 + R(3,4), R(11,4));
        shouldEqual(2 - R(3,4), R(5,4));
        shouldEqual(2 * R(3,4), R(3,2));
        shouldEqual(2 / R(3,4), R(8, 3));
        should(!(2 == R(3,4)));
        should(2 != R(3,4));
        should(2 > R(3,4));
        should(2 < R(19,4));
        should(2 >= R(3,4));
        should(2 <= R(19,4));
    }

    void testConversion()
    {
        shouldEqual(rational_cast<R>(R(3,2)), R(3,2));
        shouldEqual(rational_cast<int>(R(3,2)), 1);
        shouldEqual(rational_cast<double>(R(3,2)), 1.5);
        shouldEqual(rational_cast<double>(1.5), 1.5);

        shouldEqual(R(Rational<short>((short)-2, (short)-4)), R(1,2));

        shouldEqual(R(3.5, 1e-4), R(7,2));
        shouldEqual(R(-3.5, 1e-4), R(-7,2));
        shouldEqual(R(0.123, 1e-4), R(123,1000));
        shouldEqual(R(-0.123, 1e-4), R(-123,1000));
        shouldEqual(R(0.123456, 1e-4), R(1235,10000));
        shouldEqual(R(0.123432, 1e-4), R(1234,10000));
        shouldEqual(R(-0.123456, 1e-4), R(-1235,10000));
        shouldEqual(R(-0.123432, 1e-4), R(-1234,10000));
    }

    void testFunctions()
    {
        shouldEqual(pow(R(1,2),2), R(1,4));
        shouldEqual(pow(R(2),-2), R(1,4));
        shouldEqual(pow(R(-1,2),2), R(1,4));
        shouldEqual(pow(R(-2),-2), R(1,4));
        shouldEqual(pow(R(-1,2),3), R(-1,8));
        shouldEqual(pow(R(-2),-3), R(-1,8));
        shouldEqual(pow(R(3),0), R(1));
        shouldEqual(pow(R(0),3), R(0));
        shouldEqual(pow(R(0),0), R(1));
        should(pow(R(0),-3).is_pinf());

        should(pow(R(1,0, false), 1).is_pinf());
        should(pow(R(-1,0, false), 1).is_ninf());
        shouldEqual(pow(R(1,0, false), -1), R(0));
        shouldEqual(pow(R(-1,0, false), -1), R(0));
        try { pow(R(1,0, false), 0); failTest("No exception thrown"); } catch(bad_rational &) {}
        try { pow(R(-1,0, false), 0); failTest("No exception thrown"); } catch(bad_rational &) {}

        shouldEqual(floor(R(2)), R(2));
        shouldEqual(floor(R(3,2)), R(1));
        shouldEqual(floor(R(1,2)), R(0));
        shouldEqual(floor(R(-1,2)), R(-1));
        shouldEqual(floor(R(1,-2)), R(-1));
        shouldEqual(floor(R(-3,2)), R(-2));
        shouldEqual(floor(R(-2)), R(-2));
        shouldEqual(floor(R(1,0,false)), R(1,0,false));
        shouldEqual(floor(R(-1,0,false)), R(-1,0,false));

        shouldEqual(ceil(R(2)), R(2));
        shouldEqual(ceil(R(3,2)), R(2));
        shouldEqual(ceil(R(1,2)), R(1));
        shouldEqual(ceil(R(-1,2)), R(0));
        shouldEqual(ceil(R(1,-2)), R(0));
        shouldEqual(ceil(R(-3,2)), R(-1));
        shouldEqual(ceil(R(-2)), R(-2));
        shouldEqual(ceil(R(1,0,false)), R(1,0,false));
        shouldEqual(ceil(R(-1,0,false)), R(-1,0,false));
    }

    void testInf()
    {
        R inf(2,0);
        R ninf(-2,0);

        should(inf.is_inf());
        should(inf.is_pinf());
        should(!inf.is_ninf());
        should(ninf.is_inf());
        should(ninf.is_ninf());
        should(!ninf.is_pinf());
        shouldEqual(inf.numerator(), 1);
        shouldEqual(ninf.numerator(), -1);

        should((inf + R(1)).is_pinf());
        should((inf + R(0)).is_pinf());
        should((inf + R(-1)).is_pinf());
        should((ninf + R(1)).is_ninf());
        should((ninf + R(0)).is_ninf());
        should((ninf + R(-1)).is_ninf());
        should((inf + 1).is_pinf());
        should((inf + 0).is_pinf());
        should((inf + (-1)).is_pinf());
        should((ninf + 1).is_ninf());
        should((ninf + 0).is_ninf());
        should((ninf + (-1)).is_ninf());
        should((inf + inf).is_pinf());
        should((ninf + ninf).is_ninf());
        shouldEqual((inf + R(3)).numerator(), 1);
        shouldEqual((ninf + R(3)).numerator(), -1);

        should((inf - R(1)).is_pinf());
        should((inf - R(0)).is_pinf());
        should((inf - R(-1)).is_pinf());
        should((ninf - R(1)).is_ninf());
        should((ninf - R(0)).is_ninf());
        should((ninf - R(-1)).is_ninf());
        should((inf - 1).is_pinf());
        should((inf - 0).is_pinf());
        should((inf - (-1)).is_pinf());
        should((ninf - 1).is_ninf());
        should((ninf - 0).is_ninf());
        should((ninf - (-1)).is_ninf());
        should((inf - ninf).is_pinf());
        should((ninf - inf).is_ninf());
        shouldEqual((inf - R(3)).numerator(), 1);
        shouldEqual((ninf - R(3)).numerator(), -1);

        should((inf * R(1)).is_pinf());
        should((inf * R(-1)).is_ninf());
        should((ninf * R(1)).is_ninf());
        should((ninf * R(-1)).is_pinf());
        should((inf * 1).is_pinf());
        should((inf * (-1)).is_ninf());
        should((ninf * 1).is_ninf());
        should((ninf * (-1)).is_pinf());
        should((inf * inf).is_pinf());
        should((inf * ninf).is_ninf());
        should((ninf * inf).is_ninf());
        should((ninf * ninf).is_pinf());
        shouldEqual((inf * R(3)).numerator(), 1);
        shouldEqual((ninf * R(3)).numerator(), -1);
        shouldEqual((inf * R(-3)).numerator(), -1);
        shouldEqual((ninf * R(-3)).numerator(), 1);

        should((inf / R(1)).is_pinf());
        should((inf / R(0)).is_pinf());
        should((inf / R(-1)).is_ninf());
        should((ninf / R(1)).is_ninf());
        should((ninf / R(0)).is_ninf());
        should((ninf / R(-1)).is_pinf());
        shouldEqual(R(1) / inf, R(0));
        shouldEqual(R(-1) / inf, R(0));
        shouldEqual(R(1) / ninf, R(0));
        shouldEqual(R(-1) / ninf, R(0));
        should((inf / 1).is_pinf());
        should((inf / 0).is_pinf());
        should((inf / (-1)).is_ninf());
        should((ninf / 1).is_ninf());
        should((ninf / 0).is_ninf());
        should((ninf / (-1)).is_pinf());

        shouldEqual(2 / inf, R(0));
        shouldEqual((-2) / inf, R(0));
        shouldEqual(2 / ninf, R(0));
        shouldEqual((-2) / ninf, R(0));
        shouldEqual((2 / inf).denominator(), 1);
        shouldEqual(((-2) / inf).denominator(), 1);
        shouldEqual((2 / ninf).denominator(), 1);
        shouldEqual(((-2) / ninf).denominator(), 1);

        shouldEqual((inf / R(3)).numerator(), 1);
        shouldEqual((ninf / R(3)).numerator(), -1);
        shouldEqual((inf / R(-3)).numerator(), -1);
        shouldEqual((ninf / R(-3)).numerator(), 1);

        should(inf == inf);
        should(!(inf != inf));
        should(!(inf < inf));
        should(inf <= inf);
        should(!(inf > inf));
        should(inf >= inf);
        should(ninf == ninf);
        should(!(ninf != ninf));
        should(!(ninf < ninf));
        should(ninf <= ninf);
        should(!(ninf > ninf));
        should(ninf >= ninf);
        should(inf != ninf);
        should(ninf != inf);
        should(inf > ninf);
        should(ninf < inf);
        should(!(inf < ninf));
        should(!(ninf > inf));

        should(inf != 0);
        should(ninf != 0);
        should(inf > 0);
        should(inf >= 0);
        should(ninf < 0);
        should(ninf <= 0);
        should(!(0 < ninf));
        should(!(0 > inf));
        should(!(0 <= ninf));
        should(!(0 >= inf));

        should(inf != R(1));
        should(ninf != R(1));
        should(inf > R(1));
        should(inf >= R(1));
        should(ninf < R(1));
        should(ninf <= R(1));
        should(!(R(1) < ninf));
        should(!(R(1) > inf));
        should(!(R(1) <= ninf));
        should(!(R(1) >= inf));

        try { inf + ninf; failTest("No exception thrown"); } catch(bad_rational &) {}
        try { ninf + inf; failTest("No exception thrown"); } catch(bad_rational &) {}
        try { inf - inf; failTest("No exception thrown"); } catch(bad_rational &) {}
        try { ninf - ninf; failTest("No exception thrown"); } catch(bad_rational &) {}
        try { inf * R(0); failTest("No exception thrown"); } catch(bad_rational &) {}
        try { ninf * R(0); failTest("No exception thrown"); } catch(bad_rational &) {}
        try { R(0) * inf; failTest("No exception thrown"); } catch(bad_rational &) {}
        try { R(0) * ninf; failTest("No exception thrown"); } catch(bad_rational &) {}
        try { inf * 0; failTest("No exception thrown"); } catch(bad_rational &) {}
        try { ninf * 0; failTest("No exception thrown"); } catch(bad_rational &) {}
        try { 0 * inf; failTest("No exception thrown"); } catch(bad_rational &) {}
        try { 0 * ninf; failTest("No exception thrown"); } catch(bad_rational &) {}
        try { inf / inf; failTest("No exception thrown"); } catch(bad_rational &) {}
        try { inf / ninf; failTest("No exception thrown"); } catch(bad_rational &) {}
        try { ninf / inf; failTest("No exception thrown"); } catch(bad_rational &) {}
        try { R(0) / R(0); failTest("No exception thrown"); } catch(bad_rational &) {}
        try { R(0) / 0; failTest("No exception thrown"); } catch(bad_rational &) {}
        try { 0 / R(0); failTest("No exception thrown"); } catch(bad_rational &) {}
    }
};

struct RationalTestSuite
: public vigra::test_suite
{
    RationalTestSuite()
    : vigra::test_suite("RationalTest")
    {
        add( testCase(&RationalTest::testGcdLcm));
        add( testCase(&RationalTest::testOStreamShifting));
        add( testCase(&RationalTest::testOperators));
        add( testCase(&RationalTest::testConversion));
        add( testCase(&RationalTest::testFunctions));
        add( testCase(&RationalTest::testInf));
    }
};

int main(int argc, char ** argv)
{
    RationalTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}
