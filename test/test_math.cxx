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
#include <vigra2/math.hxx>
// #include <vigra2/array_vector.hxx>
// #include <vigra2/fixedpoint.hxx>
// #include <vigra2/linear_algebra.hxx>
// #include <vigra2/singular_value_decomposition.hxx>
// #include <vigra2/regression.hxx>
// #include <vigra2/random.hxx>
// #include <vigra2/quaternion.hxx>
// #include <vigra2/timing.hxx>

using namespace vigra;

struct MathTest
{
    void testSpecialIntegerFunctions()
    {
        for(int32_t i = 0; i < 1024; ++i)
        {
            shouldEqual(sqrti(i), (int32_t)floor(sqrt((double)i)));
        }

        shouldEqual(roundi(0.0), 0);
        shouldEqual(roundi(1.0), 1);
        shouldEqual(roundi(1.1), 1);
        shouldEqual(roundi(1.6), 2);
        shouldEqual(roundi(-1.0), -1);
        shouldEqual(roundi(-1.1), -1);
        shouldEqual(roundi(-1.6), -2);

        uint32_t roundPower2[] = {0, 1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 0xffff, 0x7fffffff, 0x80000000, 0x80000001, 0xffffffff};
        uint32_t floorResult[] = {0, 1, 2, 2, 4, 4, 4, 8, 8, 8, 16, 0x8000, 0x40000000, 0x80000000, 0x80000000, 0x80000000};
        uint32_t ceilResult[] = {0, 1, 2, 4, 4, 8, 8, 8, 16, 16, 16, 0x10000, 0x80000000, 0x80000000, 0, 0};
        for(unsigned int i = 0; i < sizeof(roundPower2) / sizeof(uint32_t); ++i)
        {
            shouldEqual(floorPower2(roundPower2[i]), floorResult[i]);
            shouldEqual(ceilPower2(roundPower2[i]), ceilResult[i]);
        }

        for(int32_t k=0; k<32; ++k)
        {
            shouldEqual(log2i(1 << k), k);
            shouldEqual(log2i((1 << k) + 1), k == 0 ? 1 : k);
            shouldEqual(log2i((1 << k) - 1), k-1);
        }

        should(even(0));
        should(!odd(0));
        should(!even(1));
        should(odd(1));
        should(even(2));
        should(!odd(2));
        should(!even(-1));
        should(odd(-1));
        should(even(-2));
        should(!odd(-2));
    }


    void testSpecialFunctions()
    {
        shouldEqualTolerance(ellipticIntegralE(M_PI / 2.0, 0.0), M_PI / 2.0, 1e-14);
        shouldEqualTolerance(ellipticIntegralF(0.3, 0.3), 0.30039919311549118, 1e-14);
        shouldEqualTolerance(ellipticIntegralE(0.3, 0.3), 0.29960175507025716, 1e-14);

        shouldEqualTolerance(erf(0.3), 0.32862675945912745, 1e-7);

        should(noncentralChi2CDFApprox(200, 0.0, 200.0) > 0.5);
        should(noncentralChi2CDFApprox(200, 0.0, 199.0) < 0.5);
        should(noncentralChi2CDF(200, 0.0, 200.0) > 0.5);
        should(noncentralChi2CDF(200, 0.0, 199.0) < 0.5);

        shouldEqualTolerance(noncentralChi2CDF(2, 2.0, 2.0), 0.34574583872316456, 1e-7);
        shouldEqualTolerance(noncentralChi2(2, 2.0, 2.0), 0.154254161276835, 1e-7);
        shouldEqualTolerance(noncentralChi2CDF(3, 2.0, 2.0), 0.22073308707450343, 1e-7);
        shouldEqualTolerance(noncentralChi2(3, 2.0, 2.0), 0.13846402271767755, 1e-7);
        shouldEqualTolerance(noncentralChi2CDFApprox(2, 2.0, 2.0), 0.34574583872316456, 1e-1);
        shouldEqualTolerance(noncentralChi2CDFApprox(3, 2.0, 2.0), 0.22073308707450343, 1e-1);

        for(double x = -4.0; x <= 4.0; x += 1.0)
        {
            shouldEqual(sin_pi(x), 0.0);
            shouldEqual(cos_pi(x+0.5), 0.0);
        }
        
        for(double x = -4.5; x <= 4.5; x += 2.0)
        {
            shouldEqual(sin_pi(x), -1.0);
            shouldEqual(cos_pi(x+0.5), 1.0);
        }
        
        for(double x = -3.5; x <= 4.5; x += 2.0)
        {
            shouldEqual(sin_pi(x), 1.0);
            shouldEqual(cos_pi(x+0.5), -1.0);
        }
        
        for(double x = -4.0; x <= 4.0; x += 0.0625)
        {
            shouldEqualTolerance(sin_pi(x), std::sin(M_PI*x), 1e-14);
            shouldEqualTolerance(cos_pi(x), std::cos(M_PI*x), 1e-14);
        }

        shouldEqualTolerance(sin_pi(0.25), 0.5*M_SQRT2, 2e-16);
        shouldEqualTolerance(cos_pi(0.25), 0.5*M_SQRT2, 2e-16);

        shouldEqual(gamma(4.0), 6.0);
        shouldEqualTolerance(gamma(0.1), 9.5135076986687306, 1e-15);
        shouldEqualTolerance(gamma(3.2), 2.4239654799353683, 1e-15);
        shouldEqualTolerance(gamma(170.2), 1.1918411166366696e+305, 1e-15);
        shouldEqualTolerance(gamma(-0.1), -10.686287021193193, 1e-14);
        shouldEqualTolerance(gamma(-3.2), 0.689056412005979, 1e-14);
        shouldEqualTolerance(gamma(-170.2), -2.6348340538196879e-307, 1e-14);
        try { gamma(0.0); failTest("No exception thrown"); } catch(ContractViolation &) {}
        try { gamma(-1.0); failTest("No exception thrown"); } catch(ContractViolation &) {}

        shouldEqual(loggamma(1.0), 0.0);
        shouldEqual(loggamma(2.0), 0.0);
        shouldEqualTolerance(loggamma(4.0e-22), 49.2705776847491144296, 1e-15);
        shouldEqualTolerance(loggamma(0.1), 2.2527126517342055401, 1e-15);
        shouldEqualTolerance(loggamma(0.3), 1.0957979948180756047, 1e-15);
        shouldEqualTolerance(loggamma(0.8), 0.15205967839983755563, 1e-15);
        shouldEqualTolerance(loggamma(1.1), -0.049872441259839757344, 1e-15);
        shouldEqualTolerance(loggamma(1.3), -0.10817480950786048655, 1e-15);
        shouldEqualTolerance(loggamma(1.8), -0.071083872914372153717, 1e-15);
        shouldEqualTolerance(loggamma(3.0), 0.69314718055994528623, 1e-15);
        shouldEqualTolerance(loggamma(3.1), 0.78737508327386251938, 1e-15);
        shouldEqualTolerance(loggamma(4.0), 1.79175946922805500081, 1e-15);
        shouldEqualTolerance(loggamma(8.0), 8.5251613610654143002, 1e-15);
        shouldEqualTolerance(loggamma(1000.0), 5905.2204232091812118261, 1e-15);
        shouldEqualTolerance(loggamma(1000.2), 5906.6018942569799037, 1e-15);
        shouldEqualTolerance(loggamma(2.8e+17), 1.096859847946237952e+19, 1e-15);
        shouldEqualTolerance(loggamma(2.9e+17), 1.1370510622188449792e+19, 1e-15);
        shouldEqualTolerance(loggamma(5.7646075230342349e+17), 2.2998295812288974848e+19, 1e-15);
        try { loggamma(0.0); failTest("No exception thrown"); } catch(ContractViolation &) {}
        try { loggamma(-1.0); failTest("No exception thrown"); } catch(ContractViolation &) {}

        double args[5] = {0.0, 1.0, 0.7, -0.7, -1.0};
        for(int i=0; i<5; ++i)
        {
            double x = args[i], x2 = x*x;
            shouldEqualTolerance(legendre(0, x), 1.0, 1e-15);
            shouldEqualTolerance(legendre(1, x), x, 1e-15);
            shouldEqualTolerance(legendre(2, x), 0.5*(3.0*x2-1.0), 1e-15);
            shouldEqualTolerance(legendre(3, x), 0.5*x*(5.0*x2-3.0), 1e-15);

            shouldEqualTolerance(legendre(0, 0, x), 1.0, 1e-15);
            shouldEqualTolerance(legendre(1, 0, x), x, 1e-15);
            shouldEqualTolerance(legendre(1, 1, x), -std::sqrt(1.0-x2), 1e-15);
            shouldEqualTolerance(legendre(2, 0, x), 0.5*(3.0*x2-1.0), 1e-15);
            shouldEqualTolerance(legendre(2, 1, x), -3.0*x*std::sqrt(1.0-x2), 1e-15);
            shouldEqualTolerance(legendre(2, 2, x), 3.0*(1.0-x2), 1e-15);
            shouldEqualTolerance(legendre(4, 2, x), 7.5*(7.0*x2-1.0)*(1.0-x2), 1e-15);
            shouldEqualTolerance(legendre(1, -1, x), -legendre(1, 1, x) / 2.0, 1e-15);
            shouldEqualTolerance(legendre(2, -1, x), -legendre(2, 1, x) / 6.0, 1e-15);
            shouldEqualTolerance(legendre(2, -2, x), legendre(2, 2, x) / 24.0, 1e-15);
        }
    }
};

// struct QuaternionTest
// {
    // typedef Quaternion<double> Q;
    // typedef Q::Vector V;

    // void testContents()
    // {
        // Q q(1.0, 2.0, 3.0, 4.0), q0, q1(-1.0), q2(q), q3(q.w(), q.v());

        // shouldEqual(q.w(), 1.0);
        // shouldEqual(q.v(), V(2.0, 3.0, 4.0));
        // shouldEqual(q0.w(), 0.0);
        // shouldEqual(q0.v(), V(0.0, 0.0, 0.0));
        // shouldEqual(q1.w(), -1.0);
        // shouldEqual(q1.v(), V(0.0, 0.0, 0.0));
        // shouldEqual(q2.w(), 1.0);
        // shouldEqual(q2.v(), V(2.0, 3.0, 4.0));
        // shouldEqual(q3.w(), 1.0);
        // shouldEqual(q3.v(), V(2.0, 3.0, 4.0));

        // shouldEqual(q[0], 1.0);
        // shouldEqual(q[1], 2.0);
        // shouldEqual(q[2], 3.0);
        // shouldEqual(q[3], 4.0);
        // shouldEqual(q.x(), 2.0);
        // shouldEqual(q.y(), 3.0);
        // shouldEqual(q.z(), 4.0);

        // should(q == q2);
        // should(q1 != q2);

        // q2 = q1;
        // shouldEqual(q2.w(), -1.0);
        // shouldEqual(q2.v(), V(0.0, 0.0, 0.0));

        // should(q != q2);
        // should(q1 == q2);

        // q3 = 10.0;
        // shouldEqual(q3.w(), 10.0);
        // shouldEqual(q3.v(), V(0.0, 0.0, 0.0));

        // q2.setW(-2.0);
        // shouldEqual(q2.w(), -2.0);
        // shouldEqual(q2.v(), V(0.0, 0.0, 0.0));

        // q2.setV(V(5.0, 6.0, 7.0));
        // shouldEqual(q2.w(), -2.0);
        // shouldEqual(q2.v(), V(5.0, 6.0, 7.0));

        // q3.setV(5.0, 6.0, 7.0);
        // shouldEqual(q3.w(), 10.0);
        // shouldEqual(q3.v(), V(5.0, 6.0, 7.0));

        // q3.setX(2.0);
        // q3.setY(3.0);
        // q3.setZ(4.0);
        // shouldEqual(q3.w(), 10.0);
        // shouldEqual(q3.v(), V(2.0, 3.0, 4.0));

        // shouldEqual(q.squaredMagnitude(), 30.0);
        // shouldEqual(squaredNorm(q), 30.0);
        // shouldEqualTolerance(q.magnitude(), std::sqrt(30.0), 1e-15);
        // shouldEqualTolerance(norm(q), std::sqrt(30.0), 1e-15);
        // shouldEqual(norm(q), abs(q));
    // }

    // void testStreamIO()
    // {
        // std::ostringstream out;
        // Q q(1.0, 2.0, 3.0, 4.0);

        // out << q;
        // shouldEqual(out.str(), "1 2 3 4");

        // std::istringstream in;
        // in.str("10 11 12 13");
        // in >> q;
        // shouldEqual(q, Q(10.0, 11.0, 12.0, 13.0));
    // }

    // void testOperators()
    // {
        // Q q(1.0, 2.0, 3.0, 4.0);

        // shouldEqual(+q, q);
        // shouldEqual(-q, Q(-1,-2,-3,-4));

        // shouldEqual(q+q, Q(2,4,6,8));
        // shouldEqual(q+2.0, Q(3,2,3,4));
        // shouldEqual(2.0+q, Q(3,2,3,4));

        // shouldEqual(Q(2,4,6,8) - q, q);
        // shouldEqual(q-2.0, Q(-1,2,3,4));
        // shouldEqual(2.0-q, Q(1,-2,-3,-4));

        // shouldEqual(Q(1,0,0,0)*Q(1,0,0,0), Q(1,0,0,0));
        // shouldEqual(Q(0,1,0,0)*Q(0,1,0,0), Q(-1,0,0,0));
        // shouldEqual(Q(0,0,1,0)*Q(0,0,1,0), Q(-1,0,0,0));
        // shouldEqual(Q(0,0,0,1)*Q(0,0,0,1), Q(-1,0,0,0));

        // shouldEqual(Q(0,1,0,0)*Q(0,0,1,0), Q(0,0,0,1));
        // shouldEqual(Q(0,0,1,0)*Q(0,1,0,0), Q(0,0,0,-1));
        // shouldEqual(Q(0,0,1,0)*Q(0,0,0,1), Q(0,1,0,0));
        // shouldEqual(Q(0,0,0,1)*Q(0,0,1,0), Q(0,-1,0,0));
        // shouldEqual(Q(0,0,0,1)*Q(0,1,0,0), Q(0,0,1,0));
        // shouldEqual(Q(0,1,0,0)*Q(0,0,0,1), Q(0,0,-1,0));

        // shouldEqual(q*q, Q(-28,4,6,8));
        // shouldEqual(q*2.0, Q(2,4,6,8));
        // shouldEqual(2.0*q, Q(2,4,6,8));

        // Q q1 = q / q;       
        // shouldEqualTolerance(q1[0], 1.0, 1e-16);
        // shouldEqualTolerance(q1[1], 0.0, 1e-16);
        // shouldEqualTolerance(q1[2], 0.0, 1e-16);
        // shouldEqualTolerance(q1[3], 0.0, 1e-16);
        // shouldEqual(Q(2,4,6,8)/2.0, q);
        // shouldEqual(60.0/q, Q(2,-4,-6,-8));

        // shouldEqualTolerance(norm(q / norm(q)), 1.0, 1e-15);
    // }

    // void testRotation()
    // {
        // Q q(1.0, 2.0, 3.0, 4.0);
        // q /= norm(q);

        // double ref[3][3] = {{-2.0/3.0,  0.4/3.0, 2.2/3.0 }, 
                            // { 2.0/3.0, -1.0/3.0, 2.0/3.0 },
                            // { 1.0/3.0,  2.8/3.0, 0.4/3.0 } };

        // Matrix<double> m(3,3), mref(3,3, (double*)ref);
        // q.fillRotationMatrix(m);
        // shouldEqualSequenceTolerance(m.begin(), m.end(), mref.begin(), 1e-15);

        // double res[3][3];
        // q.fillRotationMatrix(res);
        // shouldEqualSequenceTolerance((double*)res, (double*)res+9, (double*)ref, 1e-15);

        // Q q1 = Q::createRotation(M_PI/2.0, V(1,0,0));
        // Q q2 = Q::createRotation(M_PI/2.0, V(0,1,0));
        // Q q3 = Q::createRotation(M_PI/2.0, V(0,0,1));
        // Q q4 = q3*(-q1)*q2*q1;

        // shouldEqualTolerance(norm(q4), 1.0, 1e-15);
        // shouldEqualTolerance(q4[0], 0.0, 1e-15);
    // }
// };

// struct FixedPointTest
// {
    // void testConstruction()
    // {
        // shouldEqual(fixedPoint(3).value, 3);
        // shouldEqual(fixedPoint(-3).value, -3);
        // shouldEqual(-fixedPoint(3).value, -3);

        // shouldEqual((FixedPoint<3,4>(3).value), 3 << 4);
        // shouldEqual((FixedPoint<3,4>(-3).value), -3 << 4);
        // shouldEqual((-FixedPoint<3,4>(3).value), -3 << 4);

        // shouldEqual((FixedPoint<3,4>(3.5).value), 56);
        // shouldEqual((FixedPoint<3,4>(-3.5).value), -56);
        // shouldEqual((-FixedPoint<3,4>(3.5).value), -56);

        // try { FixedPoint<1, 8>(3.75); failTest("No exception thrown"); } catch(PreconditionViolation &) {}

        // shouldEqual((NumericTraits<FixedPoint<1, 8> >::zero()).value, 0);
        // shouldEqual((NumericTraits<FixedPoint<1, 8> >::one()).value, 1 << 8);
        // shouldEqual((NumericTraits<FixedPoint<1, 8> >::max()).value, (1 << 9) - 1);
        // shouldEqual((NumericTraits<FixedPoint<1, 8> >::min()).value, -((1 << 9) - 1));

        // FixedPoint<2, 8> v(3.75);
        // shouldEqual((FixedPoint<2, 8>(v).value), 15 << 6);
        // shouldEqual((FixedPoint<3, 10>(v).value), 15 << 8);
        // shouldEqual((FixedPoint<2, 2>(v).value), 15);
        // shouldEqual((FixedPoint<2, 0>(v).value), 4);

        // shouldEqual((FixedPoint<2, 8>(-v).value), -15 << 6);
        // shouldEqual((FixedPoint<3, 10>(-v).value), -15 << 8);
        // shouldEqual((FixedPoint<2, 2>(-v).value), -15);
        // shouldEqual((FixedPoint<2, 0>(-v).value), -4);

        // shouldEqual(fixed_point_cast<double>(v), 3.75);
        // should((frac(v) == FixedPoint<0, 8>(0.75)));
        // should((dual_frac(v) == FixedPoint<0, 8>(0.25)));
        // should(floor(v) == 3);
        // should(ceil(v) == 4);
        // should(round(v) == 4);
        // should(abs(v) == v);
        // should((frac(-v) == FixedPoint<0, 8>(0.25)));
        // should((dual_frac(-v) == FixedPoint<0, 8>(0.75)));
        // should(floor(-v) == -4);
        // should(ceil(-v) == -3);
        // should(round(-v) == -4);
        // should(abs(-v) == v);
        // should(norm(-v) == v);
        // should(squaredNorm(-v) == v*v);

        // FixedPoint<3, 10> v1;
        // shouldEqual((v1 = v).value, 15 << 8);
        // shouldEqual((v1 = -v).value, -15 << 8);

        // FixedPoint<2, 0> v2;
        // shouldEqual((v2 = v).value, 4);
        // shouldEqual((v2 = -v).value, -4);
    // }

    // void testComparison()
    // {
        // FixedPoint<3, 8> v1(3.75), v2(4);
        // FixedPoint<2, 2> v3(3.75);
        // should(v1 == v1);
        // should(v1 == v3);
        // should(!(v1 != v1));
        // should(!(v1 != v3));
        // should(v1 <= v1);
        // should(v1 <= v3);
        // should(!(v1 < v1));
        // should(!(v1 < v3));
        // should(v1 >= v1);
        // should(v1 >= v3);
        // should(!(v1 > v1));
        // should(!(v1 > v3));

        // should(v2 != v1);
        // should(v2 != v3);
        // should(!(v2 == v1));
        // should(!(v2 == v3));
        // should(!(v2 <= v1));
        // should(!(v2 <= v3));
        // should(!(v2 < v1));
        // should(!(v2 < v3));
        // should(v2 >= v1);
        // should(v2 >= v3);
        // should(v2 > v1);
        // should(v2 > v3);
    // }

    // void testArithmetic()
    // {
        // FixedPoint<1, 16> t1(0.75), t2(0.25);
        // signed char v1 = 1, v2 = 2, v4 = 4, v8 = 8;

        // should((FixedPoint<1, 16>(t1) += t1) == (FixedPoint<1, 16>(1.5)));
        // should((FixedPoint<1, 16>(t1) -= t1) == (FixedPoint<1, 16>(0.0)));
        // should((FixedPoint<2, 16>(t1) *= t1) == (FixedPoint<1, 16>(9.0 / 16.0)));

        // should(--t1 == (FixedPoint<1, 16>(-0.25)));
        // should(t1 == (FixedPoint<1, 16>(-0.25)));
        // should(++t1 == (FixedPoint<1, 16>(0.75)));
        // should(t1 == (FixedPoint<1, 16>(0.75)));
        // should(t1++ == (FixedPoint<1, 16>(0.75)));
        // should(t1 == (FixedPoint<1, 16>(1.75)));
        // should(t1-- == (FixedPoint<1, 16>(1.75)));
        // should(t1 == (FixedPoint<1, 16>(0.75)));

        // shouldEqual((t1 * fixedPoint(v1)).value, 3 << 14);
        // shouldEqual((t2 * fixedPoint(v1)).value, 1 << 14);
        // shouldEqual((-t1 * fixedPoint(v1)).value, -3 << 14);
        // shouldEqual((-t2 * fixedPoint(v1)).value, -1 << 14);
        // shouldEqual((t1 * -fixedPoint(v1)).value, -3 << 14);
        // shouldEqual((t2 * -fixedPoint(v1)).value, -1 << 14);

        // shouldEqual((FixedPoint<8, 2>(t1 * fixedPoint(v1))).value, 3);
        // shouldEqual((FixedPoint<8, 2>(t2 * fixedPoint(v1))).value, 1);
        // shouldEqual((FixedPoint<8, 2>(-t1 * fixedPoint(v1))).value, -3);
        // shouldEqual((FixedPoint<8, 2>(-t2 * fixedPoint(v1))).value, -1);

        // shouldEqual(floor(t1 * fixedPoint(v1) + t2 * fixedPoint(v2)), 1);
        // shouldEqual(ceil(t1 * fixedPoint(v1) + t2 * fixedPoint(v2)), 2);
        // shouldEqual(round(t1 * fixedPoint(v1) + t2 * fixedPoint(v2)), 1);
        // shouldEqual(floor(t1 * fixedPoint(v4) + t2 * fixedPoint(v8)), 5);
        // shouldEqual(ceil(t1 * fixedPoint(v4) + t2 * fixedPoint(v8)), 5);
        // shouldEqual(round(t1 * fixedPoint(v4) + t2 * fixedPoint(v8)), 5);

        // shouldEqual(floor(t1 * -fixedPoint(v1) - t2 * fixedPoint(v2)), -2);
        // shouldEqual(ceil(t1 * -fixedPoint(v1) - t2 * fixedPoint(v2)), -1);
        // shouldEqual(round(t1 * -fixedPoint(v1) - t2 * fixedPoint(v2)), -1);
        // shouldEqual(floor(t1 * -fixedPoint(v4) - t2 * fixedPoint(v8)), -5);
        // shouldEqual(ceil(t1 * -fixedPoint(v4) - t2 * fixedPoint(v8)), -5);
        // shouldEqual(round(t1 * -fixedPoint(v4) - t2 * fixedPoint(v8)), -5);

        // double d1 = 1.0 / 3.0, d2 = 1.0 / 7.0;
        // FixedPoint<1, 24> r1(d1), r2(d2);
        // FixedPoint<2, 24> r3;
        // add(r1, r2, r3);
        // shouldEqual(r3.value, (FixedPoint<2, 24>(d1 + d2)).value);
        // sub(r1, r2, r3);
        // shouldEqual(r3.value, (FixedPoint<2, 24>(d1 - d2)).value);
        // mul(r1, r2, r3);
        // shouldEqual(r3.value >> 2, (FixedPoint<2, 24>(d1 * d2)).value >> 2);

        // for(int i = 0; i < 1024; ++i)
        // {
            // FixedPoint<4,5> fv1(i, FPNoShift);
            // FixedPoint<5,4> fv2(i, FPNoShift);
            // FixedPoint<5,5> fv3(i, FPNoShift);
            // FixedPoint<6,6> fv4(i, FPNoShift);
            // shouldEqual(fixed_point_cast<double>(sqrt(fv1)), floor(sqrt((double)fv1.value)) / 8.0);
            // shouldEqual(fixed_point_cast<double>(sqrt(fv2)), floor(sqrt((double)fv2.value)) / 4.0);
            // shouldEqual(fixed_point_cast<double>(sqrt(fv3)), floor(sqrt((double)fv3.value)) / 4.0);
            // shouldEqual(fixed_point_cast<double>(sqrt(fv4)), floor(sqrt((double)fv4.value)) / 8.0);
        // }
    // }
// };

// struct FixedPoint16Test
// {
    // void testConstruction()
    // {
        // shouldEqual((FixedPoint16<3>(3).value), 3 << 12);
        // shouldEqual((FixedPoint16<3>(-3).value), -3 << 12);
        // shouldEqual((-FixedPoint16<3>(3).value), -3 << 12);

        // shouldEqual((FixedPoint16<8>(3).value), 3 << 7);

        // shouldEqual((FixedPoint16<3>(3.5).value), 7 << 11);
        // shouldEqual((FixedPoint16<3>(-3.5).value), -(7 << 11));
        // shouldEqual((-FixedPoint16<3>(3.5).value), -(7 << 11));

        // shouldEqual((NumericTraits<FixedPoint16<4> >::zero()).value, 0);
        // shouldEqual((NumericTraits<FixedPoint16<4> >::one()).value, 1 << 11);
        // shouldEqual((NumericTraits<FixedPoint16<4> >::max()).value, (1 << 15) - 1);
        // shouldEqual((NumericTraits<FixedPoint16<4> >::min()).value, -(1 << 15));

        // shouldEqual((FixedPoint16<1, FPOverflowSaturate>(3.75).value), (1 << 15)-1);
        // shouldEqual((FixedPoint16<1, FPOverflowSaturate>(-3.75).value), -(1 << 15));
        // try { FixedPoint16<1, FPOverflowError>(3.75); failTest("No exception thrown"); } 
        // catch(PreconditionViolation &) {}
        // try { FixedPoint16<1, FPOverflowError>(-3.75); failTest("No exception thrown"); } 
        // catch(PreconditionViolation &) {}

        // FixedPoint16<4> v(3.75);
        // shouldEqual((v.value), 15 << 9);
        // shouldEqual(((-v).value), -15 << 9);
        // shouldEqual((FixedPoint16<4>(v).value), 15 << 9);
        // shouldEqual((FixedPoint16<6>(v).value), 15 << 7);
        // shouldEqual((FixedPoint16<13>(v).value), 15);
        // shouldEqual((FixedPoint16<15>(v).value), 4);

        // shouldEqual((FixedPoint16<4>(-v).value), -15 << 9);
        // shouldEqual((FixedPoint16<6>(-v).value), -15 << 7);
        // shouldEqual((FixedPoint16<13>(-v).value), -15);
        // shouldEqual((FixedPoint16<15>(-v).value), -4);

        // shouldEqual(fixed_point_cast<double>(v), 3.75);
        // shouldEqual(fixed_point_cast<double>(-v), -3.75);
        // shouldEqual(frac(v), FixedPoint16<4>(0.75));
        // shouldEqual(dual_frac(v), FixedPoint16<4>(0.25));
        // shouldEqual(frac(-v), FixedPoint16<4>(0.25));
        // shouldEqual(dual_frac(-v), FixedPoint16<4>(0.75));
        // shouldEqual(floor(v), 3);
        // shouldEqual(ceil(v), 4);
        // shouldEqual(floor(-v), -4);
        // shouldEqual(ceil(-v), -3);
        // shouldEqual(round(v), 4);
        // shouldEqual(round(-v), -4);
        // shouldEqual(abs(v), v);
        // shouldEqual(abs(-v), v);
        // shouldEqual(norm(-v), v);
        // shouldEqual(squaredNorm(-v), v*v);

        // FixedPoint16<2> v1;
        // shouldEqual((v1 = v).value, 15 << 11);
        // shouldEqual((v1 = -v).value, -15 << 11);

        // FixedPoint16<15> v2;
        // shouldEqual((v2 = v).value, 4);
        // shouldEqual((v2 = -v).value, -4);
    // }

    // void testComparison()
    // {
        // FixedPoint16<4> v1(3.75), v2(4);
        // FixedPoint16<2> v3(3.75);
        // should(v1 == v1);
        // should(v1 == v3);
        // should(!(v1 != v1));
        // should(!(v1 != v3));
        // should(v1 <= v1);
        // should(v1 <= v3);
        // should(!(v1 < v1));
        // should(!(v1 < v3));
        // should(v1 >= v1);
        // should(v1 >= v3);
        // should(!(v1 > v1));
        // should(!(v1 > v3));

        // should(v2 != v1);
        // should(v2 != v3);
        // should(!(v2 == v1));
        // should(!(v2 == v3));
        // should(!(v2 <= v1));
        // should(!(v2 <= v3));
        // should(!(v2 < v1));
        // should(!(v2 < v3));
        // should(v2 >= v1);
        // should(v2 >= v3);
        // should(v2 > v1);
        // should(v2 > v3);
    // }

    // void testArithmetic()
    // {
        // typedef FixedPoint16<1> FP1;
        // typedef FixedPoint16<2> FP2;
        // typedef FixedPoint16<7> FP7;
        // typedef FixedPoint16<8> FP8;
        // typedef FixedPoint16<13> FP13;
        // typedef FixedPoint16<15> FP15;
        
        // FP1 t0(0), t1(0.75), t2(0.25);
        // signed char v1 = 1, v2 = 2, v4 = 4, v8 = 8;

        // shouldEqual(FP1(t1) += t1, FP1(1.5));
        // shouldEqual(FP1(t1) -= t1, FP1(0.0));
        // shouldEqual(FP1(t1) -= t2, FP1(0.5));
        // shouldEqual(FP2(t1) *= t1, FP1(9.0 / 16.0));
        // shouldEqual(FP2(t1) /= t2, FP2(3));
        // shouldEqual(FP2(t1) /= t0, NumericTraits<FP2>::max());
        // shouldEqual(FP2(-t1) /= t0, NumericTraits<FP2>::min());
        
        // FP2 res;
        // shouldEqual(add(t1, t1, res), FP2(1.5));
        // shouldEqual(sub(t1, t1, res), FP2(0));
        // shouldEqual(sub(t1, t2, res), FP2(0.5));
        // shouldEqual(mul(t1, t1, res), FP2(9.0 / 16.0));
        // shouldEqual(div(t1, t2, res), FP2(3));
        // shouldEqual(div(t1, t0, res), NumericTraits<FP2>::max());
        // shouldEqual(div(-t1, t0, res), NumericTraits<FP2>::min());

        // shouldEqual(--t1, FP1(-0.25));
        // shouldEqual(t1, FP1(-0.25));
        // shouldEqual(++t1, FP1(0.75));
        // shouldEqual(t1, FP1(0.75));
        // shouldEqual(t1++, FP1(0.75));
        // shouldEqual(t1, FP1(1.75));
        // shouldEqual(t1--, FP1(1.75));
        // shouldEqual(t1, FP1(0.75));

        // shouldEqual((t1 * FP7(v1)).value, 3 << 6);
        // shouldEqual((t2 * FP7(v1)).value, 1 << 6);
        // shouldEqual((-t1 * FP7(v1)).value, -3 << 6);
        // shouldEqual((-t2 * FP7(v1)).value, -1 << 6);
        // shouldEqual((t1 * -FP7(v1)).value, -3 << 6);
        // shouldEqual((t2 * -FP7(v1)).value, -1 << 6);

        // shouldEqual((FixedPoint16<2, FPOverflowSaturate>(t1*FP7(v8)).value), (1 << 15)-1);
        // shouldEqual((FixedPoint16<2, FPOverflowSaturate>(t1*FP7(-v8)).value), -(1 << 15));
        // try { FixedPoint16<2, FPOverflowError>(t1*FP7(v8)); failTest("No exception thrown"); } 
        // catch(PreconditionViolation &) {}
        // try { FixedPoint16<2, FPOverflowError>(t1*FP7(-v8)); failTest("No exception thrown"); } 
        // catch(PreconditionViolation &) {}

        // shouldEqual((FP13(t1 * FP7(v1))).value, 3);
        // shouldEqual((FP13(t2 * FP7(v1))).value, 1);
        // shouldEqual((FP13(-t1 * FP7(v1))).value, -3);
        // shouldEqual((FP13(-t2 * FP7(v1))).value, -1);

        // shouldEqual((t1 * FP7(v4) + t2 * FP7(v8)).value, 5 << 8);
        // shouldEqual((t1 * FP7(v1) + t2 * FP7(v2)).value, 5 << 6);

        // shouldEqual(FP7(6) / FP7(3), FP7(2));
        // shouldEqual(FP7(0.75) / FP7(0.25), FP7(3));
        // shouldEqual(FP7(12) / FP7(48), FP7(0.25));
        // shouldEqual(FP1(0.25) / FP7(2), FP7(0.125));
        // shouldEqual(FP7(10) / FP1(0.25), FP7(40));
        // shouldEqual(FP7(10) / t0, NumericTraits<FP7>::max());
        // shouldEqual(FP7(-10) / t0, NumericTraits<FP7>::min());

        // shouldEqual(floor(t1 * FP7(v1) + t2 * FP7(v2)), 1);
        // shouldEqual(ceil(t1 * FP7(v1) + t2 * FP7(v2)), 2);
        // shouldEqual(round(t1 * FP7(v1) + t2 * FP7(v2)), 1);
        // shouldEqual(floor(t1 * FP7(v4) + t2 * FP7(v8)), 5);
        // shouldEqual(ceil(t1 * FP7(v4) + t2 * FP7(v8)), 5);
        // shouldEqual(round(t1 * FP7(v4) + t2 * FP7(v8)), 5);

        // shouldEqual(floor(t1 * -FP7(v1) - t2 * FP7(v2)), -2);
        // shouldEqual(ceil(t1 * -FP7(v1) - t2 * FP7(v2)), -1);
        // shouldEqual(round(t1 * -FP7(v1) - t2 * FP7(v2)), -1);
        // shouldEqual(floor(t1 * -FP7(v4) - t2 * FP7(v8)), -5);
        // shouldEqual(ceil(t1 * -FP7(v4) - t2 * FP7(v8)), -5);
        // shouldEqual(round(t1 * -FP7(v4) - t2 * FP7(v8)), -5);

        // double d1 = 1.0 / 3.0, d2 = 1.0 / 7.0;
        // FP1 r1(d1), r2(d2);
        // FP2 r3;
        // add(r1, r2, r3);
        // shouldEqual(r3.value, FP2(d1 + d2).value);
        // sub(r1, r2, r3);
        // shouldEqual(r3.value, FP2(d1 - d2).value);
        // mul(r1, r2, r3);
        // shouldEqual(r3.value >> 2, FP2(d1 * d2).value >> 2);

        // shouldEqual(sqrt(FP7(4)).value, 1 << 12);
        // shouldEqual(sqrt(FP8(4)).value, 1 << 12);
        // shouldEqual(hypot(FP8(3), FP8(4)), FP8(5));
        // shouldEqual(hypot(FP8(-3), FP8(-4)), FP8(5));
        // shouldEqual(fixed_point_cast<double>(sqrt(FP7(4))), 2.0);
        // shouldEqual(fixed_point_cast<double>(sqrt(FP2(2.25))), 1.5);
        // shouldEqual(fixed_point_cast<double>(sqrt(FP8(6.25))), 2.5);

        // for(int i = 0; i < 1024; ++i)
        // {
            // FixedPoint16<11> fv1(i, FPNoShift);
            // FixedPoint16<10> fv2(i, FPNoShift);
            // FixedPoint16<9>  fv3(i, FPNoShift);
            // FixedPoint16<8>  fv4(i, FPNoShift);
            // shouldEqual(fixed_point_cast<double>(sqrt(fv1)), floor(sqrt((double)(i << 14))) / 512.0);
            // shouldEqual(fixed_point_cast<double>(sqrt(fv2)), floor(sqrt((double)(i << 15))) / 1024.0);
            // shouldEqual(fixed_point_cast<double>(sqrt(fv3)), floor(sqrt((double)(i << 14))) / 1024.0);
            // shouldEqual(fixed_point_cast<double>(sqrt(fv4)), floor(sqrt((double)(i << 15))) / 2048.0);
        // }

        // shouldEqual(atan2(FP1(0), FP1(1)), FP2(0));
        // shouldEqual(atan2(FP1(0), FP1(-1)), FP2(M_PI));
        // shouldEqual(atan2(FP1(1), FP1(0)), FP2(0.5*M_PI));
        // shouldEqual(atan2(FP1(-1), FP1(0)), FP2(-0.5*M_PI));
        
        // for(int i = -179; i < 180; ++i)
        // {
            // double angle = M_PI*i/180.0;
            // double c = std::cos(angle), s = std::sin(angle);
            // FP2 a = atan2(FP1(s), FP1(c));
            // should(abs(i-fixed_point_cast<double>(a)/M_PI*180.0) < 0.3);
            // a = atan2(FP15(30000.0*s), FP15(30000.0*c));
            // should(abs(i-fixed_point_cast<double>(a)/M_PI*180.0) < 0.3);
        // }
    // }
// };

// struct LinalgTest
// {
    // typedef Matrix<double> Matrix;
    // typedef Matrix::difference_type Shape;

    // unsigned int size, iterations;
    // RandomMT19937 random_;

    // LinalgTest()
    // : size(50),
      // iterations(5),
      // random_(23098349)
    // {}

    // void testOStreamShifting()
    // {
        // Matrix a = random_matrix (size, size);
        // std::ostringstream out;
        // out << a;
        // out << "Testing.." << a << 42 << std::endl;
    // }

    // double random_double ()
    // {
        // double ret = 2.0 * random_.uniform53() - 1.0;
        // return ret;
    // }

    // Matrix random_matrix(unsigned int rows, unsigned int cols)
    // {
        // Matrix ret (rows, cols);
        // for (unsigned int i = 0; i < rows; ++i)
            // for (unsigned int j = 0; j < cols; ++j)
                // ret (i, j) = random_double ();
        // return ret;
    // }

    // Matrix random_symmetric_matrix(unsigned int rows)
    // {
        // Matrix ret (rows, rows);
        // for (unsigned int i = 0; i < rows; ++i)
            // for (unsigned int j = i; j < rows; ++j)
                // ret (j, i) = ret (i, j) = random_double ();
        // return ret;
    // }

    // void testMatrix()
    // {
        // double data[] = {1.0, 5.0,
                         // 3.0, 2.0,
                         // 4.0, 7.0};
        // double tref[] = {1.0, 3.0, 4.0,
                         // 5.0, 2.0, 7.0};
        // double tref2[] = {1.0, 3.0,
                          // 5.0, 2.0};
        // double idref[] = {1.0, 0.0, 0.0,
                          // 0.0, 1.0, 0.0,
                          // 0.0, 0.0, 1.0};
        // std::string sref("1.0000 5.0000 \n3.0000 2.0000 \n4.0000 7.0000 \n");
        // unsigned int r = 3, c = 2;

        // Matrix a(r, c, data), zero(r, c);
        // shouldEqual(a.rowCount(), r);
        // shouldEqual(a.columnCount(), c);
        // shouldEqual(a.elementCount(), r*c);
        // shouldEqual(a.squaredNorm(), 104.0);
        // shouldEqual(a.norm(), std::sqrt(104.0));
        // shouldEqual(a.squaredNorm(), squaredNorm(a));
        // shouldEqual(a.norm(), norm(a));
        // shouldEqual(rowCount(a), r);
        // shouldEqual(columnCount(a), c);

        // for(unsigned int i=0, k=0; i<r; ++i)
            // for(unsigned int j=0; j<c; ++j, ++k)
                // shouldEqual(zero(i,j), 0.0);

        // Matrix one = zero + Matrix(r,c).init(1.0);
        // for(unsigned int i=0, k=0; i<r; ++i)
            // for(unsigned int j=0; j<c; ++j, ++k)
                // shouldEqual(one(i,j), 1.0);

        // std::stringstream s;
        // s << std::setprecision(4) << a;
        // shouldEqual(s.str(), sref);

        // for(unsigned int i=0, k=0; i<r; ++i)
        // {
            // Matrix::view_type ar = a.rowVector(i);
            // shouldEqual(rowCount(ar), 1);
            // shouldEqual(columnCount(ar), c);
            // Matrix::view_type ar1 = rowVector(a, i);
            // shouldEqual(rowCount(ar1), 1);
            // shouldEqual(columnCount(ar1), c);
            // for(unsigned int j=0; j<c; ++j, ++k)
            // {
                // shouldEqual(a(i,j), data[k]);
                // shouldEqual(ar(0, j), data[k]);
                // shouldEqual(ar1(0, j), data[k]);
            // }
        // }

        // Matrix aa(r, c, tref, ColumnMajor);
        // shouldEqual(aa.rowCount(), r);
        // shouldEqual(aa.columnCount(), c);
        // for(unsigned int i=0, k=0; i<r; ++i)
            // for(unsigned int j=0; j<c; ++j, ++k)
                // shouldEqual(aa(i,j), a(i,j));

        // Matrix b = a;
        // shouldEqual(b.rowCount(), r);
        // shouldEqual(b.columnCount(), c);
        // shouldEqualSequence(a.begin(), a.end(), b.begin());

        // b.init(0.0);
        // should(b == zero);

        // Matrix::iterator ib = b.begin();
        // b = a;
        // shouldEqual(ib, b.begin());
        // shouldEqualSequence(a.begin(), a.end(), b.begin());

        // b = 4.0 + a;
        // for(unsigned int i=0, k=0; i<r; ++i)
            // for(unsigned int j=0; j<c; ++j, ++k)
                // shouldEqual(b(i,j), 4.0+data[k]);
        // b = a + 3.0;
        // for(unsigned int i=0, k=0; i<r; ++i)
            // for(unsigned int j=0; j<c; ++j, ++k)
                // shouldEqual(b(i,j), data[k]+3.0);
        // b += 4.0;
        // for(unsigned int i=0, k=0; i<r; ++i)
            // for(unsigned int j=0; j<c; ++j, ++k)
                // shouldEqual(b(i,j), 7.0+data[k]);
        // b += a;
        // for(unsigned int i=0, k=0; i<r; ++i)
            // for(unsigned int j=0; j<c; ++j, ++k)
                // shouldEqual(b(i,j), 7.0+2.0*data[k]);


        // b = 4.0 - a;
        // for(unsigned int i=0, k=0; i<r; ++i)
            // for(unsigned int j=0; j<c; ++j, ++k)
                // shouldEqual(b(i,j), 4.0-data[k]);
        // b = a - 3.0;
        // for(unsigned int i=0, k=0; i<r; ++i)
            // for(unsigned int j=0; j<c; ++j, ++k)
                // shouldEqual(b(i,j), data[k]-3.0);
        // b -= 4.0;
        // for(unsigned int i=0, k=0; i<r; ++i)
            // for(unsigned int j=0; j<c; ++j, ++k)
                // shouldEqual(b(i,j), data[k]-7.0);
        // b -= a;
        // for(unsigned int i=0, k=0; i<r; ++i)
            // for(unsigned int j=0; j<c; ++j, ++k)
                // shouldEqual(b(i,j), -7.0);

        // b = 4.0 * a;
        // for(unsigned int i=0, k=0; i<r; ++i)
            // for(unsigned int j=0; j<c; ++j, ++k)
                // shouldEqual(b(i,j), 4.0*data[k]);
        // b = a * 3.0;
        // for(unsigned int i=0, k=0; i<r; ++i)
            // for(unsigned int j=0; j<c; ++j, ++k)
                // shouldEqual(b(i,j), data[k]*3.0);
        // b *= 4.0;
        // for(unsigned int i=0, k=0; i<r; ++i)
            // for(unsigned int j=0; j<c; ++j, ++k)
                // shouldEqual(b(i,j), data[k]*12.0);
        // b *= a;
        // for(unsigned int i=0, k=0; i<r; ++i)
            // for(unsigned int j=0; j<c; ++j, ++k)
                // shouldEqual(b(i,j), data[k]*data[k]*12.0);

        // b = 4.0 / a;
        // for(unsigned int i=0, k=0; i<r; ++i)
            // for(unsigned int j=0; j<c; ++j, ++k)
                // shouldEqual(b(i,j), 4.0/data[k]);
        // b = a / 3.0;
        // for(unsigned int i=0, k=0; i<r; ++i)
            // for(unsigned int j=0; j<c; ++j, ++k)
                // shouldEqual(b(i,j), data[k] / 3.0);
        // b /= 4.0;
        // for(unsigned int i=0, k=0; i<r; ++i)
            // for(unsigned int j=0; j<c; ++j, ++k)
                // shouldEqual(b(i,j), data[k] / 12.0);
        // b /= a;
        // for(unsigned int i=0, k=0; i<r; ++i)
            // for(unsigned int j=0; j<c; ++j, ++k)
                // shouldEqualTolerance(b(i,j), 1.0 / 12.0, 1e-12);

        // b = a + a;
        // for(unsigned int i=0, k=0; i<r; ++i)
            // for(unsigned int j=0; j<c; ++j, ++k)
                // shouldEqual(b(i,j), 2.0 * data[k]);

        // b = a - a;
        // for(unsigned int i=0; i<r; ++i)
            // for(unsigned int j=0; j<c; ++j)
                // shouldEqual(b(i,j), 0.0);

        // b = -a;
        // for(unsigned int i=0, k=0; i<r; ++i)
            // for(unsigned int j=0; j<c; ++j, ++k)
                // shouldEqual(b(i,j), -data[k]);

        // b = a * pointWise(a);
        // for(unsigned int i=0, k=0; i<r; ++i)
            // for(unsigned int j=0; j<c; ++j, ++k)
                // shouldEqual(b(i,j), data[k] * data[k]);

        // b = a / pointWise(a);
        // for(unsigned int i=0; i<r; ++i)
            // for(unsigned int j=0; j<c; ++j)
                // shouldEqual(b(i,j), 1.0);

        // b = pow(a, 2);
        // for(unsigned int i=0, k=0; i<r; ++i)
            // for(unsigned int j=0; j<c; ++j, ++k)
                // shouldEqual(b(i,j), data[k] * data[k]);

        // b = sqrt(a);
        // for(unsigned int i=0, k=0; i<r; ++i)
            // for(unsigned int j=0; j<c; ++j, ++k)
                // shouldEqual(b(i,j), sqrt(data[k]));

        // b = sq(a);
        // for(unsigned int i=0, k=0; i<r; ++i)
            // for(unsigned int j=0; j<c; ++j, ++k)
                // shouldEqual(b(i,j), sq(data[k]));

        // b = sign(a);
        // for(unsigned int i=0, k=0; i<r; ++i)
            // for(unsigned int j=0; j<c; ++j, ++k)
                // shouldEqual(b(i,j), sign(data[k]));

        // Matrix at = transpose(a);
        // shouldEqual(at.rowCount(), c);
        // shouldEqual(at.columnCount(), r);
        // for(unsigned int i=0, k=0; i<c; ++i)
        // {
            // Matrix::view_type ac = a.columnVector(i);
            // shouldEqual(rowCount(ac), r);
            // shouldEqual(columnCount(ac), 1);
            // Matrix::view_type ac1 = columnVector(a, i);
            // shouldEqual(rowCount(ac1), r);
            // shouldEqual(columnCount(ac1), 1);
            // for(unsigned int j=0; j<r; ++j, ++k)
            // {
                // shouldEqual(at(i,j), tref[k]);
                // shouldEqual(ac(j,0), tref[k]);
                // shouldEqual(ac1(j,0), tref[k]);
            // }
            // shouldEqual(ac, subVector(ac, 0, r));
            // shouldEqual(a.subarray(Shape(1, i), Shape(r-1, i+1)), subVector(ac, 1, r-1));
        // }

        // double sn = squaredNorm(columnVector(a, 0));
        // shouldEqual(sn, 26.0);
        // shouldEqual(sn, dot(columnVector(a, 0), columnVector(a, 0)));
        // shouldEqual(sn, dot(rowVector(at, 0), columnVector(a, 0)));
        // shouldEqual(sn, dot(columnVector(a, 0), rowVector(at, 0)));
        // shouldEqual(sn, dot(rowVector(at, 0), rowVector(at, 0)));
        // shouldEqual(0.0, dot(a.subarray(Shape(0,0), Shape(1,0)), a.subarray(Shape(0,0), Shape(0,1))));
        // shouldEqual(0.0, dot(a.subarray(Shape(0,0), Shape(0,1)), a.subarray(Shape(0,0), Shape(0,1))));
        // shouldEqual(0.0, dot(a.subarray(Shape(0,0), Shape(1,0)), a.subarray(Shape(0,0), Shape(1,0))));
        // shouldEqual(0.0, dot(a.subarray(Shape(0,0), Shape(0,1)), a.subarray(Shape(0,0), Shape(1,0))));

        // Matrix a2(c, c, data);
        // a2 = a2.transpose();
        // for(unsigned int i=0, k=0; i<c; ++i)
            // for(unsigned int j=0; j<c; ++j, ++k)
                // shouldEqual(a2(i,j), tref2[k]);
        
        // shouldEqual(trace(a2), 3.0);

        // Matrix id = identityMatrix<double>(r);
        // shouldEqual(id.rowCount(), r);
        // shouldEqual(id.columnCount(), r);
        // for(unsigned int i=0, k=0; i<r; ++i)
            // for(unsigned int j=0; j<r; ++j, ++k)
                // shouldEqual(id(i,j), idref[k]);

        // shouldEqual(trace(id), 3.0);

        // Matrix d = diagonalMatrix(Matrix(r, 1, data));
        // shouldEqual(d.rowCount(), r);
        // shouldEqual(d.columnCount(), r);
        // for(unsigned int i=0, k=0; i<r; ++i)
            // for(unsigned int j=0; j<r; ++j, ++k)
                // shouldEqual(d(i,j), idref[k]*data[i]);

        // Matrix e(r*c, 1, data);
        // shouldEqual(dot(transpose(e), e), e.squaredNorm());

        // double dc1[] = {1.0, 1.0, 1.0},
               // dc2[] = {1.2, 2.4, 3.6};
        // Matrix c1(3,1, dc1), c2(3,1, dc2);
        // Matrix cr = cross(c1, c2);
        // shouldEqualTolerance(cr(0,0), 1.2, 1e-12);
        // shouldEqualTolerance(cr(1,0), -2.4, 1e-12);
        // shouldEqualTolerance(cr(2,0), 1.2, 1e-12);

        // Matrix f(1, r*c - 1, tref);
        // Matrix g = outer(e, f);
        // shouldEqual(g.rowCount(), e.rowCount());
        // shouldEqual(g.columnCount(), f.columnCount());
        // for(int i=0; i<g.rowCount(); ++i)
            // for(int j=0; j<g.columnCount(); ++j)
                // shouldEqual(g(i,j), data[i]*tref[j]);

        // Matrix g1 = outer(e);
        // shouldEqual(g1.rowCount(), e.rowCount());
        // shouldEqual(g1.columnCount(), e.rowCount());
        // for(int i=0; i<g1.rowCount(); ++i)
            // for(int j=0; j<g1.columnCount(); ++j)
                // shouldEqual(g1(i,j), data[i]*data[j]);

        // Matrix g2 = outer(TinyArray<double, 6>(data));
        // shouldEqual(g2.rowCount(), 6);
        // shouldEqual(g2.columnCount(), 6);
        // for(int i=0; i<g2.rowCount(); ++i)
            // for(int j=0; j<g2.columnCount(); ++j)
                // shouldEqual(g2(i,j), data[i]*data[j]);

        // Matrix h = transpose(a) * a;
        // shouldEqual(h.rowCount(), c);
        // shouldEqual(h.columnCount(), c);
        // for(int i=0; i<(int)c; ++i)
            // for(int j=0; j<(int)c; ++j)
                // shouldEqual(h(i,j), dot(rowVector(at, i), columnVector(a, j)));

        // should(isSymmetric(random_symmetric_matrix(10)));
        // should(!isSymmetric(random_matrix(10, 10)));

        // Matrix tm(2, 2, tref2);
        // TinyArray<double, 2> tv(1.0, 2.0), tvrref(7.0, 9.0), tvlref(11.0, 7.0);
        // shouldEqual(tm * tv, tvrref);
        // shouldEqual(tv * tm, tvlref);

        // Matrix rep = repeatMatrix(a, 2, 4);
        // shouldEqual(rowCount(rep), 2*r);
        // shouldEqual(columnCount(rep), 4*c);

        // for(unsigned int l=0; l<4; ++l)
            // for(unsigned int k=0; k<2; ++k)
                // for(unsigned int j=0; j<c; ++j)
                    // for(unsigned int i=0; i<r; ++i)
                        // shouldEqual(rep(k*r+i, l*c+j), a(i,j));

        // double columnSum[] = {8.0, 14.0};
        // double rowSum[] = {6.0, 5.0, 11.0};
        // Matrix matColumnSum = Matrix(1, 2, columnSum);
        // Matrix matRowSum = Matrix(3, 1, rowSum);
        // shouldEqualSequence(matColumnSum.data(), matColumnSum.data()+2, a.sum(0).data());
        // shouldEqualSequence(matRowSum.data(), matRowSum.data()+3, a.sum(1).data());

        // double columnMean[] = {8/3.0, 14/3.0};
        // double rowMean[] = {3.0, 2.5, 5.5};
        // Matrix matColumnMean = Matrix(1, 2, columnMean);
        // Matrix matRowMean = Matrix(3, 1, rowMean);
        // shouldEqualSequence(matColumnMean.data(), matColumnMean.data()+2, a.mean(0).data());
        // shouldEqualSequence(matRowMean.data(), matRowMean.data()+3, a.mean(1).data());  
    // }

    // void testArgMinMax()
    // {
        // using namespace functor;

        // double data[] = {1.0, 5.0,
                         // 3.0, 2.0,
                        // -2.0, 4.0};
        // unsigned int r = 3, c = 2;
        // Matrix minmax(r, c, data);

        // shouldEqual(argMin(minmax), 2);
        // shouldEqual(argMax(minmax), 3);
        // shouldEqual(argMinIf(minmax, Arg1() > Param(0.0)), 0);
        // shouldEqual(argMinIf(minmax, Arg1() > Param(5.0)), -1);
        // shouldEqual(argMaxIf(minmax, Arg1() < Param(5.0)), 5);
        // shouldEqual(argMaxIf(minmax, Arg1() < Param(-2.0)), -1);
    // }

    // void testColumnAndRowStatistics()
    // {
        // double epsilon = 1e-11;

        // Matrix rowMean(size, 1), columnMean(1, size);
        // Matrix rowStdDev(size, 1), columnStdDev(1, size);
        // Matrix rowNorm(size, 1), columnNorm(1, size);
        // Matrix rowCovariance(size, size), columnCovariance(size, size);

        // for(unsigned int i = 0; i < iterations; ++i)
        // {
            // Matrix a = random_matrix (size, size);

            // rowStatistics(a, rowMean, rowStdDev, rowNorm);
            // columnStatistics(a, columnMean, columnStdDev, columnNorm);

            // for(unsigned int k=0; k<size; ++k)
            // {
                // double rm = 0.0, cm = 0.0, rn = 0.0, cn = 0.0, rs = 0.0, cs = 0.0;
                // for(unsigned int l=0; l<size; ++l)
                // {
                    // rm += a(k, l);
                    // cm += a(l, k);
                    // rn += sq(a(k, l));
                    // cn += sq(a(l, k));
                // }
                // rm /= size;
                // cm /= size;
                // rn = std::sqrt(rn);
                // cn = std::sqrt(cn);

                // shouldEqualTolerance(rm, rowMean(k,0), epsilon);
                // shouldEqualTolerance(cm, columnMean(0,k), epsilon);
                // shouldEqualTolerance(rn, rowNorm(k,0), epsilon);
                // shouldEqualTolerance(cn, columnNorm(0,k), epsilon);

                // for(unsigned int l=0; l<size; ++l)
                // {
                    // rs += sq(a(k, l) - rm);
                    // cs += sq(a(l, k) - cm);
                // }
                // rs = std::sqrt(rs / (size-1));
                // cs = std::sqrt(cs / (size-1));

                // shouldEqualTolerance(rs, rowStdDev(k,0), epsilon);
                // shouldEqualTolerance(cs, columnStdDev(0,k), epsilon);
            // }

            // covarianceMatrixOfRows(a, rowCovariance);
            // covarianceMatrixOfColumns(a, columnCovariance);
            // Matrix rowCovarianceRef(size, size), columnCovarianceRef(size, size);
            // for(unsigned int k=0; k<size; ++k)
            // {
                // for(unsigned int l=0; l<size; ++l)
                // {
                    // for(unsigned int m=0; m<size; ++m)
                    // {
                        // rowCovarianceRef(l, m) += (a(l, k) - rowMean(l, 0)) * (a(m, k) - rowMean(m, 0));
                        // columnCovarianceRef(l, m) += (a(k, l) - columnMean(0, l)) * (a(k, m) - columnMean(0, m));
                    // }
                // }
            // }
            // rowCovarianceRef /= (size-1);
            // columnCovarianceRef /= (size-1);

            // shouldEqualSequenceTolerance(rowCovariance.data(), rowCovariance.data()+size*size, rowCovarianceRef.data(), epsilon);
            // shouldEqualSequenceTolerance(columnCovariance.data(), columnCovariance.data()+size*size, columnCovarianceRef.data(), epsilon);
        // }
    // }

    // void testColumnAndRowPreparation()
    // {
        // using ZeroMean;
        // using UnitVariance;
        // using UnitNorm;
        // using UnitSum;

        // double epsilon = 1e-11;

        // Matrix rowMean(size, 1), columnMean(1, size);
        // Matrix rowStdDev(size, 1), columnStdDev(1, size);
        // Matrix rowNorm(size, 1), columnNorm(1, size);

        // Matrix rowPrepared(size, size), columnPrepared(size, size);
        // Matrix rowMeanPrepared(size, 1), columnMeanPrepared(1, size);
        // Matrix rowStdDevPrepared(size, 1), columnStdDevPrepared(1, size);
        // Matrix rowNormPrepared(size, 1), columnNormPrepared(1, size);
        // Matrix rowOffset(size, 1), columnOffset(1, size);
        // Matrix rowScaling(size, 1), columnScaling(1, size);

        // Matrix zeroRowRef(size,1), zeroColRef(1, size);
        // Matrix oneRowRef(size,1), oneColRef(1, size);
        // oneRowRef.init(1.0);
        // oneColRef.init(1.0);

        // {
            // Matrix a = random_matrix (size, size);

            // columnStatistics(a, columnMean, columnStdDev, columnNorm);

            // prepareColumns(a, columnPrepared, columnOffset, columnScaling, UnitSum);
            // shouldEqualSequence(zeroColRef.data(), zeroColRef.data()+size, columnOffset.data());
            // columnScaling *= columnMean;
            // columnScaling *= size;
            // shouldEqualSequenceTolerance(oneColRef.data(), oneColRef.data()+size, columnScaling.data(), epsilon);
            // columnStatistics(columnPrepared, columnMeanPrepared, columnStdDevPrepared, columnNormPrepared);
            // columnMeanPrepared *= size;
            // shouldEqualSequenceTolerance(oneColRef.data(), oneColRef.data()+size, columnMeanPrepared.data(), epsilon);

            // prepareColumns(a, columnPrepared, columnOffset, columnScaling, ZeroMean);
            // columnStatistics(columnPrepared, columnMeanPrepared, columnStdDevPrepared, columnNormPrepared);
            // shouldEqualSequenceTolerance(zeroColRef.data(), zeroColRef.data()+size, columnMeanPrepared.data(), epsilon);
            // shouldEqualSequenceTolerance(columnStdDev.data(), columnStdDev.data()+size, columnStdDevPrepared.data(), epsilon);

            // Matrix ap = columnPrepared / pointWise(repeatMatrix(columnScaling, size, 1)) + repeatMatrix(columnOffset, size, 1);
            // shouldEqualSequenceTolerance(a.data(), a.data()+size*size, ap.data(), epsilon);

            // prepareColumns(a, columnPrepared, columnOffset, columnScaling, UnitNorm);
            // columnStatistics(columnPrepared, columnMeanPrepared, columnStdDevPrepared, columnNormPrepared);
            // shouldEqualSequenceTolerance(oneColRef.data(), oneColRef.data()+size, columnNormPrepared.data(), epsilon);

            // ap = columnPrepared / pointWise(repeatMatrix(columnScaling, size, 1)) + repeatMatrix(columnOffset, size, 1);
            // shouldEqualSequenceTolerance(a.data(), a.data()+size*size, ap.data(), epsilon);

            // prepareColumns(a, columnPrepared, columnOffset, columnScaling, UnitVariance);
            // columnStatistics(columnPrepared, columnMeanPrepared, columnStdDevPrepared, columnNormPrepared);
            // columnMeanPrepared /= columnScaling;
            // shouldEqualSequenceTolerance(columnMean.data(), columnMean.data()+size, columnMeanPrepared.data(), epsilon);
            // shouldEqualSequenceTolerance(oneColRef.data(), oneColRef.data()+size, columnStdDevPrepared.data(), epsilon);

            // ap = columnPrepared / pointWise(repeatMatrix(columnScaling, size, 1)) + repeatMatrix(columnOffset, size, 1);
            // shouldEqualSequenceTolerance(a.data(), a.data()+size*size, ap.data(), epsilon);

            // prepareColumns(a, columnPrepared, columnOffset, columnScaling, ZeroMean | UnitVariance);
            // columnStatistics(columnPrepared, columnMeanPrepared, columnStdDevPrepared, columnNormPrepared);
            // shouldEqualSequenceTolerance(zeroColRef.data(), zeroColRef.data()+size, columnMeanPrepared.data(), epsilon);
            // shouldEqualSequenceTolerance(oneColRef.data(), oneColRef.data()+size, columnStdDevPrepared.data(), epsilon);

            // ap = columnPrepared / pointWise(repeatMatrix(columnScaling, size, 1)) + repeatMatrix(columnOffset, size, 1);
            // shouldEqualSequenceTolerance(a.data(), a.data()+size*size, ap.data(), epsilon);

            // prepareColumns(a, columnPrepared, columnOffset, columnScaling, ZeroMean | UnitNorm);
            // columnStatistics(columnPrepared, columnMeanPrepared, columnStdDevPrepared, columnNormPrepared);
            // shouldEqualSequenceTolerance(zeroColRef.data(), zeroColRef.data()+size, columnMeanPrepared.data(), epsilon);
            // shouldEqualSequenceTolerance(oneColRef.data(), oneColRef.data()+size, columnNormPrepared.data(), epsilon);

            // ap = columnPrepared / pointWise(repeatMatrix(columnScaling, size, 1)) + repeatMatrix(columnOffset, size, 1);
            // shouldEqualSequenceTolerance(a.data(), a.data()+size*size, ap.data(), epsilon);

            // rowStatistics(a, rowMean, rowStdDev, rowNorm);

            // prepareRows(a, rowPrepared, rowOffset, rowScaling, UnitSum);
            // shouldEqualSequence(zeroRowRef.data(), zeroRowRef.data()+size, rowOffset.data());
            // rowScaling *= rowMean;
            // rowScaling *= size;
            // shouldEqualSequenceTolerance(oneRowRef.data(), oneRowRef.data()+size, rowScaling.data(), epsilon);
            // rowStatistics(rowPrepared, rowMeanPrepared, rowStdDevPrepared, rowNormPrepared);
            // rowMeanPrepared *= size;
            // shouldEqualSequenceTolerance(oneRowRef.data(), oneRowRef.data()+size, rowMeanPrepared.data(), epsilon);

            // prepareRows(a, rowPrepared, rowOffset, rowScaling, ZeroMean);
            // rowStatistics(rowPrepared, rowMeanPrepared, rowStdDevPrepared, rowNormPrepared);
            // shouldEqualSequenceTolerance(zeroRowRef.data(), zeroRowRef.data()+size, rowMeanPrepared.data(), epsilon);
            // shouldEqualSequenceTolerance(rowStdDev.data(), rowStdDev.data()+size, rowStdDevPrepared.data(), epsilon);

            // ap = rowPrepared / pointWise(repeatMatrix(rowScaling, 1, size)) + repeatMatrix(rowOffset, 1, size);
            // shouldEqualSequenceTolerance(a.data(), a.data()+size*size, ap.data(), epsilon);

            // prepareRows(a, rowPrepared, rowOffset, rowScaling, UnitNorm);
            // rowStatistics(rowPrepared, rowMeanPrepared, rowStdDevPrepared, rowNormPrepared);
            // shouldEqualSequenceTolerance(oneRowRef.data(), oneRowRef.data()+size, rowNormPrepared.data(), epsilon);

            // ap = rowPrepared / pointWise(repeatMatrix(rowScaling, 1, size)) + repeatMatrix(rowOffset, 1, size);
            // shouldEqualSequenceTolerance(a.data(), a.data()+size*size, ap.data(), epsilon);

            // prepareRows(a, rowPrepared, rowOffset, rowScaling, UnitVariance);
            // rowStatistics(rowPrepared, rowMeanPrepared, rowStdDevPrepared, rowNormPrepared);
            // rowMeanPrepared /= rowScaling;
            // shouldEqualSequenceTolerance(rowMean.data(), rowMean.data()+size, rowMeanPrepared.data(), epsilon);
            // shouldEqualSequenceTolerance(oneRowRef.data(), oneRowRef.data()+size, rowStdDevPrepared.data(), epsilon);

            // ap = rowPrepared / pointWise(repeatMatrix(rowScaling, 1, size)) + repeatMatrix(rowOffset, 1, size);
            // shouldEqualSequenceTolerance(a.data(), a.data()+size*size, ap.data(), epsilon);

            // prepareRows(a, rowPrepared, rowOffset, rowScaling, ZeroMean | UnitVariance);
            // rowStatistics(rowPrepared, rowMeanPrepared, rowStdDevPrepared, rowNormPrepared);
            // shouldEqualSequenceTolerance(zeroRowRef.data(), zeroRowRef.data()+size, rowMeanPrepared.data(), epsilon);
            // shouldEqualSequenceTolerance(oneRowRef.data(), oneRowRef.data()+size, rowStdDevPrepared.data(), epsilon);

            // ap = rowPrepared / pointWise(repeatMatrix(rowScaling, 1, size)) + repeatMatrix(rowOffset, 1, size);
            // shouldEqualSequenceTolerance(a.data(), a.data()+size*size, ap.data(), epsilon);

            // prepareRows(a, rowPrepared, rowOffset, rowScaling, ZeroMean | UnitNorm);
            // rowStatistics(rowPrepared, rowMeanPrepared, rowStdDevPrepared, rowNormPrepared);
            // shouldEqualSequenceTolerance(zeroRowRef.data(), zeroRowRef.data()+size, rowMeanPrepared.data(), epsilon);
            // shouldEqualSequenceTolerance(oneRowRef.data(), oneRowRef.data()+size, rowNormPrepared.data(), epsilon);

            // ap = rowPrepared / pointWise(repeatMatrix(rowScaling, 1, size)) + repeatMatrix(rowOffset, 1, size);
            // shouldEqualSequenceTolerance(a.data(), a.data()+size*size, ap.data(), epsilon);
        // }

        // {
            // Matrix a(size, size, 2.0), aref(size, size, 1.0/std::sqrt((double)size));

            // prepareColumns(a, columnPrepared, columnOffset, columnScaling, ZeroMean | UnitVariance);
            // shouldEqualSequence(a.data(), a.data()+size*size, columnPrepared.data());

            // prepareColumns(a, columnPrepared, columnOffset, columnScaling, ZeroMean | UnitNorm);
            // shouldEqualSequenceTolerance(aref.data(), aref.data()+size*size, columnPrepared.data(), epsilon);
            // Matrix ap = columnPrepared / pointWise(repeatMatrix(columnScaling, size, 1)) + repeatMatrix(columnOffset, size, 1);
            // shouldEqualSequenceTolerance(a.data(), a.data()+size*size, ap.data(), epsilon);

            // prepareRows(a, rowPrepared, rowOffset, rowScaling, ZeroMean | UnitVariance);
            // shouldEqualSequence(a.data(), a.data()+size*size, rowPrepared.data());

            // prepareRows(a, rowPrepared, rowOffset, rowScaling, ZeroMean | UnitNorm);
            // shouldEqualSequenceTolerance(aref.data(), aref.data()+size*size, rowPrepared.data(), epsilon);
            // ap = rowPrepared / pointWise(repeatMatrix(rowScaling, 1, size)) + repeatMatrix(rowOffset, 1, size);
            // shouldEqualSequenceTolerance(a.data(), a.data()+size*size, ap.data(), epsilon);
        // }
    // }

    // void testCholesky()
    // {
        // double epsilon = 1e-11;
        // Matrix idref = identityMatrix<double>(size);

        // for(unsigned int i = 0; i < iterations; ++i)
        // {
            // Matrix a = random_matrix (size, size);
            // a = transpose(a) * a; // make a symmetric positive definite matrix
            // Matrix l(size, size);
            // choleskyDecomposition (a, l);
            // Matrix ch = l * transpose(l);
            // shouldEqualSequenceTolerance(ch.data(), ch.data()+size*size, a.data(), epsilon);
        // }
    // }

    // void testQR()
    // {
        // double epsilon = 1e-11;
        // Matrix idref = identityMatrix<double>(size);

        // for(unsigned int i = 0; i < iterations; ++i)
        // {
            // Matrix a = random_matrix (size, size);
            // Matrix r(size, size);
            // Matrix q(size, size);
            // qrDecomposition (a, q, r);
            // Matrix id = transpose(q) * q;
            // shouldEqualSequenceTolerance(id.data(), id.data()+size*size, idref.data(), epsilon);
            // Matrix qr = q * r;
            // shouldEqualSequenceTolerance(qr.data(), qr.data()+size*size, a.data(), epsilon);
        // }
    // }

    // void testLinearSolve()
    // {
        // double epsilon = 1e-11;
        // int size = 50;

        // for(unsigned int i = 0; i < iterations; ++i)
        // {
            // Matrix a = random_matrix (size, size);
            // Matrix b = random_matrix (size, 1);
            // Matrix x(size, 1);

            // should(linearSolve (a, b, x, "QR"));
            // Matrix ax = a * x;
            // shouldEqualSequenceTolerance(ax.data(), ax.data()+size, b.data(), epsilon);

            // should(linearSolve(a, b, x, "SVD"));
            // ax = a * x;
            // shouldEqualSequenceTolerance(ax.data(), ax.data()+size, b.data(), epsilon);

            // should(linearSolve(a, b, x, "NE"));
            // ax = a * x;
            // shouldEqualSequenceTolerance(ax.data(), ax.data()+size, b.data(), epsilon);

            // Matrix c = transpose(a) * a; // make a symmetric positive definite matrix
            // Matrix d = transpose(a) * b; 
            // should(linearSolve (c, d, x, "Cholesky"));
            // ax = c * x;
            // shouldEqualSequenceTolerance(ax.data(), ax.data()+size, d.data(), epsilon);
        // }

        // size = 4;
        // Matrix a = random_matrix (size, size);
        // Matrix b = random_matrix (size, 1);
        // Matrix x(size, 1);

        // TinyArray<double, 4> vb(b.data()), vx;
        // should(linearSolve (a, b, x));
        // should(linearSolve (a, vb, vx));
        // shouldEqualSequenceTolerance(x.data(), x.data()+size, vx.data(), epsilon);
    // }

    // void testUnderdetermined()
    // {
        // // test singular matrix
        // Matrix a = identityMatrix<Matrix::value_type> (size);
        // a(0,0) = 0;
        // Matrix b = random_matrix (size, 1);
        // Matrix x(size, 1);
        // should(!linearSolve (a, b, x, "Cholesky"));
        // should(!linearSolve (a, b, x, "QR"));
        // should(!linearSolve (a, b, x, "SVD"));

        // {
            // // square, rank-deficient system (compute minimum norm solution)
            // double mdata[] = {1.0,  3.0,  7.0,
                             // -1.0,  4.0,  4.0,
                              // 1.0, 10.0, 18.0};
            // double rhsdata[] = { 5.0, 2.0, 12.0};
            // double refdata[] = { 0.3850, -0.1103, 0.7066 };

            // Matrix m(3,3,mdata), rhs(3,1,rhsdata), xx(3,1);

            // shouldEqual(linearSolveQR(m, rhs, xx), 2u);
            // shouldEqualSequenceTolerance(refdata, refdata+3, xx.data(), 1e-3);
        // }
        // {
            // // underdetermined, full-rank system (compute minimum norm solution)
            // double mdata[] = {2.0, -3.0, 1.0, -6.0,
                              // 4.0,  1.0, 2.0,  9.0,
                              // 3.0,  1.0, 1.0,  8.0};
            // double rhsdata[] = { -7.0, -7.0, -8.0};
            // double refdata[] = { -3.26666666666667, 3.6, 5.13333333333333, -0.86666666666667 };

            // Matrix m(3,4,mdata), rhs(3,1,rhsdata), xx(4,1);

            // shouldEqual(linearSolveQR(m, rhs, xx), 3u);
            // shouldEqualSequenceTolerance(refdata, refdata+4, xx.data(), 1e-12);
        // }
        // {
            // // underdetermined, rank-deficient, consistent system (compute minimum norm solution)
            // double mdata[] = {1.0,  3.0, 3.0, 2.0,
                              // 2.0,  6.0, 9.0, 5.0,
                             // -1.0, -3.0, 3.0, 0.0};
            // double rhsdata[] = { 1.0, 5.0, 5.0};
            // double refdata[] = { -0.211009, -0.633027, 0.963303, 0.110092 };

            // Matrix m(3,4,mdata), rhs(3,1,rhsdata), xx(4,1);

            // shouldEqual(linearSolveQR(m, rhs, xx), 2u);
            // shouldEqualSequenceTolerance(refdata, refdata+4, xx.data(), 1e-5);
        // }
        // {
            // // underdetermined, rank-deficient, inconsistent system (compute minimum norm least squares solution)
            // double mdata[] = {2.0, 1.0,  7.0, -7.0,
                             // -3.0, 4.0, -5.0, -6.0,
                              // 1.0, 1.0,  4.0, -5.0};
            // double rhsdata[] = { 2.0, 3.0, 2.0};
            // double refdata[] = { -0.0627, 0.1561, -0.0321, -0.3427 };

            // Matrix m(3,4,mdata), rhs(3,1,rhsdata), xx(4,1);

            // shouldEqual(linearSolveQR(m, rhs, xx), 2u);
            // shouldEqualSequenceTolerance(refdata, refdata+4, xx.data(), 1e-3);
        // }
    // }

    // void testOverdetermined()
    // {
        // double epsilon = 1e-11;

        // unsigned int n = 5;
        // unsigned int size = 1000;
        // double noiseStdDev = 0.1;

        // Matrix A(size, n), xs(n,1), xq(n,1), xn(n,1), r(size, 1);

        // for(unsigned int iter=0; iter<iterations; ++iter)
        // {
            // // set up a linear regression problem for a polynomial of degree n
            // Matrix weights = random_matrix (n, 1);
            // Matrix v = random_matrix (size, 1);

            // // init rhs with Gaussian noise with zero mean and noiseStdDev
            // Matrix rhs = 0.5*noiseStdDev*random_matrix (size, 1);
            // for(unsigned int k=1; k<12; ++k)
                // rhs += 0.5*noiseStdDev*random_matrix (size, 1);

            // for(unsigned int k=0; k<size; ++k)
            // {
                // for(unsigned int l=0; l<n; ++l)
                // {
                    // A(k,l) = std::pow(v(k,0), double(l));
                    // rhs(k,0) += weights(l,0)*A(k,l);
                // }
            // }

            // shouldEqual(linearSolve(A, rhs, xs, "SVD"), true);

            // // check that solution is indeed a minimum by
            // // testing for zero derivative of the objective
            // Matrix derivative = abs(transpose(A)*(A*xs - rhs));
            // int absIndex = argMax(derivative);
            // shouldEqualTolerance(derivative(absIndex,0), 0.0, epsilon);

            // shouldEqual(linearSolveQR(A, rhs, xq), n);
            // shouldEqualSequenceTolerance(xs.data(), xs.data()+n, xq.data(), epsilon);

            // shouldEqual(linearSolve(A, rhs, xn, "ne"), true);
            // shouldEqualSequenceTolerance(xs.data(), xs.data()+n, xn.data(), epsilon);
        // }
    // }

    // void testIncrementalLinearSolve()
    // {
        // double epsilon = 1e-11;
        // int size = 50;

        // for(unsigned int i = 0; i < iterations; ++i)
        // {
            // Matrix a = random_matrix (size, size);
            // Matrix b = random_matrix (size, 1);
            // Matrix x(size, 1);

            // should(linearSolve(a, b, x, "QR"));

            // {
                // Matrix r(a), qtb(b), px(size,1), xx(size,1);
                // std::vector<unsigned int> permutation(size);

                // for(int k=0; k<size; ++k)
                // {
                    // // use Givens steps for a change (Householder steps like
                    // //    should(linalg::detail::qrColumnHouseholderStep(k, r, qtb));
                    // // work as well, but are already extensively tested within the QR algorithm)
                    // should(linalg::detail::qrGivensStepImpl(k, r, qtb));
                    // permutation[k] = k;
                // }

                // for(int k=0; k<size; ++k)
                // {
                    // int i = random_.uniformInt(size), j = random_.uniformInt(size);
                    // if(i==j) continue;

                    // linalg::detail::upperTriangularCyclicShiftColumns(i, j, r, qtb, permutation);
                // }
                // should(linalg::linearSolveUpperTriangular(r, qtb, px));
                // linalg::detail::inverseRowPermutation(px, xx, permutation);

                // shouldEqualSequenceTolerance(x.data(), x.data()+size, xx.data(), epsilon);
            // }

            // {
                // Matrix r(a), qtb(b), px(size,1), xx(size,1);
                // std::vector<unsigned int> permutation(size);

                // for(int k=0; k<size; ++k)
                // {
                    // // use Givens steps for a change (Householder steps like
                    // //    should(linalg::detail::qrColumnHouseholderStep(k, r, qtb));
                    // // work as well, but are already extensively tested within the QR algorithm)
                    // should(linalg::detail::qrGivensStepImpl(k, r, qtb));
                    // permutation[k] = k;
                // }

                // for(int k=0; k<size; ++k)
                // {
                    // int i = random_.uniformInt(size), j = random_.uniformInt(size);
                    // linalg::detail::upperTriangularSwapColumns(i, j, r, qtb, permutation);
                // }
                // should(linalg::linearSolveUpperTriangular(r, qtb, px));
                // linalg::detail::inverseRowPermutation(px, xx, permutation);

                // shouldEqualSequenceTolerance(x.data(), x.data()+size, xx.data(), epsilon);
            // }
        // }
    // }

    // void testInverse()
    // {
        // double epsilon = 1e-11;
        // Matrix idref = identityMatrix<double>(size);

        // for(unsigned int i = 0; i < iterations; ++i)
        // {
            // Matrix a = random_matrix (size, size);
            // Matrix id = a * inverse(a);
            // shouldEqualSequenceTolerance(id.data(), id.data()+size*size, idref.data(), epsilon);
            // id = inverse(a) * a;
            // shouldEqualSequenceTolerance(id.data(), id.data()+size*size, idref.data(), epsilon);
            // id = inverse(idref * a) * a; // test inverse(const TemporaryMatrix<T> &v)
            // shouldEqualSequenceTolerance(id.data(), id.data()+size*size, idref.data(), epsilon);
        // }

        // double data[] = { 1.0, 0.0, 0.0, 0.0 };
        // Matrix singular(2, 2, data);
        // try {
            // inverse(singular);
            // failTest("inverse(singular) didn't throw an exception.");
        // }
        // catch(PreconditionViolation & c)
        // {
            // std::string expected("\nPrecondition violation!\ninverse(): matrix is not invertible.");
            // std::string message(c.what());
            // should(0 == expected.compare(message.substr(0,expected.size())));
        // }

        // // test pseudo-inverse
        // double data2[] = { 0.,  1.,  0.,  0.,  0.,
                           // 0.,  0.,  1.,  0.,  2.,
                           // 2.,  0.,  0.,  3.,  0. };
        // double refdata[] = {  0.0, 0.0, 0.15384615384615388,
                              // 1.0, 0.0, 0.0,
                              // 0.0, 0.2, 0.0,
                              // 0.0, 0.0, 0.23076923076923081,
                              // 0.0, 0.4, 0.0 };
            
        // Matrix m(3, 5, data2), piref(5, 3, refdata), pitref(transpose(piref));
        // Matrix pi = inverse(m);
        // shouldEqual(pi.shape(), Shape(5, 3));
        // shouldEqualSequenceTolerance(piref.data(), piref.data()+15, pi.data(), 1e-15);

        // Matrix pit = inverse(transpose(m));
        // shouldEqual(pit.shape(), Shape(3, 5));
        // shouldEqualSequenceTolerance(pitref.data(), pitref.data()+15, pit.data(), 1e-15);
    // }

    // void testSymmetricEigensystem()
    // {
        // double epsilon = 1e-8;
        // Matrix idref = identityMatrix<double>(size);

        // for(unsigned int i = 0; i < iterations; ++i)
        // {
            // Matrix a = random_symmetric_matrix (size);
            // Matrix ew(size, 1);
            // Matrix ev(size, size);
            // should(symmetricEigensystem(a, ew, ev));
            // Matrix id = ev * transpose(ev);
            // shouldEqualSequenceTolerance(id.data(), id.data()+size*size, idref.data(), epsilon);
            // Matrix ae = ev * diagonalMatrix(ew) * transpose(ev);
            // shouldEqualSequenceTolerance(ae.data(), ae.data()+size*size, a.data(), epsilon);
        // }
    // }

    // void testSymmetricEigensystemAnalytic()
    // {
        // double epsilon = 1e-8;

        // int size = 2;
        // for(unsigned int i = 0; i < iterations; ++i)
        // {
            // Matrix a = random_symmetric_matrix (size);
            // Matrix ew(size, 1), ewref(size, 1);
            // Matrix ev(size, size);
            // symmetricEigensystem(a, ewref, ev);
            // symmetric2x2Eigenvalues(
                // a(0,0), a(0,1),
                // a(1,1),
                // &ew(0,0), &ew(1,0));
            // shouldEqualSequenceTolerance(ew.data(), ew.data()+size, ewref.data(), epsilon);
        // }

        // size = 3;
        // for(unsigned int i = 0; i < iterations; ++i)
        // {
            // Matrix a = random_symmetric_matrix (size);
            // Matrix ew(size, 1), ewref(size, 1);
            // Matrix ev(size, size);
            // symmetricEigensystem(a, ewref, ev);
            // symmetric3x3Eigenvalues<double>(
                // a(0,0), a(0,1), a(0,2),
                // a(1,1), a(1,2),
                // a(2,2),
                // &ew(0,0), &ew(1,0), &ew(2,0));
            // shouldEqualSequenceTolerance(ew.data(), ew.data()+size, ewref.data(), epsilon);
        // }
    // }

    // void testNonsymmetricEigensystem()
    // {
        // double epsilon = 1e-8;
        // Matrix idref = identityMatrix<double>(size);

        // for(unsigned int i = 0; i < iterations; ++i)
        // {
            // Matrix a = random_matrix (size, size);
            // Matrix<std::complex<double> > ew(size, 1);
            // Matrix ev(size, size);
            // should(nonsymmetricEigensystem(a, ew, ev));

            // Matrix ewm(size, size);
            // for(unsigned int k = 0; k < size; k++)
            // {
                // ewm(k, k) = ew(k, 0).real();
                // if(ew(k, 0).imag() > 0.0)
                // {
                    // ewm(k, k+1) = ew(k, 0).imag();
                // }
                // else if(ew(k, 0).imag() < 0.0)
                // {
                    // ewm(k, k-1) = ew(k, 0).imag();
                // }
            // }
            // Matrix ae = ev * ewm * inverse(ev);
            // shouldEqualSequenceTolerance(ae.data(), ae.data()+size*size, a.data(), epsilon);
        // }
    // }

    // void testDeterminant()
    // {
        // double ds2[] = {1, 2, 2, 1};
        // double dns2[] = {1, 2, 3, 1};
        // Matrix ms2(Shape(2,2), ds2);
        // Matrix mns2(Shape(2,2), dns2);
        // double eps = 1e-12;
        // shouldEqualTolerance(determinant(ms2), -3.0, eps);
        // shouldEqualTolerance(determinant(mns2), -5.0, eps);
        // shouldEqualTolerance(logDeterminant(transpose(ms2)*ms2), std::log(9.0), eps);
        // shouldEqualTolerance(logDeterminant(transpose(mns2)*mns2), std::log(25.0), eps);

        // double ds3[] = {1, 2, 3, 2, 3, 1, 3, 1, 2};
        // double dns3[] = {1, 2, 3, 5, 3, 1, 3, 1, 2};
        // Matrix ms3(Shape(3,3), ds3);
        // Matrix mns3(Shape(3,3), dns3);
        // shouldEqualTolerance(determinant(ms3), -18.0, eps);
        // shouldEqualTolerance(determinant(mns3), -21.0, eps);
        // shouldEqualTolerance(determinant(transpose(ms3)*ms3, "Cholesky"), 324.0, eps);
        // shouldEqualTolerance(determinant(transpose(mns3)*mns3, "Cholesky"), 441.0, eps);
        // shouldEqualTolerance(logDeterminant(transpose(ms3)*ms3), std::log(324.0), eps);
        // shouldEqualTolerance(logDeterminant(transpose(mns3)*mns3), std::log(441.0), eps);
    // }

    // void testSVD()
    // {
        // unsigned int m = 6, n = 4;
        // Matrix a(m, n);
        // for(unsigned int i1= 0; i1 < m; i1++)
            // for(unsigned int i2= 0; i2 < n; i2++)
                // a(i1, i2)= random_double();
        // Matrix u(m, n);
        // Matrix v(n, n);
        // Matrix S(n, 1);

        // unsigned int rank = singularValueDecomposition(a, u, S, v);
        // shouldEqual(rank, n);

        // double eps = 1e-11;

        // shouldEqualToleranceMessage(norm(a-u*diagonalMatrix(S)*transpose(v)), 0.0, eps, VIGRA_TOLERANCE_MESSAGE);
        // shouldEqualToleranceMessage(norm(identityMatrix<double>(4) - transpose(u)*u), 0.0, eps, VIGRA_TOLERANCE_MESSAGE);
        // shouldEqualToleranceMessage(norm(identityMatrix<double>(4) - transpose(v)*v), 0.0, eps, VIGRA_TOLERANCE_MESSAGE);
        // shouldEqualToleranceMessage(norm(identityMatrix<double>(4) - v*transpose(v)), 0.0, eps, VIGRA_TOLERANCE_MESSAGE);
    // }
// };

// struct RandomTest
// {
    // void testTT800()
    // {
        // const unsigned int n = 50;
        // unsigned int iref[n] = {
            // 3169973338U, 2724982910U,  347012937U, 1735893326U, 2282497071U,
            // 3975116866U,   62755666U,  500522132U,  129776071U, 1978109378U,
            // 4040131704U, 3800592193U, 3057303977U, 1468369496U,  370579849U,
            // 3630178833U,   51910867U,  819270944U,  476180518U,  190380673U,
            // 1370447020U, 1620916304U,  663482756U, 1354889312U, 4000276916U,
             // 868393086U, 1441698743U, 1086138563U, 1899869374U, 3717419747U,
            // 2455034041U, 2617437696U, 1595651084U, 4148285605U, 1860328467U,
             // 928897371U,  263340857U, 4091726170U, 2359987311U, 1669697327U,
            // 1882626857U, 1635656338U,  897501559U, 3233276032U,  373770970U,
            // 2950632840U, 2706386845U, 3294066568U, 3819538748U, 1902519841U };

        // RandomTT800 random;
        // for(unsigned int k=0; k<n; ++k)
            // shouldEqual(random(), iref[k]);

        // double fref[n] = {
              // 0.738067,   0.634460,   0.080795,   0.404169,   0.531435,
              // 0.925529,   0.014611,   0.116537,   0.030216,   0.460564,
              // 0.940666,   0.884894,   0.711834,   0.341881,   0.086282,
              // 0.845217,   0.012086,   0.190751,   0.110869,   0.044326,
              // 0.319082,   0.377399,   0.154479,   0.315460,   0.931387,
              // 0.202189,   0.335672,   0.252886,   0.442348,   0.865529,
              // 0.571607,   0.609420,   0.371516,   0.965848,   0.433141,
              // 0.216276,   0.061314,   0.952679,   0.549477,   0.388757,
              // 0.438333,   0.380831,   0.208966,   0.752806,   0.087025,
              // 0.686998,   0.630130,   0.766960,   0.889306,   0.442965 };
        // RandomTT800 randomf;
        // for(unsigned int k=0; k<n; ++k)
            // should(abs(randomf.uniform() - fref[k]) < 2e-6);

        // RandomTT800 randomr(RandomSeed);
    // }

    // void testMT19937()
    // {
        // const unsigned int n = 20, skip = 960, ilen = 4;
        // unsigned int first[n] = {
             // 956529277U, 3842322136U, 3319553134U, 1843186657U, 2704993644U,
             // 595827513U,  938518626U, 1676224337U, 3221315650U, 1819026461U,
            // 2401778706U, 2494028885U,  767405145U, 1590064561U, 2766888951U,
            // 3951114980U, 2568046436U, 2550998890U, 2642089177U,  568249289U };
        // unsigned int last[n] = {
            // 2396869032U, 1982500200U, 2649478910U,  839934727U, 3814542520U,
             // 918389387U,  995030736U, 2017568170U, 2621335422U, 1020082601U,
              // 24244213U, 2575242697U, 3941971804U,  922591409U, 2851763435U,
            // 2055641408U, 3695291669U, 2040276077U, 4118847636U, 3528766079U };

        // RandomMT19937 random(0xDEADBEEF);
        // for(unsigned int k=0; k<n; ++k)
            // shouldEqual(random(), first[k]);
        // for(unsigned int k=0; k<skip; ++k)
            // random();
        // for(unsigned int k=0; k<n; ++k)
            // shouldEqual(random(), last[k]);

        // for(unsigned int k=0; k<skip; ++k)
            // should(random.uniformInt(31) < 31);

        // random.seed(0xDEADBEEF);
        // for(unsigned int k=0; k<n; ++k)
            // shouldEqual(random(), first[k]);
        // for(unsigned int k=0; k<skip; ++k)
            // random();
        // for(unsigned int k=0; k<n; ++k)
            // shouldEqual(random(), last[k]);

        // unsigned int firsta[n] = {
            // 1067595299U,  955945823U,  477289528U, 4107218783U, 4228976476U,
            // 3344332714U, 3355579695U,  227628506U,  810200273U, 2591290167U,
            // 2560260675U, 3242736208U,  646746669U, 1479517882U, 4245472273U,
            // 1143372638U, 3863670494U, 3221021970U, 1773610557U, 1138697238U };
        // unsigned int lasta[n] = {
             // 123599888U,  472658308U, 1053598179U, 1012713758U, 3481064843U,
            // 3759461013U, 3981457956U, 3830587662U, 1877191791U, 3650996736U,
             // 988064871U, 3515461600U, 4089077232U, 2225147448U, 1249609188U,
            // 2643151863U, 3896204135U, 2416995901U, 1397735321U, 3460025646U };

        // unsigned int init[ilen] = {0x123, 0x234, 0x345, 0x456};
        // RandomMT19937 randoma(init, ilen);
        // for(unsigned int k=0; k<n; ++k)
            // shouldEqual(randoma(), firsta[k]);
        // for(unsigned int k=0; k<skip; ++k)
            // randoma();
        // for(unsigned int k=0; k<n; ++k)
            // shouldEqual(randoma(), lasta[k]);

        // double ref53[n] = {
            // 0.76275444, 0.98670464, 0.27933125, 0.94218739, 0.78842173,
            // 0.92179002, 0.54534773, 0.38107717, 0.65286910, 0.22765212,
            // 0.74557914, 0.54708246, 0.42043117, 0.19189126, 0.70259889,
            // 0.77408120, 0.04605807, 0.69398269, 0.61711170, 0.10133577};
        // for(unsigned int k=0; k<n; ++k)
            // should(abs(randoma.uniform53()-ref53[k]) < 2e-8);

        // randoma.seed(init, ilen);
        // for(unsigned int k=0; k<n; ++k)
            // shouldEqual(randoma(), firsta[k]);
        // for(unsigned int k=0; k<skip; ++k)
            // randoma();
        // for(unsigned int k=0; k<n; ++k)
            // shouldEqual(randoma(), lasta[k]);
    // }

    // void testRandomFunctors()
    // {
        // const unsigned int n = 50;
        // unsigned int iref[n] = {
            // 3169973338U, 2724982910U,  347012937U, 1735893326U, 2282497071U,
            // 3975116866U,   62755666U,  500522132U,  129776071U, 1978109378U,
            // 4040131704U, 3800592193U, 3057303977U, 1468369496U,  370579849U,
            // 3630178833U,   51910867U,  819270944U,  476180518U,  190380673U,
            // 1370447020U, 1620916304U,  663482756U, 1354889312U, 4000276916U,
             // 868393086U, 1441698743U, 1086138563U, 1899869374U, 3717419747U,
            // 2455034041U, 2617437696U, 1595651084U, 4148285605U, 1860328467U,
             // 928897371U,  263340857U, 4091726170U, 2359987311U, 1669697327U,
            // 1882626857U, 1635656338U,  897501559U, 3233276032U,  373770970U,
            // 2950632840U, 2706386845U, 3294066568U, 3819538748U, 1902519841U };
        // double fref[n] = {
              // 0.738067,   0.634460,   0.080795,   0.404169,   0.531435,
              // 0.925529,   0.014611,   0.116537,   0.030216,   0.460564,
              // 0.940666,   0.884894,   0.711834,   0.341881,   0.086282,
              // 0.845217,   0.012086,   0.190751,   0.110869,   0.044326,
              // 0.319082,   0.377399,   0.154479,   0.315460,   0.931387,
              // 0.202189,   0.335672,   0.252886,   0.442348,   0.865529,
              // 0.571607,   0.609420,   0.371516,   0.965848,   0.433141,
              // 0.216276,   0.061314,   0.952679,   0.549477,   0.388757,
              // 0.438333,   0.380831,   0.208966,   0.752806,   0.087025,
              // 0.686998,   0.630130,   0.766960,   0.889306,   0.442965 };
        // double nref[n] = {
            // 1.35298, 0.764158, -0.757076, -0.173069, 0.0586711,
            // 0.794212, -0.483372, -0.0405762, 1.27956, -0.955101,
            // -1.5062, -1.02069, -0.871562, -0.465495, -0.799888,
            // -1.20286, -0.170944, 1.08383, 1.26832, 1.93807,
            // -0.098183, 0.355986, -0.336965, -1.42996, 0.966012,
            // -2.17195, -1.05422, -2.03724, -0.769992, 0.668851,
            // -0.570259, 0.258217, 0.632492, 1.29755, 0.96869,
            // -0.141918, -0.836236, -0.62337, 0.116509, -0.0314471,
            // 0.402451, -1.20504, -0.140861, -0.0765263, 1.06057,
            // 2.57671, 0.0299117, 0.471425, 1.59464, 1.37346};

        // RandomTT800 random1;
        // UniformRandomFunctor<RandomTT800> f1(random1);
        // for(unsigned int k=0; k<n; ++k)
            // should(abs(f1() - fref[k]) < 2e-6);

        // RandomTT800 random2;
        // UniformIntRandomFunctor<RandomTT800> f2(4, 34, random2, true);
        // for(unsigned int k=0; k<n; ++k)
            // shouldEqual(f2(), iref[k] % 31 + 4);

        // RandomTT800 random3;
        // UniformIntRandomFunctor<RandomTT800> f3(random3);
        // for(unsigned int k=0; k<n; ++k)
            // shouldEqual(f3(32), iref[k] % 32);

        // RandomTT800 random4;
        // NormalRandomFunctor<RandomTT800> f4(random4);
        // for(unsigned int k=0; k<n; ++k)
            // shouldEqualTolerance(f4(), nref[k], 1e-5);
    // }
// };


struct MathTestSuite
: public test_suite
{
    MathTestSuite()
    : test_suite("MathTest")
    {
        add( testCase(&MathTest::testSpecialIntegerFunctions));
        add( testCase(&MathTest::testSpecialFunctions));

        // add( testCase(&QuaternionTest::testContents));
        // add( testCase(&QuaternionTest::testStreamIO));
        // add( testCase(&QuaternionTest::testOperators));
        // add( testCase(&QuaternionTest::testRotation));

        // add( testCase(&LinalgTest::testOStreamShifting));
        // add( testCase(&LinalgTest::testMatrix));
        // add( testCase(&LinalgTest::testArgMinMax));
        // add( testCase(&LinalgTest::testColumnAndRowStatistics));
        // add( testCase(&LinalgTest::testColumnAndRowPreparation));
        // add( testCase(&LinalgTest::testCholesky));
        // add( testCase(&LinalgTest::testQR));
        // add( testCase(&LinalgTest::testLinearSolve));
        // add( testCase(&LinalgTest::testUnderdetermined));
        // add( testCase(&LinalgTest::testOverdetermined));
        // add( testCase(&LinalgTest::testIncrementalLinearSolve));
        // add( testCase(&LinalgTest::testInverse));
        // add( testCase(&LinalgTest::testSymmetricEigensystem));
        // add( testCase(&LinalgTest::testNonsymmetricEigensystem));
        // add( testCase(&LinalgTest::testSymmetricEigensystemAnalytic));
        // add( testCase(&LinalgTest::testDeterminant));
        // add( testCase(&LinalgTest::testSVD));

        // add( testCase(&FixedPointTest::testConstruction));
        // add( testCase(&FixedPointTest::testComparison));
        // add( testCase(&FixedPointTest::testArithmetic));
        // add( testCase(&FixedPoint16Test::testConstruction));
        // add( testCase(&FixedPoint16Test::testComparison));
        // add( testCase(&FixedPoint16Test::testArithmetic));

        // add( testCase(&RandomTest::testTT800));
        // add( testCase(&RandomTest::testMT19937));
        // add( testCase(&RandomTest::testRandomFunctors));
    }
};

int main(int argc, char ** argv)
{
  try
  {
    MathTestSuite test;

    int failed = test.run(testsToBeExecuted(argc, argv));

    std::cerr << test.report() << std::endl;

    return (failed != 0);
  }
  catch(std::exception & e)
  {
    std::cerr << "Unexpected exception: " << e.what() << "\n";
    return 1;
  }
}
