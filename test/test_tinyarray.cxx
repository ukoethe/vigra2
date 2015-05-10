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

#include <muli/tinyarray.hxx>
#include <muli/algorithm.hxx>

using namespace muli;

static float di[] = {1, 2, 4, 5, 8, 10 };
static float df[] = {1.2f, 2.4f, 3.6f, 4.8f, 8.1f, 9.7f };

template <class BVector, class IVector, class FVector, int SIZE>
struct TinyArrayTest
{
    typedef BVector BV;
    typedef IVector IV;
    typedef FVector FV;

    BV bv0, bv1, bv3;
    IV iv0, iv1, iv3;
    FV fv0, fv1, fv3;

    template <class VECTOR>
    void printVector(VECTOR const & v)
    {
        std::cerr << "(";
        for(unsigned int i=0; i<v.size(); ++i)
            std::cerr << (float)v[i] << ", ";
        std::cerr << ")\n";
    }

    template <class VECTOR, class VALUE>
    bool equalValue(VECTOR const & v, VALUE const & vv)
    {
        for(unsigned int i=0; i<v.size(); ++i)
            if(v[i] != vv)
                return false;
        return true;
    }

    template <class VECTOR1, class VECTOR2>
    bool equalVector(VECTOR1 const & v1, VECTOR2 const & v2)
    {
        for(unsigned int i=0; i<v1.size(); ++i)
            if(v1[i] != v2[i])
                return false;
        return true;
    }

    template <class ITER1, class ITER2>
    bool equalIter(ITER1 i1, ITER1 i1end, ITER2 i2)
    {
        if(i1end - i1 != SIZE)
            return false;
        for(; i1<i1end; ++i1, ++i2)
            if(*i1 != *i2)
                return false;
        return true;
    }

    TinyArrayTest()
    : bv0((unsigned char)0), bv1((unsigned char)1), bv3()
    , iv0(0), iv1(1), iv3()
    , fv0(0.0), fv1(1.0), fv3()
    {
        bv3.init(df, df+SIZE); // check that float inputs are correctly rounded
        iv3.init(di, di+SIZE);
        fv3.init(df, df+SIZE);
    }

    void testConstruction()
    {
        should(bv0.size() == SIZE);
        should(iv0.size() == SIZE);
        should(fv0.size() == SIZE);

        should(equalValue(bv0, 0));
        should(equalValue(iv0, 0));
        should(equalValue(fv0, 0.0f));

        should(equalValue(bv1, 1));
        should(equalValue(iv1, 1));
        should(equalValue(fv1, 1.0f));

        should(equalIter(bv3.begin(), bv3.end(), di));
        should(equalIter(iv3.begin(), iv3.end(), di));
        should(equalIter(fv3.begin(), fv3.end(), df));

        should(!equalVector(bv3, fv3));
        should(!equalVector(iv3, fv3));

        BV bv(fv3);
        should(equalIter(bv3.begin(), bv3.end(), bv.begin()));
        should(equalVector(bv3, bv));

        BV bv4(bv3.begin());
        should(equalIter(bv3.begin(), bv3.end(), bv4.begin()));
        should(equalVector(bv3, bv4));

        BV bv5(bv3.begin(), ReverseCopy);
        should(equalIter(bv3.begin(), bv3.end(),
                         std::reverse_iterator<typename BV::iterator>(bv5.end())));

        FV fv(iv3);
        should(equalIter(iv3.begin(), iv3.end(), fv.begin()));
        should(equalVector(iv3, fv));

        fv = fv3;
        should(equalIter(fv3.begin(), fv3.end(), fv.begin()));
        should(equalVector(fv3, fv));

        fv = bv3;
        should(equalIter(bv3.begin(), bv3.end(), fv.begin()));
        should(equalVector(bv3, fv));

        TinyArray<double, 5> fv5;
        fv5.init(fv3.begin(), fv3.end());
        shouldEqual(fv5[0], fv3[0]);
        shouldEqual(fv5[1], fv3[1]);
        shouldEqual(fv5[2], fv3[2]);
        shouldEqual(fv5[3], SIZE <= 3 ? 0.0 : fv3[3]);
        shouldEqual(fv5[4], SIZE <= 4 ? 0.0 : fv3[4]);
        
        shouldEqual(iv3, (iv3.template subarray<0,SIZE>()));
        shouldEqual(2, (iv3.template subarray<0,2>().size()));
        shouldEqual(iv3[0], (iv3.template subarray<0,2>()[0]));
        shouldEqual(iv3[1], (iv3.template subarray<0,2>()[1]));
        shouldEqual(2, (iv3.template subarray<1,3>().size()));
        shouldEqual(iv3[1], (iv3.template subarray<1,3>()[0]));
        shouldEqual(iv3[2], (iv3.template subarray<1,3>()[1]));
        shouldEqual(1, (iv3.template subarray<1,2>().size()));
        shouldEqual(iv3[1], (iv3.template subarray<1,2>()[0]));
        
        for(int k=0; k<SIZE; ++k)
        {
            IV iv = IV::template unitVector<SIZE>(k);
            shouldEqual(iv[k], 1);
            iv[k] = 0;
            should(!iv.any());
        }

        IV seq = IV::linearSequence(), seq_ref;
        linearSequence(seq_ref.begin(), seq_ref.end());
        shouldEqual(seq, seq_ref);

        seq = IV::linearSequence(2);
        linearSequence(seq_ref.begin(), seq_ref.end(), 2);
        shouldEqual(seq, seq_ref);

        seq = IV::linearSequence(20, -1);
        linearSequence(seq_ref.begin(), seq_ref.end(), 20, -1);
        shouldEqual(seq, seq_ref);

        IV r = reverse(iv3);
        for(int k=0; k<SIZE; ++k)
            shouldEqual(iv3[k], r[SIZE-1-k]);

        shouldEqual(transpose(r, IV::linearSequence(SIZE-1, -1)), iv3);

        typedef TinyArray<typename FV::value_type, SIZE-1> FV1;
        FV1 fv10(fv3.begin());
        shouldEqual(fv10, fv3.dropIndex(SIZE-1));
        FV1 fv11(fv3.begin()+1);
        shouldEqual(fv11, fv3.dropIndex(0));
    }

    void testComparison()
    {
        should(bv0 == bv0);
        should(iv0 == iv0);
        should(fv0 == fv0);
        should(iv0 == bv0);
        should(iv0 == fv0);
        should(fv0 == bv0);

        should(bv3 == bv3);
        should(iv3 == iv3);
        should(fv3 == fv3);
        should(iv3 == bv3);
        should(iv3 != fv3);
        should(fv3 != bv3);

        should(bv0 < bv1);

        should(allLess(bv0, bv1));
        should(!allLess(bv1, bv3));
        should(allGreater(bv1, bv0));
        should(!allGreater(bv3, bv1));
        should(allLessEqual(bv0, bv1));
        should(allLessEqual(bv1, bv3));
        should(!allLessEqual(bv3, bv1));
        should(allGreaterEqual(bv1, bv0));
        should(allGreaterEqual(bv3, bv1));
        should(!allGreaterEqual(bv1, bv3));

        should(closeAtTolerance(fv3, fv3));
        
        should(!bv0.any() && !bv0.all() && bv1.any() && bv1.all());
        should(!iv0.any() && !iv0.all() && iv1.any() && iv1.all());
        should(!fv0.any() && !fv0.all() && fv1.any() && fv1.all());
        IV iv;
        iv = IV(); iv[0] = 1;
        should(iv.any() && !iv.all());
        iv = IV(); iv[1] = 1;
        should(iv.any() && !iv.all());
        iv = IV(); iv[SIZE-1] = 1;
        should(iv.any() && !iv.all());
        iv = IV(1); iv[0] = 0;
        should(iv.any() && !iv.all());
        iv = IV(1); iv[1] = 0;
        should(iv.any() && !iv.all());
        iv = IV(1); iv[SIZE-1] = 0;
        should(iv.any() && !iv.all());
    }

    void testArithmetic()
    {
        IV ivm3 = -iv3;
        FV fvm3 = -fv3;

        int mi[] = { -1, -2, -4, -5, -8, -10};
        float mf[] = { -1.2f, -2.4f, -3.6f, -4.8f, -8.1f, -9.7f };

        should(equalIter(ivm3.begin(), ivm3.end(), mi));
        should(equalIter(fvm3.begin(), fvm3.end(), mf));

        IV iva3 = abs(ivm3);
        FV fva3 = abs(fvm3);
        should(equalVector(iv3, iva3));
        should(equalVector(fv3, fva3));

        int fmi[] = { -2, -3, -4, -5, -9, -10 };
        int fpi[] = { 1, 2, 3, 4, 8, 9 };
        int ri[] = { 1, 2, 4, 5, 8, 10};
        IV ivi3 = floor(fvm3);
        should(equalIter(ivi3.begin(), ivi3.end(), fmi));
        ivi3 = -ceil(fv3);
        should(equalIter(ivi3.begin(), ivi3.end(), fmi));
        ivi3 = round(fv3);
        should(equalIter(ivi3.begin(), ivi3.end(), ri));
        ivi3 = floor(fv3);
        should(equalIter(ivi3.begin(), ivi3.end(), fpi));
        ivi3 = -ceil(fvm3);
        should(equalIter(ivi3.begin(), ivi3.end(), fpi));
        ivi3 = -round(fvm3);
        should(equalIter(ivi3.begin(), ivi3.end(), ri));

        shouldEqual(clipLower(iv3), iv3);
        shouldEqual(clipLower(iv3, 11), IV(11));
        shouldEqual(clipUpper(iv3, 0), IV(0));
        shouldEqual(clipUpper(iv3, 11), iv3);
        shouldEqual(clip(iv3, 0, 11), iv3);
        shouldEqual(clip(iv3, 11, 12), IV(11));
        shouldEqual(clip(iv3, -1, 0), IV(0));
        shouldEqual(clip(iv3, IV(0), IV(11)), iv3);
        shouldEqual(clip(iv3, IV(11), IV(12)), IV(11));
        shouldEqual(clip(iv3, IV(-1), IV(0)), IV(0));

        should(squaredNorm(bv1) == SIZE);
        should(squaredNorm(iv1) == SIZE);
        should(squaredNorm(fv1) == (float)SIZE);

        float expectedSM = 1.2f*1.2f + 2.4f*2.4f + 3.6f*3.6f;
        if(SIZE == 6)
            expectedSM += 4.8f*4.8f + 8.1f*8.1f + 9.7f*9.7f;
        shouldEqualTolerance(squaredNorm(fv3), expectedSM, 1e-7f);

        shouldEqual(dot(bv3, bv3), squaredNorm(bv3));
        shouldEqual(dot(iv3, bv3), squaredNorm(iv3));
        shouldEqual(dot(fv3, fv3), squaredNorm(fv3));

        // shouldEqual(squaredNorm(bv3), bv3.squaredNorm());
        // shouldEqual(squaredNorm(iv3), iv3.squaredNorm());
        // shouldEqual(squaredNorm(fv3), fv3.squaredNorm());

        // shouldEqual(norm(bv3), bv3.norm());
        // shouldEqual(norm(iv3), iv3.norm());
        // shouldEqual(norm(fv3), fv3.norm());

        TinyArray<IV, 3> ivv(iv3, iv3, iv3);
        shouldEqual(squaredNorm(ivv), 3*squaredNorm(iv3));
        shouldEqual(norm(ivv), sqrt(3.0*squaredNorm(iv3)));

        shouldEqualTolerance(sqrt(dot(bv3, bv3)), norm(bv3), 0.0);
        shouldEqualTolerance(sqrt(dot(iv3, bv3)), norm(iv3), 0.0);
        shouldEqualTolerance(sqrt(dot(fv3, fv3)), norm(fv3), 0.0f);

        BV bv = bv3;
        bv[2] = 200;
        int expectedSM2 = 40005;
        if(SIZE == 6)
            expectedSM2 += 189;
        should(dot(bv, bv) == expectedSM2);
        should(squaredNorm(bv) == expectedSM2);

        should(equalVector(bv0 + 1.0, fv1));
        should(equalVector(1.0 + bv0, fv1));
        should(equalVector(bv1 - 1.0, fv0));
        should(equalVector(1.0 - bv1, fv0));
        should(equalVector(bv3 - iv3, bv0));
        should(equalVector(fv3 - fv3, fv0));
        BV bvp = (bv3 + bv3)*0.5;
        FV fvp = (fv3 + fv3)*0.5;
        should(equalVector(bvp, bv3));
        should(equalVector(fvp, fv3));
        bvp = 2.0*bv3 - bv3;
        fvp = 2.0*fv3 - fv3;
        should(equalVector(bvp, bv3));
        should(equalVector(fvp, fv3));

        IV ivp = bv + bv;
        int ip1[] = {2, 4, 400, 10, 16, 20};
        should(equalIter(ivp.begin(), ivp.end(), ip1));
        should(equalVector(bv0 - iv1, -iv1));

        bvp = bv3 / 2.0;
        fvp = bv3 / 2.0;
        int ip[] = {1, 1, 2, 3, 4, 5}; // half-integers are rounded upwards
        float fp[] = {0.5, 1.0, 2.0, 2.5, 4.0, 5.0};
        should(equalIter(bvp.begin(), bvp.end(), ip));
        should(equalIter(fvp.begin(), fvp.end(), fp));
        fvp = fv3 / 2.0;
        float fp1[] = {0.6f, 1.2f, 1.8f, 2.4f, 4.05f, 4.85f};
        should(equalIter(fvp.begin(), fvp.end(), fp1));
        shouldEqual(2.0 / fv1, 2.0 * fv1);        
        float fp2[] = {1.0f, 0.5f, 0.25f, 0.2f, 0.125f, 0.1f};
        fvp = 1.0 / bv3;
        should(equalIter(fvp.begin(), fvp.end(), fp2));

        int ivsq[] = { 1, 4, 16, 25, 64, 100 };
        ivp = iv3*iv3;
        should(equalIter(ivp.begin(), ivp.end(), ivsq));
        shouldEqual(iv3 * iv1, iv3);
        shouldEqual(iv0 * iv3, iv0);
        shouldEqual(iv3 / iv3, iv1);
        shouldEqual(iv3 % iv3, iv0);
        shouldEqual(iv3 % (iv3+iv1), iv3);

        float minRef[] = { 1.0f, 2.0f, 3.6f, 4.8f, 8.0f, 9.7f };
        auto minRes = min(iv3, fv3);
        shouldEqualSequence(minRef, minRef+SIZE, minRes.cbegin());
        IV ivmin = floor(fv3);
        ivmin[1] = 3;
        int minRef2[] = { 1, 2, 3, 4, 8, 9 };
        auto minRes2 = min(iv3, ivmin);
        shouldEqualSequence(minRef2, minRef2+SIZE, minRes2.cbegin());
        shouldEqual(min(iv3), di[0]);
        shouldEqual(min(fv3), df[0]);
        shouldEqual(max(iv3), di[SIZE-1]);
        shouldEqual(max(fv3), df[SIZE-1]);

        float maxRef[] = { 1.2f, 2.4f, 4.0f, 5.0f, 8.1f, 10.0f };
        shouldEqualSequence(maxRef, maxRef+SIZE, max(iv3, fv3).begin());
        IV ivmax = floor(fv3);
        ivmax[1] = 3;
        int maxRef2[] = { 1, 3, 4, 5, 8, 10 };
        shouldEqualSequence(maxRef2, maxRef2+SIZE, max(iv3, ivmax).begin());
        
        shouldEqual(sum(iv3), SIZE == 3 ? 7 : 30);
        shouldEqual(sum(fv3), SIZE == 3 ? 7.2f : 29.8f);
        shouldEqual(prod(iv3), SIZE == 3 ? 8 : 3200);
        shouldEqual(prod(fv3), SIZE == 3 ? 10.368f : 3910.15f);

        float cumsumRef[] = {1.2f, 3.6f, 7.2f, 12.0f, 20.1f, 29.8f };
        shouldEqualSequenceTolerance(cumsumRef, cumsumRef+3, cumsum(fv3).begin(), 1e-6);
        float cumprodRef[] = {1.2f, 2.88f, 10.368f, 49.7664f, 403.108f, 3910.15f };
        shouldEqualSequenceTolerance(cumprodRef, cumprodRef+3, cumprod(fv3).begin(), 1e-6);
    }

    void testCross()
    {
        shouldEqual(cross(bv3, bv3), IV(0));
        shouldEqual(cross(iv3, bv3), IV(0));
        shouldEqualTolerance(cross(fv3, fv3), FV(0.0), FV(1e-6f));

        FV cr = cross(fv1, fv3);
        shouldEqualTolerance(cr[0], 1.2, 1e-6f);
        shouldEqualTolerance(cr[1], -2.4, 1e-6f);
        shouldEqualTolerance(cr[2], 1.2, 1e-6f);
    }

    void testOStreamShifting()
    {
        std::ostringstream out;
        out << iv3;
        std::string expected("{1, 2, 4}");
        shouldEqual(expected, out.str());
        out << "Testing.." << fv3 << 42;
        out << bv3 << std::endl;
    }
    
    void test()
    {
        using Array1 = TinyArray<int, 6> ;
        using Array = TinyArray<int, 2, 3> ;
        std::cerr << Array::static_ndim << " ndim\n";
        std::cerr << Array::static_size << " size\n";
        Array a = {4,5,6,7,8,9};
        std::cerr << a[1] << " a[1]\n";
        std::cerr << a[{0,1}] << " a[{0,1}]\n";
        std::cerr << a[Array::index_type(0,1)] << " a[Array::index_type(0,1)]\n";
        std::cerr << a(0,1) << " a(0,1)\n";
        
        std::cerr << TinySymmetricView<int, 3>(a.data()).shape() << " symshape\n";
        std::cerr << TinySymmetricView<int, 3>(a.data()) << " symview\n";
        
        Array::AsType<float> b = a;
        std::cerr << b << " b\n";
        
        a = {0,2,3,4,5,6};
        std::cerr << a << " a\n";
        std::cerr << squaredNorm(a) << " squaredNorm(a)\n";
        std::cerr << norm(a) << " norm(a)\n";
        std::cerr << a.minimum() << " a.minimum()\n";
        std::cerr << a.maximum() << " a.maximum()\n";
        std::cerr << a.all() << " a.all()\n";
        std::cerr << a.any() << " a.any()\n";
        
        std::cerr << (a+b) << " a+b\n";
        std::cerr << (a+4) << " a+4\n";
        std::cerr << (4+a) << " 4+a\n";
        std::cerr << (-a) << " -a\n";
        std::cerr << sqrt(a) << " sqrt(a)\n";
        std::cerr << sum(a) << " sum(a)\n";
        std::cerr << mean(a) << " mean(a)\n";
        std::cerr << cumsum(a) << " cumsum(a)\n";
        
        a.swap(b);
        std::cerr << b << " b\n";
        std::cerr << a << " a\n";
        std::cerr << (a==a) << " a==a\n";
        std::cerr << (a==b) << " a==b\n";
        std::cerr << isZero(a) << " isZero(a)\n";
        std::cerr << closeAtTolerance(b, b) << " closeAtTolerance(b, b)\n";
        std::cerr << min(a,b) << " min(a,b)\n";
        std::cerr << max(a,b) << " max(a,b)\n";
        
        Array c(a.begin(), ReverseCopy);
        std::cerr << c << " c\n";
        std::cerr << reverse(c) << " reverse(c)\n";
        
        Array1 d = {1,2,3,4,5,6};
        std::cerr << d << " d\n";
        std::cerr << d.subarray<1,4>() << " d.subarray<1,4>()\n";
        std::cerr << d.dropIndex(2) << " d.dropIndex(2)\n";
        std::cerr << dot(d, d) << " dot(d, d)\n";
        
        std::cerr << Array::tinyEye<3>() << " tinyEye<3>()\n";
        std::cerr << Array::unitVector<3>(1) << " Array::unitVector<3>(1)\n";
        std::cerr << Array::linearSequence() << " Array::linearSequence()\n";
        
        TinyArray<int, 3> f = {1,2,3};
        std::cerr << cross(f,f) << " cross(f,f)\n";
        
        TinyArray<float, 3> g = {1.f,2.f,3.f};
        std::cerr << cross(g,g) << " cross(g,g)\n";
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

struct TinyArrayTestSuite
: public muli::test_suite
{
    typedef TinyArrayTest<TinyArray<unsigned char, 3>,
                          TinyArray<int, 3>,
                          TinyArray<float, 3>, 3> Tests;
    TinyArrayTestSuite()
    : muli::test_suite("TinyArrayTest")
    {
        add( testCase(&Tests::test));
        add( testCase(&Tests::testConstruction));
        add( testCase(&Tests::testComparison));
        add( testCase(&Tests::testArithmetic));
        add( testCase(&Tests::testCross));
        add( testCase(&Tests::testOStreamShifting));
    }
};

int main(int argc, char ** argv)
{
    TinyArrayTestSuite test;

    int failed = test.run(muli::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}
