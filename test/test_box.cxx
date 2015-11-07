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
#include <vigra2/box.hxx>

using namespace vigra;

template <int DIM>
struct IBoxTest
{
    typedef Box<DIM> IBox;
    typedef typename IBox::EraseResult IBox1;

    IBox rect1_1;
    IBox emptyRect;
    IBox bigRect;

    typedef typename IBox::Vector IPoint;
    typedef typename IBox1::Vector IPoint1;

    IBoxTest()
        : rect1_1(IPoint{1, 1}, IPoint{2, 2}),
          bigRect(IPoint{10, 11})
    {
    }

    void testProperties()
    {
        shouldEqual(rect1_1.shape()[0], 1);
        shouldEqual(rect1_1.shape()[1], 1);
        should(!rect1_1.empty());

        should(emptyRect.ndim() == 0 || emptyRect.shape()[0] <= 0);
        should(emptyRect.ndim() == 0 || emptyRect.shape()[1] <= 0);
        should(emptyRect.empty());

        shouldEqual(bigRect.shape()[0], 10);
        shouldEqual(bigRect.shape()[1], 11);
        should(!bigRect.empty());

        should(rect1_1 != emptyRect);
        should(bigRect != emptyRect);
        should(bigRect != rect1_1);

        bigRect = rect1_1;
        should(bigRect == rect1_1);
    }

    void testContains()
    {
        should(!emptyRect.contains(IPoint{0, 0}));
        should(!emptyRect.contains(IPoint{0, 1}));
        should(!emptyRect.contains(IPoint{0, 2}));
        should(!emptyRect.contains(IPoint{1, 0}));
        should(!emptyRect.contains(IPoint{1, 1}));
        should(!emptyRect.contains(IPoint{1, 2}));
        should(!emptyRect.contains(IPoint{2, 0}));
        should(!emptyRect.contains(IPoint{2, 1}));
        should(!emptyRect.contains(IPoint{2, 2}));

        should( emptyRect.contains(emptyRect));
        should(!emptyRect.contains(rect1_1));
        should(!emptyRect.contains(bigRect));

        should(!rect1_1.contains(IPoint{0, 0}));
        should(!rect1_1.contains(IPoint{0, 1}));
        should(!rect1_1.contains(IPoint{0, 2}));
        should(!rect1_1.contains(IPoint{1, 0}));
        should( rect1_1.contains(IPoint{1, 1}));
        should(!rect1_1.contains(IPoint{1, 2}));
        should(!rect1_1.contains(IPoint{2, 0}));
        should(!rect1_1.contains(IPoint{2, 1}));
        should(!rect1_1.contains(IPoint{2, 2}));

        should( rect1_1.contains(emptyRect));
        should( rect1_1.contains(rect1_1));
        should(!rect1_1.contains(bigRect));

        should(bigRect.contains(IPoint{0, 0}));
        should(bigRect.contains(IPoint{0, 1}));
        should(bigRect.contains(IPoint{0, 2}));
        should(bigRect.contains(IPoint{1, 0}));
        should(bigRect.contains(IPoint{1, 1}));
        should(bigRect.contains(IPoint{1, 2}));
        should(bigRect.contains(IPoint{2, 0}));
        should(bigRect.contains(IPoint{2, 1}));
        should(bigRect.contains(IPoint{2, 2}));

        should( bigRect.contains(emptyRect));
        should( bigRect.contains(rect1_1));
        should( bigRect.contains(bigRect));
    }

    void testIntersection()
    {
        should(!emptyRect.intersects(emptyRect));
        should(!emptyRect.intersects(rect1_1));
        should(!emptyRect.intersects(bigRect));
        should(!rect1_1.intersects(emptyRect));
        should( rect1_1.intersects(rect1_1));
        should( rect1_1.intersects(bigRect));
        should(!bigRect.intersects(emptyRect));
        should( bigRect.intersects(rect1_1));
        should( bigRect.intersects(bigRect));

        should(!bigRect.intersects(IBox(IPoint{3, -3}, IPoint{3, 3})));
        should( bigRect.intersects(IBox(IPoint{3, -3}, IPoint{4, 3})));
        should( bigRect.intersects(IBox(IPoint{3, -3}, IPoint{14, 3})));

        should((rect1_1 & emptyRect).empty());
        should(!(rect1_1 & bigRect).empty());
        should((rect1_1 & bigRect) == rect1_1);
    }

    void testUnion()
    {
        should(!(rect1_1 | emptyRect).empty());
        should((rect1_1 | emptyRect) == rect1_1);
        should((rect1_1 | bigRect) == bigRect);
        rect1_1 |= IPoint{3, 3};
        shouldEqual(rect1_1.lower(), (IPoint{1, 1}));
        shouldEqual(rect1_1.upper(), (IPoint{4, 4}));
    }

    void testShapes()
    {
        shouldEqual(rect1_1.shape(), (IPoint{1, 1}));
        shouldEqual(bigRect.shape(), (IPoint{10, 11}));
        emptyRect.set(IPoint{0, 0}, IPoint{10, 11});
        should(bigRect == emptyRect);
        emptyRect.addShape(IPoint{-4, -7});
        shouldEqual(emptyRect.shape(), (IPoint{6, 4}));
        emptyRect.setShape(bigRect.shape());
        should(bigRect == emptyRect);

        IBox1 rect1D = bigRect.erase(0);
        shouldEqual(rect1D.ndim(), 1);
        shouldEqual(rect1D.lower(), IPoint1{0});
        shouldEqual(rect1D.upper(), IPoint1{11});

        rect1D = bigRect.erase(1);
        shouldEqual(rect1D.ndim(), 1);
        shouldEqual(rect1D.lower(), IPoint1{0});
        shouldEqual(rect1D.upper(), IPoint1{10});

        shouldEqual(rect1D.insert(1, 0, 11), bigRect);
    }

    void testScaling()
    {
        shouldEqual((rect1_1 * 4).shape(), (IPoint{4, 4}));
        shouldEqual((rect1_1 * 4).upper(), (IPoint{8, 8}));
        IBox r2(rect1_1);
        r2 *= 5;
        should(rect1_1 * 5 == r2);
        r2 /= 5;
        should(rect1_1 == r2);
    }
};

struct FBoxTest
{
    typedef Box<2, float> FBox;
    typedef Box<1, float> FBox1;

    FBox rect1_1;
    FBox emptyRect;
    FBox bigRect;

    typedef FBox::Vector FPoint;
    typedef FBox1::Vector FPoint1;

    FBoxTest()
        : rect1_1(FPoint(1, 1), FPoint(2, 2)),
          bigRect(FPoint(10, 11))
    {
    }

    void testProperties()
    {
        shouldEqual(rect1_1.shape()[0], 1);
        shouldEqual(rect1_1.shape()[1], 1);
        should(!rect1_1.empty());

        should(emptyRect.shape()[0] <= 0);
        should(emptyRect.shape()[1] <= 0);
        should(emptyRect.empty());

        shouldEqual(bigRect.shape()[0], 10);
        shouldEqual(bigRect.shape()[1], 11);
        should(!bigRect.empty());

        should(rect1_1 != emptyRect);
        should(bigRect != emptyRect);
        should(bigRect != rect1_1);

        bigRect = rect1_1;
        should(bigRect == rect1_1);
    }

    void testContains()
    {
        should(!emptyRect.contains(FPoint(0, 0)));
        should(!emptyRect.contains(FPoint(0, 1)));
        should(!emptyRect.contains(FPoint(0, 2)));
        should(!emptyRect.contains(FPoint(1, 0)));
        should(!emptyRect.contains(FPoint(1, 1)));
        should(!emptyRect.contains(FPoint(1, 2)));
        should(!emptyRect.contains(FPoint(2, 0)));
        should(!emptyRect.contains(FPoint(2, 1)));
        should(!emptyRect.contains(FPoint(2, 2)));

        should( emptyRect.contains(emptyRect));
        should(!emptyRect.contains(rect1_1));
        should(!emptyRect.contains(bigRect));

        should(!rect1_1.contains(FPoint(0, 0)));
        should(!rect1_1.contains(FPoint(0, 1)));
        should(!rect1_1.contains(FPoint(0, 2)));
        should(!rect1_1.contains(FPoint(1, 0)));
        should( rect1_1.contains(FPoint(1, 1)));
        should( rect1_1.contains(FPoint(1, 2)));
        should(!rect1_1.contains(FPoint(1, 2.1f)));
        should(!rect1_1.contains(FPoint(2, 0)));
        should( rect1_1.contains(FPoint(2, 1)));
        should(!rect1_1.contains(FPoint(2.1f, 1)));
        should( rect1_1.contains(FPoint(2, 2)));

        should( rect1_1.contains(emptyRect));
        should( rect1_1.contains(rect1_1));
        should(!rect1_1.contains(bigRect));

        should(bigRect.contains(FPoint(0, 0)));
        should(bigRect.contains(FPoint(0, 1)));
        should(bigRect.contains(FPoint(0, 2)));
        should(bigRect.contains(FPoint(1, 0)));
        should(bigRect.contains(FPoint(1, 1)));
        should(bigRect.contains(FPoint(1, 2)));
        should(bigRect.contains(FPoint(2, 0)));
        should(bigRect.contains(FPoint(2, 1)));
        should(bigRect.contains(FPoint(2, 2)));

        should( bigRect.contains(emptyRect));
        should( bigRect.contains(rect1_1));
        should( bigRect.contains(bigRect));
    }

    void testIntersection()
    {
        should(!emptyRect.intersects(emptyRect));
        should(!emptyRect.intersects(rect1_1));
        should(!emptyRect.intersects(bigRect));
        should(!rect1_1.intersects(emptyRect));
        should( rect1_1.intersects(rect1_1));
        should( rect1_1.intersects(bigRect));
        should(!bigRect.intersects(emptyRect));
        should( bigRect.intersects(rect1_1));
        should( bigRect.intersects(bigRect));

        should( bigRect.intersects(FBox(FPoint(3, -3), FPoint(4, 3))));
        should( bigRect.intersects(FBox(FPoint(3, -3), FPoint(14, 3))));

        should((rect1_1 & emptyRect).empty());
        should(!(rect1_1 & bigRect).empty());
        should((rect1_1 & bigRect) == rect1_1);
    }

    void testUnion()
    {
        should(!(rect1_1 | emptyRect).empty());
        should((rect1_1 | emptyRect) == rect1_1);
        should((rect1_1 | bigRect) == bigRect);
        rect1_1 |= FPoint(3, 3);
        shouldEqual(rect1_1.lower(), FPoint(1, 1));
        shouldEqual(rect1_1.upper(), FPoint(3, 3));
    }

    void testShapes()
    {
        shouldEqual(rect1_1.shape(), FPoint(1, 1));
        shouldEqual(bigRect.shape(), FPoint(10, 11));
        emptyRect.setLower(FPoint(0, 0));
        emptyRect.setShape(FPoint(10, 11));
        should(bigRect == emptyRect);
        emptyRect.addShape(FPoint(-4, -7));
        shouldEqual(emptyRect.shape(), FPoint(6, 4));
        emptyRect.setShape(bigRect.shape());
        should(bigRect == emptyRect);

        FBox1 rect1D = bigRect.erase(0);
        shouldEqual(rect1D.ndim(), 1);
        shouldEqual(rect1D.lower(), FPoint1{0});
        shouldEqual(rect1D.upper(), FPoint1{11});

        rect1D = bigRect.erase(1);
        shouldEqual(rect1D.ndim(), 1);
        shouldEqual(rect1D.lower(), FPoint1{0});
        shouldEqual(rect1D.upper(), FPoint1{10});

        shouldEqual(rect1D.insert(1, 0, 11), bigRect);
    }

    void testScaling()
    {
        shouldEqual((rect1_1 * 4).shape(), FPoint(4, 4));
        shouldEqual((rect1_1 * 4).upper(), FPoint(8, 8));
        FBox r2(rect1_1);
        r2 *= 5;
        should(rect1_1 * 5 == r2);
        r2 /= 5;
        should(rect1_1 == r2);
    }
};

struct BoxTestSuite
: public test_suite
{
    BoxTestSuite()
    : test_suite("BoxTestSuite")
    {
        add(testCase(&IBoxTest<runtime_size>::testProperties));
        add(testCase(&IBoxTest<runtime_size>::testContains));
        add(testCase(&IBoxTest<runtime_size>::testIntersection));
        add(testCase(&IBoxTest<runtime_size>::testUnion));
        add(testCase(&IBoxTest<runtime_size>::testShapes));
        add(testCase(&IBoxTest<runtime_size>::testScaling));

        add(testCase(&IBoxTest<2>::testProperties));
        add(testCase(&IBoxTest<2>::testContains));
        add(testCase(&IBoxTest<2>::testIntersection));
        add(testCase(&IBoxTest<2>::testUnion));
        add(testCase(&IBoxTest<2>::testShapes));
        add(testCase(&IBoxTest<2>::testScaling));

        add(testCase(&FBoxTest::testProperties));
        add(testCase(&FBoxTest::testContains));
        add(testCase(&FBoxTest::testIntersection));
        add(testCase(&FBoxTest::testUnion));
        add(testCase(&FBoxTest::testShapes));
        add(testCase(&FBoxTest::testScaling));
    }
};


int main(int argc, char ** argv)
{
    BoxTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}
