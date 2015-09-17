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
#include <vigra2/loop.hxx>

using namespace vigra;

struct PrintArgument
{
    template <class T, class U>
    void operator()(T const & t, U const & u) const
    {
        std::cerr << t << " = ";
        for(int i=0; i<u.size(); ++i)
            std::cerr << *(u[i].data_) << ", ";
        std::cerr << "\n";
    }
};

struct LoopTest
{
    LoopTest()
    {}
    
    void testArray()
    {
        Shape<> a{1,2,3}, b{1,2,3}, c = a, d = a + b;
        shouldEqual(a, b);
        shouldEqual(a, c);
        c.init(2,4,6);
        shouldEqual(d, c);
        should(all(d));
        
        for(auto i: range(4))
            std::cerr << i << "\n";
    }
    
    void test()
    {
        {
            LoopIndex i, j;
            Loop l{i, j};
            shouldEqual(i.id_, l.indices_[0].id_);
            shouldEqual(j.id_, l.indices_[1].id_);
        }
        {
            LoopIndex i(3), j(2);
            Loop l{i, j};
            l.run(PrintArgument());
        }
        {
            LoopIndex i(-1,-4,-1), j(2);
            Loop l{i, j};
            l.run(PrintArgument());
        }
        {
            MArray<int> a{2,3}, b{3,2};
            LoopIndex i(2), j, _;
            Loop{i, j}
              .add(a(i, j))
              .add(b(j, i))
              // .add(a(_, j))
              .run(PrintArgument())
            ;
        }
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

struct LoopTestSuite
: public vigra::test_suite
{
    LoopTestSuite()
    : vigra::test_suite("LoopTestSuite")
    {
        add( testCase(&LoopTest::test));
        add( testCase(&LoopTest::testArray));
    }
};

int main(int argc, char ** argv)
{
    LoopTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}
