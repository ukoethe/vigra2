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
#include <vigra2/axistags.hxx>

using namespace vigra;

struct AxisTagsTest
{
    AxisTagsTest()
    {
    }
    
    void test()
    {
        using namespace tags;
        shouldEqual(to_string(AxisInfo(axis_x)),  std::string("AxisInfo: 'x' (type: Space)"));
        shouldEqual(to_string(AxisInfo(axis_y)),  std::string("AxisInfo: 'y' (type: Space)"));
        shouldEqual(to_string(AxisInfo(axis_z)),  std::string("AxisInfo: 'z' (type: Space)"));
        shouldEqual(to_string(AxisInfo(axis_t)),  std::string("AxisInfo: 't' (type: Time)"));
        shouldEqual(to_string(AxisInfo(axis_c)),  std::string("AxisInfo: 'c' (type: Channels)"));
        shouldEqual(to_string(AxisInfo(axis_n)),  std::string("AxisInfo: 'n' (type: Space)"));
        shouldEqual(to_string(AxisInfo(axis_e)),  std::string("AxisInfo: 'e' (type: Edge)"));
        shouldEqual(to_string(AxisInfo(axis_fx)), std::string("AxisInfo: 'x' (type: Space Frequency)"));
        shouldEqual(to_string(AxisInfo(axis_fy)), std::string("AxisInfo: 'y' (type: Space Frequency)"));
        shouldEqual(to_string(AxisInfo(axis_fz)), std::string("AxisInfo: 'z' (type: Space Frequency)"));
        shouldEqual(to_string(AxisInfo(axis_ft)), std::string("AxisInfo: 't' (type: Time Frequency)"));
        shouldEqual(to_string(AxisInfo()), std::string("AxisInfo: '?' (type: none)"));
        
        AxisTag sorted[] = { axis_c, axis_n, axis_x, axis_y, axis_z, axis_t, 
                             axis_fx, axis_fy, axis_fz, axis_ft, axis_e, axis_unknown};
        
        for(int i=0; i<11; ++i)
            should(AxisInfo(sorted[i]) < AxisInfo(sorted[i+1]));
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

struct AxisTagsTestSuite
: public vigra::test_suite
{
    AxisTagsTestSuite()
    : vigra::test_suite("AxisTagsTestSuite")
    {
        add( testCase(&AxisTagsTest::test));
    }
};

int main(int argc, char ** argv)
{
    AxisTagsTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}
