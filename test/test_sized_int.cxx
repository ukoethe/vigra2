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
#include <type_traits>
#include <vigra2/unittest.hxx>
#include <vigra2/sized_int.hxx>

struct SizedIntTest
{
    SizedIntTest()
    {
    }
    
    void test()
    {
        shouldEqual(sizeof(vigra::int8_t), 1u);
        shouldEqual(sizeof(vigra::int16_t), 2u);
        shouldEqual(sizeof(vigra::int32_t), 4u);
        shouldEqual(sizeof(vigra::int64_t), 8u);
        shouldEqual(sizeof(vigra::uint8_t), 1u);
        shouldEqual(sizeof(vigra::uint16_t), 2u);
        shouldEqual(sizeof(vigra::uint32_t), 4u);
        shouldEqual(sizeof(vigra::uint64_t), 8u);
        
        should(std::is_integral<vigra::int8_t>::value);
        should(std::is_integral<vigra::int16_t>::value);
        should(std::is_integral<vigra::int32_t>::value);
        should(std::is_integral<vigra::int64_t>::value);
        should(std::is_integral<vigra::uint8_t>::value);
        should(std::is_integral<vigra::uint16_t>::value);
        should(std::is_integral<vigra::uint32_t>::value);
        should(std::is_integral<vigra::uint64_t>::value);
    }
};

struct SizedIntTestSuite
: public vigra::test_suite
{
    SizedIntTestSuite()
    : vigra::test_suite("SizedIntTest")
    {
        add( testCase(&SizedIntTest::test));
    }
};

int main(int argc, char ** argv)
{
    SizedIntTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}
