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

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <typeinfo>
#include <iostream>
#include <string>
#include <muli/unittest.hxx>
#include <muli/error.hxx>

struct ErrorTest
{
    void testPrecondition()
    {
        try
        {
            muli_precondition(0, "Intentional error");
            failTest("no exception thrown");
        }
        catch(muli::ContractViolation & c)
        {
            std::string expected("\nPrecondition violation!\nIntentional error");
            std::string message(c.what());
            should(0 == expected.compare(message.substr(0,expected.size())));
        }
        try
        {
            muli_assert(0, "Intentional error");
            failTest("no exception thrown");
        }
        catch(muli::ContractViolation & c)
        {
            std::string expected("\nPrecondition violation!\nIntentional error");
            std::string message(c.what());
            should(0 == expected.compare(message.substr(0,expected.size())));
        }
    }

    void testPostcondition()
    {
        try
        {
            muli_postcondition(0, "Intentional error");
            failTest("no exception thrown");
        }
        catch(muli::ContractViolation & c)
        {
            std::string expected("\nPostcondition violation!\nIntentional error");
            std::string message(c.what());
            should(0 == expected.compare(message.substr(0,expected.size())));
        }
    }

    void testInvariant()
    {
        try
        {
            muli_invariant(0, "Intentional error");
            failTest("no exception thrown");
        }
        catch(muli::ContractViolation & c)
        {
            std::string expected("\nInvariant violation!\nIntentional error");
            std::string message(c.what());
            should(0 == expected.compare(message.substr(0,expected.size())));
        }
    }
};

struct ErrorTestSuite
: public muli::test_suite
{
    ErrorTestSuite()
    : muli::test_suite("ErrorTest")
    {
        add( testCase(&ErrorTest::testPrecondition));
        add( testCase(&ErrorTest::testPostcondition));
        add( testCase(&ErrorTest::testInvariant));
    }
};

int main(int argc, char ** argv)
{
    ErrorTestSuite test;

    int failed = test.run(muli::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}
