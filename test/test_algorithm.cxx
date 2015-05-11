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
#include <vigra2/algorithm.hxx>

using namespace vigra;

struct AlgorithmTest
{
    AlgorithmTest()
    {
    }
    
    void testArgMinMax()
    {
        double data[] = {1.0, 5.0,
                         3.0, 2.0,
                        -2.0, 4.0};
        double *end = data + 6;

        shouldEqual(argMin(data, end), data+4);
        shouldEqual(argMax(data, end), data+1);
        shouldEqual(argMinIf(data, end, [](double x) { return x > 0.0; }), data);
        shouldEqual(argMinIf(data, end, [](double x) { return x > 5.0; }), end);
        shouldEqual(argMaxIf(data, end, [](double x) { return x < 5.0; }), data+5);
        shouldEqual(argMaxIf(data, end, [](double x) { return x < -2.0; }), end);
    }

    void testAlgorithms()
    {
        static const int size = 6;
        int index[size];

        linearSequence(index, index+size);
        int indexref[size] = {0, 1, 2, 3, 4, 5};
        shouldEqualSequence(index, index+size, indexref);

        linearSequence(index, index+size, 5, 5);
        int indexref2[size] = {5, 10, 15, 20, 25, 30};
        shouldEqualSequence(index, index+size, indexref2);

        double data[size] = {1.0, 5.0,
                         3.0, 2.0,
                        -2.0, 4.0};

        indexSort(data, data+size, index, std::greater<double>());
        int sortref[size] = {1, 5, 2, 3, 0, 4};
        shouldEqualSequence(index, index+size, sortref);

        indexSort(data, data+size, index);
        int sortref2[size] = {4, 0, 3, 2, 5, 1};
        shouldEqualSequence(index, index+size, sortref2);

        double res[size];
        applyPermutation(index, index+size, data, res);
        double ref[size] = {-2.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
        shouldEqualSequence(res, res+size, ref);

        int inverse[size];
        inversePermutation(index, index+size, inverse);
        int inverseref[size] = {1, 5, 3, 2, 0, 4};
        shouldEqualSequence(inverse, inverse+size, inverseref);

        applyPermutation(inverse, inverse+size, ref, res);
        shouldEqualSequence(res, res+size, data);
    }

    void testChecksum()
    {
        std::string s("");
        uint32_t crc = checksum(s.c_str(), s.size());
        shouldEqual(crc, 0u);

        s = "hello world";
        crc = checksum(s.c_str(), s.size());
        shouldEqual(crc, 222957957u);

        s = "hallo world";
        crc = checksum(s.c_str(), s.size());
        shouldEqual(crc, 77705727u);

        int split = 5;
        std::string s1 = s.substr(0, split), s2 = s.substr(split);
        crc = checksum(s1.c_str(), s1.size());
        crc = concatenateChecksum(crc, s2.c_str(), s2.size());
        shouldEqual(crc, 77705727u);

        const int size = 446;
        char t[size+1] =  
            "Lorem ipsum dolor sit amet, consectetur adipisicing elit, "
            "sed do eiusmod tempor incididunt ut labore et dolore magna "
            "aliqua. Ut enim ad minim veniam, quis nostrud exercitation "
            "ullamco laboris nisi ut aliquip ex ea commodo consequat. "
            "Duis aute irure dolor in reprehenderit in voluptate velit "
            "esse cillum dolore eu fugiat nulla pariatur. Excepteur "
            "sint occaecat cupidatat non proident, sunt in culpa qui "
            "officia deserunt mollit anim id est laborum.";

        crc = checksum(t, size);
        shouldEqual(crc, 2408722991u);

        for(split = 64; split < 80; ++split) // check alignment
        {
            crc = checksum(t, split);
            crc = concatenateChecksum(crc, t+split, size-split);
            shouldEqual(crc, 2408722991u);
        }
    }
};

struct AlgorithmTestSuite
: public vigra::test_suite
{
    AlgorithmTestSuite()
    : vigra::test_suite("AlgorithmTest")
    {
        add( testCase(&AlgorithmTest::testArgMinMax));
        add( testCase(&AlgorithmTest::testAlgorithms));
        add( testCase(&AlgorithmTest::testChecksum));
    }
};

int main(int argc, char ** argv)
{
    AlgorithmTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}
