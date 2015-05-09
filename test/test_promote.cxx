#include <muli/unittest.hxx>
#include <muli/promote.hxx>

void test()
{
    should((std::is_same<unsigned int, muli::Promote<unsigned int, short> >::value));
    should((std::is_same<int, muli::Promote<unsigned char> >::value));
}

struct PromoteTestSuite
: public muli::test_suite
{
    PromoteTestSuite()
    : muli::test_suite("PromoteTest")
    {
        add( testCase(&test));
    }
};



int main(int argc, char ** argv)
{
    PromoteTestSuite test;

    int failed = test.run(muli::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}
