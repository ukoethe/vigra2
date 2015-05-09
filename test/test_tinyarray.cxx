#include <muli/unittest.hxx>
#include <muli/tinyarray.hxx>
#include <typeinfo>

using namespace muli;

void test()
{
    using Array = TinyArray<int, 2, 3> ;
    std::cerr << Array::static_ndim << " ndim\n";
    std::cerr << Array::static_size << " size\n";
    Array a = {4,5,6,7,8,9};
    std::cerr << a[1] << " a[1]\n";
    std::cerr << a[{0,1}] << " a[{0,1}]\n";
    std::cerr << a[Array::index_type(0,1)] << " a[Array::index_type(0,1)]\n";
    std::cerr << a(0,1) << " a(0,1)\n";
    
    Array::AsType<float> b = a;
    std::cerr << b << " b\n";
    
    a = {0,2,3,4,5,6};
    std::cerr << a << " a\n";
    std::cerr << a.squaredNorm() << " a.squaredNorm()\n";
    std::cerr << a.minimum() << " a.minimum()\n";
    std::cerr << a.maximum() << " a.maximum()\n";
    std::cerr << a.all() << " a.all()\n";
    std::cerr << a.any() << " a.any()\n";
    
    std::cerr << (a+b) << " a+b\n";
    
    
}

struct TinyArrayTestSuite
: public muli::test_suite
{
    TinyArrayTestSuite()
    : muli::test_suite("TinyArrayTest")
    {
        add( testCase(&test));
    }
};



int main(int argc, char ** argv)
{
    TinyArrayTestSuite test;

    int failed = test.run(muli::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}
