#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/numpy.hpp>

char const* greet()
{
   return "hello, world";
}

BOOST_PYTHON_MODULE(hello_ext)
{
    using namespace boost::python;
    def("greet", greet);
}
