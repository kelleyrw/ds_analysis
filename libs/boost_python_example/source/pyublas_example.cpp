#include <pyublas/numpy.hpp>

pyublas::numpy_vector<double> doublify(const pyublas::numpy_vector<double> x)
{
    return 2*x;
}

BOOST_PYTHON_MODULE(pyublas_example)
{
    boost::python::def("doublify", doublify);
}
