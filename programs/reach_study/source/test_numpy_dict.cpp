#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
//#include "localsolver.h"
#include "NumPyArrayData.hpp" 
#include <boost/numpy.hpp>
#include <boost/scoped_array.hpp>
#include <cmath>
//#include <boost/python/object.hpp>
//#include <boost/python/extract.hpp>
//#include <boost/python/list.hpp>
//#include <boost/python/dict.hpp>


//using namespace localsolver;
using namespace std;
namespace bp = boost::python;
namespace np = boost::numpy;

// pass a python object, but this should be a python dictionary of numpy values.
void pass_dict_np(boost::python::object const &pydict)
{
    boost::python::extract<boost::python::dict> cppdict_ext(pydict);
    if(!cppdict_ext.check())
    {
        throw std::runtime_error("PassObj::pass_dict: type error: not a python dict.");
    }

    boost::python::dict const cppdict = cppdict_ext();
    boost::python::list const keylist = cppdict.keys();

    // careful with boost name. there already have a conflict.
    int const len = bp::len(keylist);
    std::cout << "len(keylist) = " << len << std::endl;
    for(int i = 0; i < len; ++i)
    {
        // operator[] is in python::boost::object
        std::string const key_str   = bp::extract<std::string>(bp::str(keylist[i]     ));
        np::ndarray const val_array = bp::extract<np::ndarray>(cppdict[keylist[i]]);
        NumPyArrayData<double>  val_array_data(val_array); 
        std::cout << "key: " << key_str << std::endl;

        int const num_dim = val_array.get_nd();
        if( num_dim == 1)
        {
            std::cout << "1D array with " 
                      <<  val_array.shape(0) 
                      << " elements\n" 
                      <<  std::endl;
        }
        else if( num_dim == 2)
        {
            std::cout << "2D array with dim (" 
                      <<  val_array.shape(0)  << ", "
                      <<  val_array.shape(1)  << ")\n"
                      <<  std::endl;
        }
        else if( num_dim == 3)
        {
            std::cout << "3D array with dim (" 
                      <<  val_array.shape(0)  << ", "
                      <<  val_array.shape(1)  << ", "
                      <<  val_array.shape(2)  << ")\n"
                      <<  std::endl;
        }
    }
}

BOOST_PYTHON_MODULE(test_numpy_dict) {
    np::initialize();  // have to put this in any module that uses Boost.NumPy
    bp::def("pass_dict_np", pass_dict_np);
}
