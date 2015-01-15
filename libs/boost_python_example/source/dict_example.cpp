#include <iostream>
#include <boost/python.hpp>
#include <boost/python/object.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/list.hpp>
#include <boost/python/dict.hpp>
#include <boost/numpy.hpp>
#include <boost/scoped_array.hpp>


// using namespace only for example
using namespace boost::python;
namespace np = boost::numpy;

// pass a python object, but this should be a python dictionary.
void pass_dict(boost::python::object const &pydict)
{
    boost::python::extract<boost::python::dict> cppdict_ext(pydict);
    if(!cppdict_ext.check())
    {
        throw std::runtime_error("PassObj::pass_dict: type error: not a python dict.");
    }

    boost::python::dict const cppdict = cppdict_ext();
    boost::python::list const keylist = cppdict.keys();

    // careful with boost name. there already have a conflict.
    int const len = boost::python::len(keylist);
    std::cout << "len(keylist) = " << len << std::endl;
    for(int i = 0; i < len; ++i)
    {
        // operator[] is in python::boost::object
        std::string const key_str = extract<std::string>(str(keylist[i]         ));
        std::string const val_str = extract<std::string>(str(cppdict[keylist[i]]));
        std::cout << "key:[" << key_str << "]->[" << val_str << "]" << std::endl;
    }
}

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
    int const len = boost::python::len(keylist);
    std::cout << "len(keylist) = " << len << std::endl;
    for(int i = 0; i < len; ++i)
    {
        // operator[] is in python::boost::object
        std::string const key_str   = extract<std::string>(str(keylist[i]     ));
        np::ndarray const val_array = extract<np::ndarray>(cppdict[keylist[i]]);
        std::cout << "key:[" << key_str << "]->[" << val_array.shape(0) << "]" << std::endl;
    }
}


BOOST_PYTHON_MODULE(dict_example)
{
    np::initialize();  // have to put this in any module that uses Boost.NumPy
    boost::python::def("pass_dict"   , pass_dict);
    boost::python::def("pass_dict_np", pass_dict_np);
}
