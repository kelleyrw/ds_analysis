#include <boost/python.hpp>
#include <boost/python/object.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/list.hpp>
#include <boost/python/dict.hpp>
#include <iostream>


/// using namespace only for example
using namespace boost::python;

/// pass a python object, but this should be a python dictionary.
/// \param[in] pydict a dictionary
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

BOOST_PYTHON_MODULE(dict_example)
{
//     np::initialize();  // have to put this in any module that uses Boost.NumPy
    boost::python::def("pass_dict", pass_dict);
}
