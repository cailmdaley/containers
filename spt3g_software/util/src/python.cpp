#include <pybindings.h>

namespace bp = boost::python;
BOOST_PYTHON_MODULE(util)
{
  bp::docstring_options docopts(true, true, false);
  G3ModuleRegistrator::CallRegistrarsFor("util");
}