#include <pybindings.h>
#include <G3Frame.h>
#include <G3.h>

namespace bp = boost::python;
BOOST_PYTHON_MODULE(todfilter)
{
  bp::import("spt3g.core");
  bp::docstring_options docopts(true, true, false);
  G3ModuleRegistrator::CallRegistrarsFor("todfilter");
}


