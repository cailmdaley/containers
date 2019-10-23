#include <pybindings.h>
#include <G3Frame.h>
#include <G3.h>

namespace bp = boost::python;
BOOST_PYTHON_MODULE(frbutils)
{
  bp::import("spt3g.core");
  bp::import("spt3g.mapmaker");
  bp::import("spt3g.todfilter");
  bp::docstring_options docopts(true, true, false);
  G3ModuleRegistrator::CallRegistrarsFor("frbutils");
}


