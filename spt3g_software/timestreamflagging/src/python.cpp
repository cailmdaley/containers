#include <pybindings.h>
#include <timestreamflagging/removeflags.h>

namespace bp = boost::python;
BOOST_PYTHON_MODULE(timestreamflagging)
{
	bp::import("spt3g.core");
	bp::docstring_options docopts(true, true, false);
	G3ModuleRegistrator::CallRegistrarsFor("timestreamflagging");
}
