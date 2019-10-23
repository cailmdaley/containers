#include <pybindings.h>
#include <mapmaker/mappingutils.h>
#include <mapmaker/pointingutils.h>

namespace bp = boost::python;
BOOST_PYTHON_MODULE(mapmaker)
{
	bp::import("spt3g.core");
	bp::docstring_options docopts(true, true, false);
	pointingutils_pybindings();
	mappingutils_pybindings();
	G3ModuleRegistrator::CallRegistrarsFor("mapmaker");
}


