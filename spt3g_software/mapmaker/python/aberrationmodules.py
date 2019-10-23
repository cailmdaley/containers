from spt3g.mapmaker import RelativisticAberrationPointing
from spt3g import core
from spt3g.core.util import Rename, Delete

@core.pipesegment_nodoc
def aberrate_alpha_delta(pipe,
        peculiar_ra,
        peculiar_dec,
        peculiar_vel,
        coord_sys,
        alpha_key,
        delta_key):
    """
    Modify the alpha-delta pointings in each scan frame to account for the
    effects of relativistic aberration in the given frame of reference.

    peculiar_ra: Right ascension of the observer's peculiar velocity
    peculiar_dec: Declination of the observer's peculiar velocity
    peculiar_vel: Magnitude of observer's peculiar velocity
    coord_sys: Coordinate system corresponding to the scan
    detector_alpha_key:  Where the detector alpha pointings is stored
    detector_delta_key: Where the detector delta pointing is stored
    """

    tmp_alpha_key="tmp_"+alpha_key
    tmp_delta_key="tmp_"+delta_key

    pipe.Add(RelativisticAberrationPointing,
            peculiar_ra=peculiar_ra,
            peculiar_dec=peculiar_dec,
            peculiar_vel=peculiar_vel,
            coord_sys=coord_sys,
            detector_alpha_key=alpha_key,
            detector_delta_key=delta_key,
            detector_alpha_out_key=tmp_alpha_key,
            detector_delta_out_key=tmp_delta_key)
    pipe.Add(Delete, keys=[alpha_key, delta_key], type=core.G3FrameType.Scan)
    pipe.Add(Rename, keys={tmp_alpha_key:alpha_key, tmp_delta_key:delta_key},
            type=core.G3FrameType.Scan)
