"""
Various utilities for coadding maps
"""
import os
from glob import glob
import argparse
from spt3g import core, coordinateutils, util
from fnmatch import fnmatch

def coadd_maps(maps_in, map_id=None, output_file=None):
    """
    Coadd `maps_in` into a single map.
    
    Parameters
    ----------
    maps_in: str or list of str
        A filepath, directory, or lists of either that
        points to maps stored in .g3(.gz) files.
        
    map_id: str or list of str
        Add maps that have an Id key matching the pattern, or one of the 
        patterns in the list. Maps matching separate patterns are added to
        separate coadds. Understands Unix shell-style wildcards (e.g. *, ?).
        
    output_file: str
        If specified, save the output map to this path.
        
    Returns
    -------
    G3Frame containing the coadded map
    """
    if not isinstance(maps_in, list):
        maps_in = [maps_in]
    if not all([isinstance(mp, str) for mp in maps_in]):
        raise TypeError("All inputs must be strings")
        
    if map_id is not None:
        if not isinstance(map_id, list):
            map_id = [map_id]
        if not all([isinstance(mid, str) for mid in map_id]):
            raise TypeError("All inputs must be strings")


    maps_to_add = []
    for pth in maps_in:
        if os.path.isdir(pth):
            now_maps = glob(os.path.join(pth,'*.g3*'))
            maps_to_add += now_maps
        elif os.path.isfile(pth):
            maps_to_add.append(pth)
        else:
            raise OSError(
                "%s is not an existing file or directory."%pth)
    
    map_out = _GrabMapFrame()
    pipe = core.G3Pipeline()
    pipe.Add(core.G3Reader, filename=maps_to_add)
    pipe.Add(lambda frame: frame.type == core.G3FrameType.Map)
    if not map_id:
        pipe.Add(util.framecombiner.MapFrameCombiner, fr_id=None)
    else:
        for id in map_id:
            pipe.Add(util.framecombiner.MapFrameCombiner, fr_id=id)
    pipe.Add(lambda frame: 'Id' not in frame or \
             fnmatch(frame['Id'], 'combined_*'))
    if output_file is not None:
        if not os.path.exists(os.path.dirname(output_file)):
            os.makedirs(os.path.dirname(output_file))
        pipe.Add(core.G3Writer, filename=output_file)    
    pipe.Add(map_out)
    pipe.Run()
    
    return map_out.map
    
class _GrabMapFrame():
    def __init__(self):
        self.map = None
    def __call__(self, frame):
        if frame.type == core.G3FrameType.Map:
            self.map = frame
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('maps_in', nargs='+')
    parser.add_argument('--output-file', '-o', default='map_coadd.g3',
                       help='output filename')
    parser.add_argument('--map-id', nargs='+', default=None,
                        help='Strings used to filter map ids.')
    args = parser.parse_args()
    
    coadd_maps(args.maps_in, map_id=args.map_id, output_file=args.output_file)
