from spt3g import core, dfmux
from spt3g.calibration import BolometerPropertiesMap
from functools import reduce
import math

@core.usefulfunc
def get_template_groups(boloprops, wiring_map = None, 
                        per_band = True, per_wafer=False, per_squid = False,
                        per_board = False,
                        per_pixel = False, include_keys = False, exclude_empty=False):
    '''
    Return a list of lists of bolometer IDs, sorted according to similarity of
    template images. By default, this just returns a list of the detectors in
    each band. Optionally, it can also split detectors by wafer.
    
    Takes a bolometer properties map as input.
    '''

    if per_squid or per_board:
        assert(not wiring_map is None)

    groups = {}
    for bolo, props in boloprops.iteritems():
        if per_band:
            if exclude_empty and (math.isnan(props.band) or props.band <= 0):
                continue
            dev_lst = [props.band / core.G3Units.GHz if not math.isnan(props.band) else 0]
        else:
            dev_lst = []
        if per_wafer:
            if exclude_empty and not props.wafer_id:
                continue
            dev_lst.append(props.wafer_id)
        if per_board:
            if not bolo in wiring_map:
                continue
            else:
                cmap = wiring_map[bolo]
                dev_lst.append(cmap.board_serial) 
        if per_squid:
            if not bolo in wiring_map:
                continue
            else:
                cmap = wiring_map[bolo]
                dev_lst.append( (cmap.board_serial, cmap.module) ) 
        if per_pixel:
            if exclude_empty and not props.pixel_id:
                continue
            dev_lst.append( props.pixel_id )
        group = str(reduce(lambda x,y: str(x)+'_'+str(y), dev_lst))
        if group == '':
            group = 'All'
        if group not in groups:
            groups[group] = []
        groups[group].append(bolo)
    if include_keys:
        return groups
    else:
        return groups.values()

def get_bolos_template_group(bolo_id, 
                             boloprops, wiring_map = None,
                             per_band = True, per_wafer=False, per_squid = False,
                             per_board = False, 
                             per_pixel = False, exclude_empty=False):
    assert(bolo_id in boloprops)
    new_bprops = BolometerPropertiesMap()
    new_bprops[bolo_id] = boloprops[bolo_id]
    return get_template_groups(new_bprops, wiring_map,
                               per_band, per_wafer, per_squid,
                               per_pixel = per_pixel,
                               per_board = per_board,
                               include_keys = True,
                               exclude_empty=exclude_empty).keys()[0]



def get_template_groups_inverse(boloprops, wiring_map = None, 
                                per_band = True, per_wafer=False, per_squid = False,
                                per_board = False,
                                per_pixel = False, exclude_empty=False):
    tgs =  get_template_groups(boloprops, wiring_map, per_band, per_wafer, per_squid,
                               per_pixel = per_pixel, per_board = per_board,
                               include_keys = True, exclude_empty=exclude_empty)
    tg_inv = {}
    for k, bolos in tgs.items():
        for bid in bolos:
            tg_inv[bid] = k
    return tg_inv


def _is_valid_ident( dev_to_bolos, bolo_to_dev):
    for k in dev_to_bolos.keys():
        for bid in dev_to_bolos[k]:
            if k != bolo_to_dev[bid]:
                return False
    for k in bolo_to_dev.keys():
        if not k in dev_to_bolos[bolo_to_dev[k]]:
            return False
    return True

def _is_valid_hierarchy( low_dev_to_bolos, high_bolo_to_dev, 
                         bolo_props = None, wiring_map = None):
    for k in low_dev_to_bolos.keys():
        for bid_root in low_dev_to_bolos[k]:
            pn = bolo_props[bid_root].physical_name
            if pn != '' and pn[0] != 'F' and (pn[-1] == 'X' or pn[-1] == 'Y'):
                break
        #print("root", bid_root)
        #bid_root = low_dev_to_bolos[k][0]
        dev = high_bolo_to_dev[bid_root]
        for bid in low_dev_to_bolos[k]:
            if high_bolo_to_dev[bid] != dev:
                pn = bolo_props[bid].physical_name
                if pn[0] == 'F' or  pn[-1] == 'H' or pn[-1] == 'J' or pn[-1] == 'D' or pn[0] == 'A':
                    continue
                print('bad hier')
                print(bid_root, bolo_props[bid_root].physical_name)
                print(bid, bolo_props[bid].physical_name)
                print()
                #return False
    return True

def validate_template_groups(bolo_props, wiring_map):
    pixel_to_bolo_map = get_template_groups(
        bolo_props, per_band = False, per_pixel = True, include_keys = True)
    bolo_to_pixel_map = get_template_groups_inverse(
        bolo_props, per_band = False, per_pixel = True)
    
    squid_to_bolo_map = get_template_groups(
        bolo_props, wiring_map, per_band = False, per_squid = True,
        include_keys = True)
    bolo_to_squid_map = get_template_groups_inverse(
        bolo_props, wiring_map, per_band = False, per_squid = True)
    
    wafer_to_bolo_map = get_template_groups(
        bolo_props, wiring_map, per_band = False, per_wafer = True,
        include_keys = True)
    bolo_to_wafer_map = get_template_groups_inverse(
        bolo_props, wiring_map, per_band = False, per_wafer = True)
    #tests identities
    if not _is_valid_ident(pixel_to_bolo_map, bolo_to_pixel_map):
        core.log_fatal("template groups pixel id is not symmetric")
    if not _is_valid_ident(squid_to_bolo_map, bolo_to_squid_map):
        core.log_fatal("template groups squid id is not symmetric")
    if not _is_valid_ident(wafer_to_bolo_map, bolo_to_wafer_map):
        core.log_fatal("template groups wafer id is not symmetric")
    #tests hierarchy
    if not _is_valid_hierarchy( pixel_to_bolo_map, bolo_to_wafer_map, bolo_props, wiring_map):
        core.log_fatal("Detectors on pixel do not share same wafer")
    if not _is_valid_hierarchy( pixel_to_bolo_map, bolo_to_squid_map, bolo_props, wiring_map):
        core.log_fatal("Detectors on pixel do not share same squid")
    if not _is_valid_hierarchy( squid_to_bolo_map, bolo_to_wafer_map, bolo_props, wiring_map):
        core.log_fatal("Detectors on squid do not share same wafer")

