from spt3g import core


def add_partner_flagging(flag_map, bolo_partner_ids):
    assert(bolo_partner_ids  != None)
    keys = flag_map.keys()
    for k in keys:
        if not bolo_partner_ids[k] in flag_map:
            flag_map[bolo_partner_ids[k]] =["partner_flagged"]

def get_constant_flagging(timestream_flag_names, timestream_flags,
                          bolometer_flag_names, bolometer_flags,
                          bolo_ids, bolo_partner_ids,
                          ignore_ts_flags = [], invert_ts_flags = [],
                          ignore_bolo_flags = ['has_time_const', 'good_angle_fit', 
                                               'good_xpol_fit'],
                          invert_bolo_flags = ['has_pointing', 'has_polcal'],
                          enforce_partner_good = True):

    ts_flag_ignore_mask = 0
    ignore_masks = [0,0]
    ignore_flag_names = [ignore_ts_flags,ignore_bolo_flags]

    invert_masks = [0,0]
    invert_flag_names = [invert_ts_flags,invert_bolo_flags]

    name_keys = [timestream_flag_names, bolometer_flag_names]
    flag_keys = [timestream_flags, bolometer_flags]

    catty_flag_name = ['FlaggedWithSptPolInsanity']

    for i in range(len(name_keys)):
        for k in name_keys[i]:
            bit_index =  name_keys[i][k]
            if k in ignore_flag_names[i]:
                ignore_masks[i] |= bit_index
            elif k in invert_flag_names[i]:
                invert_masks[i] |= bit_index        
    flags = {}
    for i in range(len(flag_keys)):
        flag_lst = flag_keys[i]
        for j in range(len(flag_lst)):
            flag_val = (flag_lst[j] ^ invert_masks[i]) & (~ignore_masks[i])
            if flag_val:
                flags[bolo_ids[j]] = catty_flag_name

    if enforce_partner_good:
        add_partner_flagging(flags, bolo_partner_ids)
    return flags

def get_scan_flags(scan_ts_flag_lst, bolo_ids, enforce_partner_good = True, 
                   bolo_partner_ids = None):
    flags = {}
    for i in range(len(bolo_ids)):
        if scan_ts_flag_lst[i]:
            flags[bolo_ids[i]] = ['flagged_in_scan']
    if enforce_partner_good:
        add_partner_flagging(flags, bolo_partner_ids)
    return flags

def convert_flag_dic_to_g3map(flags):
    od = core.G3MapVectorString()
    for k in flags:
        od[k] = core.G3VectorString(flags[k])
    return od

