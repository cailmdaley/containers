import os
import sys
import glob

# If running in jupyter-notebook:
#sys.path.append("/home/joshuasobrin/spt3g/spt3g_software/build")

from spt3g import core
from spt3g.std_processing import combining_maps

# Grab all files in 150 band
s1_maps_all = sorted(glob.glob("/spt/data/onlinemaps/ra0hdec-44.75/*_150GHz_tonly.g3.gz"))
s2_maps_all = sorted(glob.glob("/spt/data/onlinemaps/ra0hdec-52.25/*_150GHz_tonly.g3.gz"))
s3_maps_all = sorted(glob.glob("/spt/data/onlinemaps/ra0hdec-59.75/*_150GHz_tonly.g3.gz"))
s4_maps_all = sorted(glob.glob("/spt/data/onlinemaps/ra0hdec-67.25/*_150GHz_tonly.g3.gz"))

exclusion_list = ['70008378', '70019142', '70095756', '70234484', '70307740', '70327253', '70344984', '70366963', '70403820', '70412736', '70482836', '70553740', '70613613', '70748993', '71131345', '71718785', '71764820', '71825113', '71853745', '71899844', '72243773', '72304242', '72313155', '72992268', '73102483', '73275784', '73582657', '74554405', '74563201', '74677973', '74698219', '74707028', '75301520', '75850013', '75914072', '75973289', '75982072', '76914614', '76957677', '76966577', '76975318', '76986953', '76995691', '77058739', '77826016', '77834811', '78183108', '78237672', '78246410', '78255309', '78284577', '78312893', '78529740', '78550148', '78601914', '78816157', '79545645', '79554557', '79574967', '79609421', '79626984', '79891185', '80725484', '80797942', '82301200', '82310114', '82335631', '82344369', '82353107', '82364506', '82795786', '82837159', '82948115', '82977151', '82985889', '82994788', '83006197', '83015021', '83023892', '83120993', '83129776', '83138559', '83147367', '83158827', '83167741', '83193111', '83201849', '83210587', '83290222', '83298960', '83591401', '83602865', '83708598', '83864978']

# Filter maps by date-range and exclusion-list
def filter_input_maps(input_list, starting_ob=70000000, ending_ob=84066968):
    
    filtered_list = []
    
    for obs in input_list:
        name = os.path.basename(obs)
        observation = int(name.split("_")[0])
        if starting_ob <= observation <= ending_ob:
            if str(observation) not in exclusion_list:
                filtered_list.append(obs)
    
    return filtered_list

s1_maps_filtered = filter_input_maps(s1_maps_all)
s2_maps_filtered = filter_input_maps(s2_maps_all)
s3_maps_filtered = filter_input_maps(s3_maps_all)
s4_maps_filtered = filter_input_maps(s4_maps_all)

print(len(s1_maps_filtered), "maps in subfield-1")
print(len(s2_maps_filtered), "maps in subfield-2")
print(len(s3_maps_filtered), "maps in subfield-3")
print(len(s4_maps_filtered), "maps in subfield-4")

# Evenly split input lists into 2 batches
def split_input_maps(input_list):
    
    index = 0
    bundle_a = []
    bundle_b = []
    
    for obs in input_list:
        if (index % 2) == 0:
            bundle_a.append(obs)
            index += 1
        else:
            bundle_b.append(obs)
            index += 1
            
    return bundle_a, bundle_b

s1_maps_a, s1_maps_b = split_input_maps(sorted(s1_maps_filtered))
s2_maps_a, s2_maps_b = split_input_maps(sorted(s2_maps_filtered))
s3_maps_a, s3_maps_b = split_input_maps(sorted(s3_maps_filtered))
s4_maps_a, s4_maps_b = split_input_maps(sorted(s4_maps_filtered))

# Combine subfield lists
input_list_a = s1_maps_a + s2_maps_a + s3_maps_a + s4_maps_a
input_list_b = s1_maps_b + s2_maps_b + s3_maps_b + s4_maps_b
print(len(input_list_a)+len(input_list_b), "150 GHz maps total")

# Do the coadding (will take about a minute per map)
output_directory = "/spt/user/joshuasobrin/clusters/2019-08-31_map_products/"
combining_maps.coadd_maps(input_list_a, map_id="150GHz", output_file=output_directory+"partial_coadd_201906_150_a.g3")
combining_maps.coadd_maps(input_list_b, map_id="150GHz", output_file=output_directory+"partial_coadd_201906_150_b.g3")