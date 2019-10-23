import os
import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
from glob import glob

from spt3g import core, calibration

def plot_weighthist(obs_id, datapath='/spt/user/ddutcher/field_obs/'):
    if not isinstance(obs_id, str):
        obs_id = str(obs_id)
        
    output_path = glob(datapath+obs_id+'/run_stats.pkl')
    if len(output_path)<1:
        print('No map run statistics found')
        return
    
    f = open(output_path[0],'rb')
    runstats = pkl.load(f)
    f.close()
    
    bolos_90, bolos_150, bolos_220 = [], [], []
    for bolo in runstats['TodWeightLists'].keys():
        freq = float(bolo.split('.')[-1])
        if freq < 2500:
            bolos_90.append(np.mean(runstats['TodWeightLists'][bolo].values()))
        elif freq > 3600:
            bolos_220.append(np.mean(runstats['TodWeightLists'][bolo].values()))
        else:
            bolos_150.append(np.mean(runstats['TodWeightLists'][bolo].values()))
            
    plt.figure()
    plt.hist(bolos_90, bins=50,alpha = 0.5, color='r', label='90GHz')
    plt.figure()
    plt.hist(bolos_150, bins=50,alpha = 0.5, color='g', label='150GHz')
    plt.figure()
    plt.hist(bolos_220, bins=50, alpha = 0.5, color='b', label='220GHz')
    plt.ylabel('# of bolos')
    plt.xlabel('Observation-averaged bolometer weights')
    plt.legend()
    plt.title('Distribution of PSD-based weights')
    plt.suptitle(obs_id)
    
    numbolos = len(runstats['TodWeightLists'].keys())
    numbelow_p5 = len(np.where(np.array(bolos_90+bolos_150+bolos_220)<0.5)[0])
    print('Bolometers with weights: %d'%numbolos)
    print('Percentage of those bolometers in plot: %.1f%%'%(100*float(numbelow_p5)/numbolos))
    
def print_flagstats(obs_id, datapath='/spt/user/ddutcher/field_obs/'):
    if not isinstance(obs_id, str):
        obs_id = str(obs_id)
        
    output_path = glob(datapath+obs_id+'/run_stats.pkl')
    if len(output_path)<1:
        print('No map run statistics found')
        return
    
    f = open(output_path[0],'rb')
    runstats = pkl.load(f)
    f.close()
    
    cut_stat={}
    max_scan = 0
    for k, d in runstats['FlagLists'].items():
        for scan, v in d.items():
            if scan > max_scan:
                max_scan = scan
            for reason in v:
                cut_stat[reason] = cut_stat.get(reason, 0) + 1
    print('Average number of bolometers flagged per scan:')
    for reason, total in cut_stat.items():
        print('%s\t%.1f'%(reason, total/(max_scan+1.0)))
        
def glitch_stuff(obs_id, datapath='/spt/user/ddutcher/field_obs/'):
    if not isinstance(obs_id, str):
        obs_id = str(obs_id)
        
    output_path = glob(datapath+obs_id+'/run_stats.pkl')
    if len(output_path)<1:
        print('No map run statistics found')
        return
    
    f = open(output_path[0],'rb')
    runstats = pkl.load(f)
    f.close()
    
    glitchy_bolos = {}
    
    for k, d in runstats['FlagLists'].items():
        for scan, v in d.items():
            if 'Glitchy' in v:
                if k not in glitchy_bolos.keys():
                    glitchy_bolos[k] = 0
                glitchy_bolos[k] += 1
                    
    bolos_90, bolos_150, bolos_220 = [], [], []    
    for bolo in glitchy_bolos.keys():
        freq = float(bolo.split('.')[-1])
        if freq < 2500:
            bolos_90.append(bolo)
        elif freq > 3600:
            bolos_220.append(bolo)
        else:
            bolos_150.append(bolo)
            
    print('Bolometers ever cut for glitches:')
    print('Glitchy 90s: %d'%len(bolos_90))
    print('Glitchy 150s: %d'%len(bolos_150))
    print('Glitchy 220s: %d'%len(bolos_220))
    
    # Something about distribution of glitches.
    # Number of glitches of each sigma per scan?
    # Focal plane position of each bolo that saw glitches each scan?

def plot_focalplane_glitches(obs_id, datapath='/spt/user/ddutcher/field_obs/'):
    '''
    This was more experimental, command line stuff.
    Don't use as canned function.
    '''
    if not isinstance(obs_id, str):
        obs_id = str(obs_id)
        
    output_path = glob(datapath+obs_id+'/run_stats.pkl')
    if len(output_path)<1:
        print('No map run statistics found')
#        return
    f = open(output_path[0],'rb')
    runstats = pkl.load(f)
    f.close()
    
    calframe = glob('/spt/user/production/calibration/calframe/ra0hdec-57.5/'+obs_id+'.g3')
    if not len(calframe) == 1:
        print('Calibration frame missing somehow?')
#        return
    cal = list(core.G3File(calframe[0]))[0]
    boloprops = cal['BolometerProperties']
    
    focalplane_x = []
    focalplane_y = []
    
    for bolo in boloprops.keys():
        focalplane_x.append(boloprops[bolo].x_offset)
        focalplane_y.append(-1 * boloprops[bolo].y_offset)
        
    glitchy_bolos = {}
    
    for k, d in runstats['FlagLists'].items():
        for scan, v in d.items():
            if v==['Glitchy']:
                if k not in glitchy_bolos.keys():
                    glitchy_bolos[k] = 0
                glitchy_bolos[k] += 1
                
 #   morethan5 = [bolo for bolo in glitchy_bolos.keys() if glitchy_bolos[bolo] > 5]
##
    import time
    import os
    from glob import glob
    import pickle as pkl
    from spt3g import core, calibration
    
    obs_id = '51832227'
    plt.ion()
    calframe = glob('/spt/user/production/calibration/calframe/ra0hdec-59.75/'+obs_id+'.g3')
    if not len(calframe) == 1:
        print('Calibration frame missing somehow?')
#        return
    cal = list(core.G3File(calframe[0]))[0]
    boloprops = cal['BolometerProperties']
    focalplane_x = []
    focalplane_y = []
    
    for bolo in boloprops.keys():
        focalplane_x.append(boloprops[bolo].x_offset)
        focalplane_y.append(-1 * boloprops[bolo].y_offset)
        
#     datadir = '/home/ddutcher/data/flags/ra0hdec-52.25/'
    run = 'oneOver20_glitch_wafCM_poly19_mhpf150_lr'
    datadir = '/spt/user/ddutcher/ra0hdec-59.75/'+run
    runstats = pkl.load(open(
        os.path.join(datadir, obs_id, run+'_'+obs_id+'_run_stats.pkl'),'rb'))
    for scan in range(0,100):      
        glitchy_x = []
        glitchy_y = []
        for bolo in runstats['FlagLists'].keys():
            if scan in runstats['FlagLists'][bolo] and \
            list(runstats['FlagLists'][bolo][scan])==['Glitchy']:
#                 if runstats['GlitchLists'][bolo][scan][0]==1 and runstats['GlitchLists'][bolo][scan][1]<5:
                glitchy_x.append(boloprops[bolo].x_offset)
                glitchy_y.append(-1*boloprops[bolo].y_offset)
        plt.figure(1)
        plt.clf()
        plt.plot(focalplane_x, focalplane_y, 'o', mfc = 'None', mec='gray', markersize=4)
        plt.plot(glitchy_x, glitchy_y, 'o', markersize=5,color='red')
        plt.title(str(scan))
        plt.draw()
        plt.pause(.5)        
    
    
    glitchy_lte5 = zeros(166)
    glitchy_mt5 = zeros(166)
    for scan in range(0,166):      
        for bolo, count in glitchy_bolos.items():
            if scan in runstats['FlagLists'][bolo] and list(
                runstats['FlagLists'][bolo][scan])==['Glitchy']:
                if count>5:
                    glitchy_mt5[scan]+=1
                else:
                    glitchy_lte5[scan]+=1
    plot(glitchy_mt5,'bo')
    plot(glitchy_lte5,'go')
    
# for scan in range(166):
#     sigma5=[]
#     sigma10=[]
#     sigma20=[]
#     for bolo, d in runstats['GlitchLists'].items():
#         if bolo in fiveormore and scan in d:
#             sigma5.append(d[scan][0])
#             sigma10.append(d[scan][1])
#             sigma20.append(d[scan][2])
#     plt.figure(1)
#     plt.clf()
#     plt.hist(sigma5,bins=29,range=(1,30))
#     plt.hist(sigma10,bins=29,range=(1,30))
#     plt.hist(sigma20,bins=29,range=(1,30))
#     plt.title(str(scan))
#     plt.pause(0.5)

def print_flags_wiki(field,
                     datapath='/home/ddutcher/data/flags/'):
    obs_list = sorted(glob(os.path.join(datapath, field,'*run_stats.pkl')))[::2]
    print('Found %s observations'%len(obs_list))
    total_stats={}
    for path in obs_list:
        obs_id = path.split('/')[-1].split('_')[-3]
        try:
            f = open(path,'rb')
            obs = pkl.load(f)
            f.close()
        except EOFError:
            continue
        print(obs_id)
        
        max_scan = 0
        obs_stats = {}
        for bolo, scanflags in obs['FlagLists'].items():
            for scan, flags in scanflags.items():
                if scan > max_scan: 
                    max_scan = scan
                for flag in flags:
                    if flag not in obs_stats:
                        obs_stats[flag] = 0
                    obs_stats[flag] += 1
        for flag, val in obs_stats.items():
            obs_stats[flag] = val/(max_scan+1.0)
        flags = sorted(obs_stats.keys())
        vals = [format(obs_stats[flag],'.1f') for flag in flags]
        total_stats[obs_id] = flags, vals
        
    all_flags = []
    for obs, fv in total_stats.items():
        for flag in fv[0]:
            if not flag in all_flags:
                all_flags.append(flag)
    all_flags=sorted(all_flags)
    print('! Obs')
    for flag in all_flags:
        print("! %s"%flag)
    print('|-')
    for obs, fv in sorted(total_stats.items()):
        to_print = "|| '''%s'''"%obs
        for flag in all_flags:
            if flag in fv[0]:
                to_print += ' || '+fv[1][fv[0].index(flag)]
            else:
                to_print += ' || 0'
        print(to_print)
        print('|-')

        
if __name__=='__main__':
    ### This is all just copy-paste junk
    results = dict()
    fields = ['44.75','52.25','59.75','67.25']
    for field in fields:
        results['ra0hdec-'+field]=dict()
        now=results['ra0hdec-'+field]
        runpkls = glob('/spt/user/ddutcher/ra0hdec-'+field+'/poly9_mhpf50/4*/*run_stats.pkl')
        for _pkl in sorted(runpkls):
            obs = _pkl.split('/')[-2]
            runstats = pkl.load(open(_pkl,'rb'))
            now[obs] = dict()
            for bolo in runstats['FlagLists'].keys():
                for scan, flags in runstats['FlagLists'][bolo].items():
                    if scan not in now[obs].keys():
                        now[obs][scan] = 0
                    if 'Glitchy' in flags:
                        now[obs][scan] += 1
        plt.figure()
        plt.hist(by_scan,range=(0,121),bins=121,align='left')
        plt.xlim(-1,121)
        plt.title(obs)