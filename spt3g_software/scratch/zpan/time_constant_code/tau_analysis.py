# This code fits for the time constants from 
# calibrator analysis outputs
import pickle as pkl
import glob
import numpy as np
import pylab as pl
from spt3g import core,gcp,dfmux
import pydfmux
import pickle
from scipy.optimize import curve_fit        
import time, os

# Collect the calibrator response from the analyzed files
def collect_response(sweepid, calibratorfiles, Waflist, sn_cut_4hz=10, sn_cut=0):
    '''
    Calibrator files: in the format of {freq: obsid}
    sn_cut_4hz: calibrator s/n cut at 4hz
    sn_cut: calibrator s/n cut for other frequencies
    '''
    el= sweepid.split('_')[-1]+'deg'
    sweep_number= sweepid.split('_')[0].split('sweep')[-1]
    label=calibratorfiles[4]                            
    power_sorted={}                                     
    s_to_n_sorted={}                                    
    powerchange={}                                      
    s_to_n={}                                           
    powerchange_reduced={}                              
    s_to_n_reduced={}                                   
    #the total number of bolometers with non zero response                   
    total_num=0                         
    #the number of bolos with response that pass the s/n cut
    num_after_sn_cut=0
    #the number of bolos which pass the chi square cut
    num_after_chisq_cut=0            
   
    # for each frequency, go through the file and find the data  
    for freq in calibratorfiles.keys():                 
        data=core.G3File(calibratorfolder+str(calibratorfiles[freq])+'.g3').next()   
        powerchange[freq]={}                             
        s_to_n[freq]={}                                  
        powerchange_reduced[freq]={}                     
        s_to_n_reduced[freq]={}                          
        for key in data['CalibratorResponsePerBolo'].keys():    
            # make a cut on 4hz data
            if freq==4 and data['CalibratorResponsePerBoloSN'][key]>0\
               and data['CalibratorResponsePerBolo'][key]>0:
                total_num=total_num+1
                       
            #make the cut                                          
            if freq==4 and data['CalibratorResponsePerBoloSN'][key]>10\
               and data['CalibratorResponsePerBolo'][key]>0:
                try:    
                    num_after_sn_cut=num_after_sn_cut+1                
                    powerchange[freq][key]=data['CalibratorResponsePerBolo'][key]            
                    s_to_n[freq][key]=data['CalibratorResponsePerBoloSN'][key]  
                except KeyError:
                    continue             
            # maybe a different cut at some other frequencies   
                                                           
            if freq!=4: 
                try:                                                 
                    powerchange[freq][key]=data['CalibratorResponsePerBolo'][key]                
                    s_to_n[freq][key]=data['CalibratorResponsePerBoloSN'][key]
                except KeyError:
                    continue                   
    print('total number of bolos with sn>0 is', total_num)
    print ('after the sn cut is', num_after_sn_cut)                                              
    # create a dictionary to store the response vs. frequency data, by bolometer name
    # result[boloname][2] is the error bar                                     
    result={}                                     
    for boloname in powerchange[sorted(calibratorfiles.keys())[0]]:                                                       
       result[boloname]=[[],[], []]                   
       for freq in sorted(calibratorfiles.keys()):                                                         
           if boloname in powerchange[freq].keys():                                                                       
               result[boloname][0].append(freq)  
               result[boloname][1].append(powerchange[freq][boloname])  
               result[boloname][2].append(powerchange[freq][boloname]/s_to_n[freq][boloname])                            
    #normalize the response to the response at the lowest frequency                                                  
    for boloname in result.keys(): 
        result[boloname][0]=np.array(result[boloname][0])                                                                    
        result[boloname][1]=np.array(result[boloname][1])  
        result[boloname][2]=np.array(result[boloname][2]) 
        result[boloname][2]=np.array(result[boloname][2])/result[boloname][1][0]
        result[boloname][1]=np.array(result[boloname][1])/result[boloname][1][0]                                             
    return result, sweep_number

#fit the result to tau                 
def fit_tau(result, Waflist, chisqcut=50):
    #create a dictionary to store time constant by wafer                                                                  
    time_constant={}               
    for waf in Waflist:                           
        time_constant[waf]={'90':[],'150':[], '220':[]} 
    # fit for the time constant             
    def model(f,t,A):                                   
        return A/np.sqrt(1+(2*np.pi*f*t)**2)                
    # store the fit result                              
    fit_results={}                                      
    #do the fit      
    chisq=[]  
    chisq_dict={}      
    for boloname in result.keys():                      
        if len(result[boloname][0])>3:               
            try:                             
                # do the fit and calculate chisquare          
                fit_result=curve_fit(model,result[boloname][0],
                                     result[boloname][1], sigma= result[boloname][2])  
                chisqfit= np.sum(((model(result[boloname][0],abs(fit_result[0][0]),\
                          abs(fit_result[0][1]))- result[boloname][1])/result[boloname][2])**2)
                chisq_dict[boloname]=chisqfit
                chisq.append(chisqfit)
                if abs(fit_result[0][0])<1 and abs(fit_result[0][0])>1e-4 and chisqfit<200: 
                    fit_results[boloname]=abs(fit_result[0])
            except:         
                continue      
    for boloname in fit_results.keys():      
        fit_result= fit_results[boloname]
        for waf in Waflist: 
            try:
                pname= physical_names[boloname]
                if pname.split('.')[0]==waf and pname.split('.')[-2]=='90':                                             
                    time_constant[waf]['90'].append(abs(fit_result[0]))  
                if pname.split('.')[0]==waf and pname.split('.')[-2]=='150':                                             
                    time_constant[waf]['150'].append(abs(fit_result[0]))  
                if pname.split('.')[0]==waf and pname.split('.')[-2]=='220':                                             
                    time_constant[waf]['220'].append(abs(fit_result[0]))  
            except:
                continue
    #print number after the chisq cut          
    num_after_chisq_cut=len(fit_results)
    print('after the chisq cut is', num_after_chisq_cut)   
    return fit_results, time_constant

# Plot the fit results and histogram the time constants
def plot_tau(response,fit_results, time_constant, Waflist, save_folder, el, label):
    #make the scatter plots     
    figurenumber=0                                      
    frequencies=np.linspace(4,70, 100)
    for waf in Waflist:  
        pl.figure(figurenumber, figsize=(10,8))   
        pl.clf()                                          
        figurenumber=figurenumber+1                       
        pl.subplot(311)                                   
        for boloname in fit_results.keys():              
            try:
                pname= physical_names[boloname]
                if pname.split('.')[0]==waf and pname.split('.')[-2]=='90':                                                      
                    pl.errorbar(result[boloname][0],result[boloname][1],
                                yerr=result[boloname][2], fmt='o', color= 'b')
                    pl.plot(frequencies, 
                            model(frequencies, fit_results[boloname][0], fit_results[boloname][1]),
                            'b', linewidth=1)                         
            except:
                continue
        pl.plot([],'bo', label='90GHz')                      
        pl.legend(loc='best')                                 
        pl.title(waf+' calibrator response vs. freq, el='+el) 
        pl.subplot(312)                                       
        for boloname in fit_results.keys():   
            try:
                pname= physical_names[boloname]
                if pname.split('.')[0]==waf and pname.split('.')[-2]=='150':                                               
                    pl.errorbar(result[boloname][0],result[boloname][1],
                                yerr=result[boloname][2], fmt='o', color= 'g')
                    pl.plot(frequencies, 
                            model(frequencies, fit_results[boloname][0], fit_results[boloname][1]),
                            'g', linewidth=1)               
            except:
                continue
        pl.plot([],'go', label='150GHz')           
        pl.legend(loc='best')                       
        pl.ylabel('Response normalized to 4Hz')     
        pl.subplot(313)                             
        for boloname in fit_results.keys():     
            try:
                pname= physical_names[boloname]    
                if pname.split('.')[0]==waf and pname.split('.')[-2]=='220':                                 
                    pl.errorbar(result[boloname][0],result[boloname][1],
                                yerr=result[boloname][2], fmt='o', color= 'r')
                    pl.plot(frequencies, 
                            model(frequencies, fit_results[boloname][0], fit_results[boloname][1]),
                            'r', linewidth=1)               
            except:
                continue
        pl.plot([],'ro', label='220GHz')           
        pl.legend(loc='best')        
        pl.xlabel('Calibrator Frequency (Hz)')       
        pl.legend(loc='best')        
        pl.savefig(save_folder+waf+'calibrator_response'+el+'_'+str(label)+'.png')                                                           
    #make the histograms of time constant                        
    figurenumber=0                                
    for waf in sorted(Waflist):                           
        pl.figure(figurenumber)                     
        pl.clf()                                    
        figurenumber=figurenumber+1 
        try:                                                                    
            pl.hist([np.array(time_constant[waf]['90'])*1000, np.array(time_constant[waf]['150'])*1000,\
                     np.array(time_constant[waf]['220'])*1000],
                    bins=np.linspace(0,17,28), alpha=0.8, color=['b', 'g', 'r'],
                    label=['90GHz', '150GHz', '220GHz'], stacked=True)
        except:
            continue                     
        data_90ghz=np.array(time_constant[waf]['90'])\
                   [np.where(np.array(time_constant[waf]['90'])<0.02)]*1000
        data_150ghz=np.array(time_constant[waf]['150'])\
                    [np.where(np.array(time_constant[waf]['150'])<0.02)]*1000
        data_220ghz=np.array(time_constant[waf]['220'])\
                    [np.where(np.array(time_constant[waf]['220'])<0.02)]*1000
        print('wafer '+waf+' 90GHz average time constant and std are',  
              np.mean(data_90ghz), np.std(data_90ghz)/np.mean(data_90ghz))
        print('wafer '+waf+' 150GHz average time constant and std are',  
              np.mean(data_150ghz), np.std(data_150ghz)/np.mean(data_150ghz))
        print('wafer '+waf+' 220GHz average time constant and std are',  
              np.mean(data_220ghz), np.std(data_220ghz)/np.mean(data_220ghz))
        pl.legend()
        pl.ylabel('Count')  
        pl.xlabel('Time constant (ms)')                        
        pl.title(waf+ ' time constant fit result, el='+el)                                      
        pl.savefig(save_folder+waf+'calibrator_time_constant'+el+'_'+str(label)+'.png')                                                            
    pl.close('all')



# The master module. A wrapper for all above functions.
def analyze_tau(all_cal_sweeps, cal_dir,
                out_dir,  year='2019', obsid_range= [0, 1e10],
                sn_cut_4hz= 10, chisqcut= 100):
    global calibratorfolder
    calibratorfolder  = cal_dir
    global outputfolder
    outputfolder = out_dir
    global physical_names

    #Waf list for 2018
    if year == '2018':
        Waflist=['w174', 'w176', 'w177', 'w188', 'w201', 'w180', 'w203', 'w172', 'w181']
        physical_names = pickle.load(open('/home/panz/useful_data/logical_to_physical_name_post_event.pkl','rb'))

    #Waf list for 2019
    if year == '2019':
        Waflist=['w172', 'w174', 'w176', 'w177', 'w180', 'w181','w188', 'w203', 'w204', 'w206']
        physical_names = pickle.load(open('/home/panz/useful_data/logical_to_physical_name_19.pkl','rb'))

    for sweepid in all_cal_sweeps.keys():  
        if int(all_cal_sweeps[sweepid][4])>obsid_range[0]\
           and int(all_cal_sweeps[sweepid][4])<obsid_range[1]:
            save_folder= outputfolder+'/'+year+'/'
            if not os.path.exists(save_folder):
                os.makedirs(save_folder)
            if not os.path.exists(save_folder+'pngs/'):
                os.makedirs(save_folder+'pngs/')
            calibratorfiles = all_cal_sweeps[sweepid]
            label = all_cal_sweeps[sweepid][4]
            el= sweepid.split('_')[-1]+'deg'
            out_fname = save_folder+'time_const'+str(label)+'_el_'+el+ '.pkl'            
            if not os.path.isfile(out_fname): 
                result, sweep_number= collect_response(sweepid, calibratorfiles,
                                                       Waflist, sn_cut_4hz=sn_cut_4hz,
                                                       sn_cut=0)
                fit_results, time_constant = fit_tau(result, Waflist, chisqcut=chisqcut)
                plot_tau(result,fit_results, time_constant, 
                         Waflist, save_folder+'pngs/', el=el, label=label)
                pickle.dump(fit_results, open(out_fname, 'wb'))
