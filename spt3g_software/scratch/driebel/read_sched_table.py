from bs4 import BeautifulSoup
from dsr_tools import Cycle
import re



def read_sched_table(target = 'latest_complete_cycle'):
    logfile = '/home/sptdaq/public_html/logtable.html'

    with open(logfile) as f:
        soup = BeautifulSoup(f, features="html5lib")

    all_sched_table = soup.find_all('table')[0]
    all_cycles = all_sched_table.find_all(class_='cycle')

        # Finding a specific time stamp:
        # target = '190509 04:19:41'
    

    if target == 'latest_complete_cycle':
        target_cycle_index = -2
    else:
        for i in range(len(all_cycles)):
            if all_cycles[i].td.contents[0].strip() == target:
                target_cycle_index = i
                break
            
    one_cadence=[]
    for text in all_cycles[target_cycle_index].next_siblings: 
        if text == all_cycles[target_cycle_index + 1]: 
            break 
        if not (text == '\n'):
            one_cadence.append(text) 

    one_cadence = [all_cycles[target_cycle_index]] + one_cadence  # Add the cycle_tune itself. one_cadence is a list of tags
    cycle_start_time = one_cadence[0].find_all('td')[0].contents[0].strip()
    cycle_stop_time = one_cadence[0].find_all('td')[1].contents[0].strip()
    schedule_end = one_cadence[-1].find_all('td')[1].contents[0].strip()
    cycle_name = one_cadence[0].find_all('td')[2].contents[0].strip()
    first_scan_name = one_cadence[1].find_all('td')[2].contents[0].strip() 
    last_scan_name = one_cadence[-1].find_all('td')[2].contents[0].strip() 

    dither_string = first_scan_name.split('_')[5].split()[2]
    dither_pattern = re.compile(r"[0-9]+")
    dither=dither_pattern.search(dither_string).group()

    num_pattern = re.compile(r"[0-9]+") # Extract the cycle num from the cycle_tune call
    
    try:
        cycle_num = num_pattern.search(cycle_name).group() #Comes as a string, which is fine
    except AttributeError:
        cycle_num = 0
    

    calib_boolean = first_scan_name.split('_')[5].split()[3] # the boolean to run calibs or not.  Determines A from C, B from D
    calib_pattern = re.compile('true')

    if first_scan_name.split('_')[4] in ['2','3']:
        cycle_type = 'A'
        if calib_pattern.match(calib_boolean):
            cycle_type = 'C'
    else:
        cycle_type = 'B'
        if calib_pattern.match(calib_boolean):
            cycle_type = 'D'

    return Cycle(cycle_num, cycle_type, dither, cycle_start_time, cycle_stop_time, schedule_end) 
    

    












