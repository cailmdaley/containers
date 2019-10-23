import datetime as dt
import pytz
from bs4 import BeautifulSoup
from last_eight_cycles import *


'''
Simple script to predict when the current Winter 2019 observation schedule is
due to complete.

It reads the latest started schedule from the schedule log table and then prints the start time (to confirm it found the right one) and the expected completion time in both UTC and Station time.

'''

def obs_predictor(cycle_start_time = 37, schedule='Q'):
    
    utc_tz = pytz.timezone('UTC')
    nz_tz = pytz.timezone('Pacific/Auckland')
    in_fmt = '%y%m%d %H:%M:%S'
    out_fmt = '%H:%M %d %b'

    if isinstance(cycle_start_time, int):
        logfile = '/home/sptdaq/public_html/logtable.html'
        with open(logfile) as f:
            soup = BeautifulSoup(f, features="html5lib")

        all_sched_table = soup.find_all('table')[0]
        all_cycles = all_sched_table.find_all(class_='cycle')
        target_cycle_index = -1
        one_cadence=[]
        for text in all_cycles[target_cycle_index].next_siblings: 
            if not (text == '\n'):
                one_cadence.append(text) 
        
        one_cadence = [all_cycles[target_cycle_index]] + one_cadence
        cycle_start_time = one_cadence[0].find_all('td')[0].contents[0].strip()
        cycle_name = one_cadence[0].find_all('td')[2].contents[0].strip()
        if len(one_cadence) > 1:
            first_scan_name = one_cadence[1].find_all('td')[2].contents[0].strip()
            calib_boolean = first_scan_name.split('_')[5].split()[3]
            calib_pattern = re.compile('true')
            if first_scan_name.split('_')[4] in ['2','3']:
                schedule = 'A'
                if calib_pattern.match(calib_boolean):
                    schedule = 'C'
            else:
                schedule = 'B'
                if calib_pattern.match(calib_boolean):
                    schedule = 'D'
        else:
            schedule = last_eight_cycles()[0] # we are in a cycle_tune.  in normal ops,
            # the current running schedule was run 8 cycles ago.

    else:
        try:
            dt.datetime.strptime(cycle_start_time, in_fmt)
        except:
            raise ValueError('Time must be formatted like in the schedule logtable')

    if schedule.upper() in ['A', 'B']:
        duration = dt.timedelta(hours=20, minutes=0)
    if schedule.upper() == 'C':
        duration = dt.timedelta(hours=19, minutes=0)
    if schedule.upper() == 'D':
        duration = dt.timedelta(hours=19, minutes=10)
    
    cycle_start = dt.datetime.strptime(cycle_start_time, in_fmt)
    cycle_start = utc_tz.localize(cycle_start)
    finish_time = cycle_start + duration

    '''print()
    print('Schedule Started: {} (NZ)'.format(dt.datetime.strftime(cycle_start.astimezone(nz_tz),out_fmt)))
    print('Cadence: {}'.format(schedule.upper()))
    print('Expected Duration: {}'.format(str(duration)))
    print('Expected Completion Time: {} (NZ)'.format(dt.datetime.strftime(finish_time.astimezone(nz_tz),out_fmt)))
    print('{} (UTC)'.format(dt.datetime.strftime(finish_time.astimezone(utc_tz),out_fmt)))
    print('{} (NZ)'.format(dt.datetime.strftime(finish_time.astimezone(nz_tz),out_fmt)))
    print()'''

    print('{:14}{:21}{:20}'.format(schedule.upper(),
                                   dt.datetime.strftime(cycle_start.astimezone(nz_tz),out_fmt),
                                   dt.datetime.strftime(finish_time.astimezone(nz_tz),out_fmt)))

    return str(dt.datetime.strftime(finish_time.astimezone(utc_tz), in_fmt))
    

if __name__ == '__main__':
    
    import argparse
    from last_eight_cycles import *
    parser = argparse.ArgumentParser()
    parser.add_argument('repeats', nargs = '?', default = 1)
    args = parser.parse_args()
    in_fmt = '%y%m%d %H:%M:%S'

    try:
        repeats = int(args.repeats)
    except:
        print('Number of repeats must be an integer')
    
    print('{:14}{:21}{:17}'.format('Cadence','Start Time (NZ)','Expected Finish (NZ)'))
    if repeats > 1:
        next_cycles = last_eight_cycles()  
        # Under normal ops, the current cadence type is always element 0 of this array
        if repeats > 8:
            next_cycles = next_cycles * ((repeats // 8) + 1)
            next_cycles = next_cycles[:repeats] 
            # next_cycles is now a list of the cadences of the right length
        
        current_cycle_end = obs_predictor() 
        repeats = repeats - 1
        next_cycles = next_cycles[1:]
        for i in range(repeats):
            current_cycle_end = obs_predictor(current_cycle_end, next_cycles[i])
    else:
        junk = obs_predictor()


