from bs4 import BeautifulSoup
from dsr_tools import Cycle
import re


def last_eight_cycles():
    logfile = '/home/sptdaq/public_html/logtable.html'

    with open(logfile) as f:
        soup = BeautifulSoup(f, features="html5lib")

    all_sched_table = soup.find_all('table')[0]
    all_cycles = all_sched_table.find_all(class_='cycle')


    cycle_list = []
    for i in range(-9,-1):
        one_cadence=[]
        for text in all_cycles[i].next_siblings: 
            if text == all_cycles[i+1]: 
                break 
            if not (text == '\n'):
                one_cadence.append(text)
        one_cadence = [all_cycles[i]] + one_cadence
        first_scan_name = one_cadence[1].find_all('td')[2].contents[0].strip()
        cycle_name = one_cadence[0].find_all('td')[2].contents[0].strip()
        cycle_start_time = one_cadence[0].find_all('td')[0].contents[0].strip()
        cycle_stop_time = one_cadence[0].find_all('td')[1].contents[0].strip()
        schedule_end = one_cadence[-1].find_all('td')[1].contents[0].strip()
        dither_string = first_scan_name.split('_')[5].split()[2]
        dither_pattern = re.compile(r"[0-9]+")
        dither=dither_pattern.search(dither_string).group()
        num_pattern = re.compile(r"[0-9]+")
        try:
            num_pattern.search(cycle_name).group() #Comes as a string, which is fine
        except AttributeError:
            cycle_num = 0
        else:
            cycle_num = num_pattern.search(cycle_name).group()

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
        cycle = Cycle(cycle_num, 
                      cycle_type, 
                      dither, 
                      cycle_start_time, 
                      cycle_stop_time, 
                      schedule_end)
        cycle_list.append(cycle)

    cycle_types = [(i.cycle_num, i.cadence) for i in cycle_list]
    return cycle_types

if __name__ == '__main__':
    a = last_eight_cycles()
    print('The last eight completed cycles were:')
    for i in a:
        print('{}: {}'.format(i[0],i[1]))






