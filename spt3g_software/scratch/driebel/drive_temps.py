import subprocess
import argparse
import re

def drive_temps():

    drives = {'P01':'sdm',
              'P02':'sdk',
              'P03':'sdl',
              'P04':'sdn',
              'P05':'sdi',
              'P06':'sdj',
              'P07':'sda',
              'P08':'sdg',
              'P09':'sdh',
              'P10':'sdb',
              'P11':'sde',
              'P12':'sdf',
              'P13':'sdc',
              'P14':'sdd',
              'P15':'sdao',
              'P16':'sdam',
              'P17':'sdan',
              'P18':'sdap',
              'P19':'sdak',
              'P20':'sdal',
              'P21':'sdac',
              'P22':'sdai',
              'P23':'sdaj',
              'P24':'sdad',
              'P25':'sdag',
              'P26':'sdah',
              'P27':'sdae',
              'P28':'sdaf',
              'P29':'sdaa',
              'P30':'sdy',
              'P31':'sdz',
              'P32':'sdab',
              'P33':'sdw',
              'P34':'sdx',
              'P35':'sdo',
              'P36':'sdu',
              'P37':'sdv',
              'P38':'sdp',
              'P39':'sds',
              'P40':'sdt',
              'P41':'sdq',
              'P42':'sdr'}

    commands = []
    for drive in drives.values():
        commands.append('smartctl -A /dev/{} | grep Current \n'.format(drive))

    ssh = subprocess.Popen(['ssh', 'root@storage01'],
                           stdout = subprocess.PIPE,
                           stderr = subprocess.PIPE,
                           stdin = subprocess.PIPE,
                           universal_newlines = True,
                           bufsize = 0)

    for command in commands:
        ssh.stdin.write(command)

    ssh.stdin.write('logout\n')
    ssh.stdin.close()

    num_pattern = re.compile(r"[0-9]+")
    temp_list = [int(num_pattern.search(line).group()) for line in ssh.stdout]
    temp_dict = dict(zip(drives.keys(), temp_list))

    return temp_dict

if __name__ == '__main__':
    temp_dict = drive_temps()
    output_string = '{:4G}'*14
    print(output_string.format(temp_dict['P29'],
                               temp_dict['P30'],
                               temp_dict['P31'],
                               temp_dict['P32'],
                               temp_dict['P33'],
                               temp_dict['P34'],
                               temp_dict['P35'],
                               temp_dict['P36'],
                               temp_dict['P37'],
                               temp_dict['P38'],
                               temp_dict['P39'],
                               temp_dict['P40'],
                               temp_dict['P41'],
                               temp_dict['P42']))
    
    print(output_string.format(temp_dict['P15'],
                               temp_dict['P16'],
                               temp_dict['P17'],
                               temp_dict['P18'],
                               temp_dict['P19'],
                               temp_dict['P20'],
                               temp_dict['P21'],
                               temp_dict['P22'],
                               temp_dict['P23'],
                               temp_dict['P24'],
                               temp_dict['P25'],
                               temp_dict['P26'],
                               temp_dict['P27'],
                               temp_dict['P28']))

    print(output_string.format(temp_dict['P01'],
                               temp_dict['P02'],
                               temp_dict['P03'],
                               temp_dict['P04'],
                               temp_dict['P05'],
                               temp_dict['P06'],
                               temp_dict['P07'],
                               temp_dict['P08'],
                               temp_dict['P09'],
                               temp_dict['P10'],
                               temp_dict['P11'],
                               temp_dict['P12'],
                               temp_dict['P13'],
                               temp_dict['P14']))
