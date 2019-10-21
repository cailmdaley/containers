import numpy, sys
from spt3g import core, gcp, std_processing

pipe = core.G3Pipeline()
pipe.Add(std_processing.ARCTimerangeReader, start_time='07-Nov-2016:08:05:00', stop_time='07-Nov-2016:08:06:00', basedir='/buffer/arc/')

registers = []
def printGCPRegisters(fr):
    if len(registers) > 0:
        return
    for key in fr.keys():
        if len(fr[key].keys()) > 0:
            for key2 in fr[key].keys():
                if len(fr[key][key2].keys()) > 0:
                    for key3 in fr[key][key2].keys():
                        registers.append(str([key])+str([key2])+str([key3]))
                else:
                    registers.append(str([key])+str([key2]))

        else:
            registers.append(str([key]))

    return

pipe.Add(printGCPRegisters)
pipe.Run(profile=False)

f = open('gcp_registers.txt', 'w')

for i in range(len(registers)):
    f.write(registers[i]+'\n')

f.close()
