import numpy as np

n_channels = 64
frequencies = np.linspace(1e6, 7e6, num=n_channels)
df = frequencies[1] - frequencies[0]
frequencies = frequencies + 2*(np.random.rand(n_channels)-0.5) * (df*0.25)
chan_str = ['chan_%d'%jchan for jchan in np.arange(1,n_channels+1)]
bolo_info = [(chan_str[jchan], frequencies[jchan]) for jchan in range(n_channels)]

f_bolos = open('dongle_bolos.csv', 'w')
f_bolos.write("bolometer\tchannel\n")
for jsquid in range(4):
    nchannel = 1
    for bolo in bolo_info:
        f_bolos.write('%s_%d\t0112/1/%d/%d\n'%(bolo[0], jsquid+1, jsquid+1, nchannel))
        nchannel = nchannel + 1
f_bolos.close()

f_bias = open('dongle_bias.json', 'w')
for jsquid in range(4):
    nchannel = 1
    for bolo in bolo_info:
        f_bias.write("%s_%d: {\n    frequency: %.1f\n}\n" % (bolo[0],jsquid+1,bolo[1]))
f_bias.close()

f_python = open('channels.py', 'w')
for jsquid in range(4):
    nchannel = 1
    for bolo in bolo_info:
        f_python.write("\'%s_%d\': %.1f,\n" % (bolo[0],jsquid+1,bolo[1]))
f_python.close()
