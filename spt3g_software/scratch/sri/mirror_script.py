import sys,os
fname = sys.argv[1:-1]
machine = sys.argv[-1]

cmd = 'rsync -trvz --delete %s %s:/home/sri/softwares/spt3g_software/scratch/sri/' %(' '.join(fname), machine)
os.system(cmd)
