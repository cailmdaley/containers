import argparse, math, os

def get_eddy_sub_script(commands,
                        mpibatch_cl_path = '/global/homes/n/nlharr/spt_code/spt3g_software/clusterutils/mpibatch_cl.x',
                        n_parallel_submits = 1,
                        walltime = '00:10:00',
                        queue = 'regular',
                        job_name  = 'TestSubmitJob',
                        
                        ):
    nodesize = 24
    use_hyperthreading = False
    openmp_num_threads = 6
    n_jobs_per_node = 4

    mppwidth = nodesize * n_parallel_submits 
    
    
    #-S mpi tasks per numa node
    #-N number per node
    #-n number of cores per a job
    header = '''#PBS -q %s
#PBS -l mppwidth=%d
#PBS -l walltime=%s
#PBS -N %s
#PBS -j oe 

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=6

export PYTHONPATH=$PYTHONPATH:%s

''' % ( queue, mppwidth, walltime, job_name, os.environ['SPT3G_SOFTWARE_BUILD_PATH']  )
    job_submit_command = '''aprun -n 4 -N 4 -d 6 -S 2 %s'''%mpibatch_cl_path

    out_str = header    
    for i, c in enumerate(commands):
        if i % n_jobs_per_node == 0:
            out_str += job_submit_command
        out_str += ' %s'%c
        if (i+1) % n_jobs_per_node == 0:
            out_str += '&\n'
        if (i+1) % (n_jobs_per_node * n_parallel_submits) == 0:
            out_str += 'wait\n'
    if out_str[-1] != '\n':
        out_str += '&\n'
    if out_str[-5:] != 'wait\n':
        out_str += 'wait\n'

    return out_str

if __name__ == '__main__':
    in_lst = ['/scratch2/scratchdirs/nlharr/frb_input/ra23h30dec-55_idf_20120828_072741_150ghz.h5',
              '/scratch2/scratchdirs/nlharr/frb_input/ra23h30dec-55_idf_20120608_184036_150ghz.h5',
              '/scratch2/scratchdirs/nlharr/frb_input/ra23h30dec-55_idf_20120609_001450_150ghz.h5',
              '/scratch2/scratchdirs/nlharr/frb_input/ra23h30dec-55_idf_20120609_043325_150ghz.h5'
              ]
    out_lst = ['/scratch2/scratchdirs/nlharr/frb_test_output/out_%d.g3'%i for i in range(4)]

    cmd_lst = []
    for i in range(len(in_lst)):
        cmd_lst.append('"python /global/homes/n/nlharr/spt_code/spt3g_software/frbutils/scripts/sptpol_idf_to_possible_events.py --input_file %s --output_file %s"'%(in_lst[i], out_lst[i]))
    print(get_eddy_sub_script(cmd_lst, walltime = '40:00'))

