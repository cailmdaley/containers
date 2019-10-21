#!/bin/bash 

#NUM_NODES=1     # number of nodes to use 
#NUM_THREADS=28  # number of cores to use
SIM_I=141          # starting index of sim number 
SIM_F=150          # ending index of sim number
CMB=2              # lens either CMB1 or CMB2 (mainly for lensing)

dir_out=/project2/chihway/yuuki/sptsz2500d_v7/       # directory of sims
lenspix=/project2/chihway/packages/lenspix/simlens   # lenspix executable
#-----------------------------------------------------------------------------------------------------
mkdir ${dir_out}
mkdir ${dir_out}/lensed_TQU${CMB}phi1
mkdir ${dir_out}/lensed_TQU${CMB}phi1/logs
mkdir ${dir_out}/configs

for i in $(eval echo "{$SIM_I..$SIM_F}")
do
echo '#!/bin/bash' > submit_job
echo "#SBATCH -t 02:00:00">> submit_job
echo "#SBATCH -o ${dir_out}/lensed_TQU${CMB}phi1/logs/outputfile_lenspix_tqu${CMB}phi1_${i}" >> submit_job
echo "#SBATCH -e ${dir_out}/lensed_TQU${CMB}phi1/logs/errorfile_lenspix_tqu${CMB}phi1_${i}" >> submit_job
echo "#SBATCH --partition=broadwl" >> submit_job
echo "#SBATCH --account=pi-chihway" >> submit_job
echo "#SBATCH --exclusive" >> submit_job
echo "#SBATCH --nodes=4" >> submit_job
echo "#SBATCH --job-name=lenspix${i}" >> submit_job
echo "export OMP_PROC_BIND=spread" >> submit_job
echo "export OMP_PLACES=threads" >> submit_job
echo "export I_MPI_PMI_LIBRARY=/software/slurm-current-$DISTARCH/lib/libpmi.so" >> submit_job
echo "cd $PWD" >> submit_job
cp lenspix/lenspix_template.ini ${dir_out}/configs/lenspix${i}.ini
echo "set_B_zero = F" >> ${dir_out}/configs/lenspix${i}.ini
echo "set_E_zero = F" >> ${dir_out}/configs/lenspix${i}.ini
echo "phi_file = ${dir_out}/input_phi1/phi_planck2018_base_plikHM_TTTEEE_lowl_lowE_lensing_seed${i}.fits" >> ${dir_out}/configs/lenspix${i}.ini
echo "TQU_file = ${dir_out}/unlensed_TQU${CMB}/TQU${CMB}_planck2018_base_plikHM_TTTEEE_lowl_lowE_lensing_seed${i}.fits" >> ${dir_out}/configs/lenspix${i}.ini
echo "out_file_root = ${dir_out}/lensed_TQU${CMB}phi1/lensedTQU${CMB}phi1_planck2018_base_plikHM_TTTEEE_lowl_lowE_lensing_seed${i}" >> ${dir_out}/configs/lenspix${i}.ini
echo "rand_seed = $i" >> ${dir_out}/configs/lenspix${i}.ini
echo source /project2/chihway/setup/setup_intel.sh >> submit_job
echo srun -n48  ${lenspix} ${dir_out}/configs/lenspix${i}.ini >> submit_job
sbatch submit_job
done

