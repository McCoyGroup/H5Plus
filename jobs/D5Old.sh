#!/bin/bash
#SBATCH --job-name=D5_old_spectrum
#SBATCH --nodes=1
#SBATCH --time=0-08:00:00
#SBATCH --mem=118G
#SBATCH --workdir=/gscratch/chem/b3m2a1/H5+
#SBATCH --partition=ckpt
#SBATCH --account=stf-ckpt
#SBATCH --mail-type=ALL

# load environment
module load mathematica_11.2

# debugging information
echo "**** Job Debugging Information ****"
echo "This job will run on $SLURM_JOB_NODELIST"
echo ""
echo "ENVIRONMENT VARIABLES"
set
echo "**********************************************"

# run job

wolframscript -file /gscratch/chem/b3m2a1/H5+/jobs/SpectrumScript.wl D5Old

exit 0
