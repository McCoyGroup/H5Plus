#!/bin/bash
#SBATCH --job-name=H5_o_orig
#SBATCH --nodes=1
#SBATCH --time=0-08:00:00
#SBATCH --mem=118G
#SBATCH --workdir=/gscratch/chem/b3m2a1/H5+
#SBATCH --partition=chem
#SBATCH --account=chem
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

wolframscript -file /gscratch/chem/b3m2a1/H5+/jobs/SpectrumScript.wl H5OldOriginalMore

exit 0
