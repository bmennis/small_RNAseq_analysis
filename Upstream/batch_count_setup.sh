#!/bin/bash -l
#SBATCH --job-name=p2
#SBATCH -o p2.out
#SBATCH --partition=compute

#SBATCH -n 12
#SBATCH -N 3
#SBATCH --mem=50G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bennis@pennstatehealth.psu.edu
echo Running on 'hostname'.
echo p1 is started 'date'.
echo "starting"

python3 /gpfs/Cores/GSC/share/03_Analysis/H087/batch_count_script_mod.py H087
