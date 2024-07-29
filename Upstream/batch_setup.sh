#!/bin/bash -l
#SBATCH --job-name=p1
#SBATCH -o p1.out
#SBATCH --partition=compute

#SBATCH -n 12
#SBATCH -N 1
#SBATCH --mem=50G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bennis@pennstatehealth.psu.edu
echo Running on 'hostname'.
echo p1 is started 'date'.
echo "starting"

## TODO: Change project number here
python3 /gpfs/Cores/GSC/share/03_Analysis/H087/batch_script_mod.py H087
