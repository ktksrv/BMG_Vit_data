#!/bin/bash

# Loop through all directories matching shear_*
for dir in Latest_shear_*; do
    if [ -d "$dir" ]; then
        echo "Submitting job in $dir"
        cd "$dir"
        sbatch LAMMPS_slurmscript.sh
        cd ..
    fi
done