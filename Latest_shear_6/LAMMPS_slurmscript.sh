#!/bin/sh
#SBATCH --job-name=PdNiP # Job name
#SBATCH -A cdsphani_nsm
#SBATCH --ntasks-per-node=48 # Number of tasks per node
#SBATCH --nodes=2
#SBATCH --time=12:00:00 # Time limit hrs:min:sec
#SBATCH --partition=small

echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"
#echo "CPU bind Info  = $SLURM_CPU_BIND"    
echo "CPU on node  = $SLURM_CPUS_ON_NODE"         

cd $SLURM_SUBMIT_DIR
. /scratch/cdskod/spack/spack/share/spack/setup-env.sh
spack load python@3.10.10 py-pip gcc intel-oneapi-mkl@2023.2.0 intel-oneapi-mpi cmake@3.26.3%gcc@12.2.0 ninja
export PYTHONPATH=${PYTHONPATH}:/scratch/cdspani/Kartikey/lammps-2Aug2023/python
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/scratch/cdspani/Kartikey/lammps-2Aug2023/build

#mpirun -n $SLURM_NTASKS python3 MC_test_6.py > /dev/null
mpirun -n $SLURM_NTASKS  python3 yield_bs.py > /dev/null
