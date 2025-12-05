#!/bin/sh
#SBATCH --job-name=PdNiP # Job name
#SBATCH -A cdsphani_nsm
#SBATCH --ntasks-per-node=48 # Number of tasks per node
#SBATCH --nodes=4
#SBATCH --time=00:05:00 # Time limit hrs:min:sec
#SBATCH --partition=debug

echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"
#echo "CPU bind Info  = $SLURM_CPU_BIND"    
echo "CPU on node  = $SLURM_CPUS_ON_NODE"         

cd $SLURM_SUBMIT_DIR
# module load openmpigcc9.2
# module load gcc9.2
# module load cuda11gcc9.2
# module load mkl/2021.2.0
#module load oneapi
#module load intelmkl20
#srun -n 2 --mpi=pmi2   /shared/LAMMPS/lammps/src/lmp_kokkos_mpi_only -k on -sf kk -in melt_quench_update > Si24C5H36_2fold_annealing.out
. /scratch/cdskod/spack/spack/share/spack/setup-env.sh
spack load python@3.10.10 py-pip gcc intel-oneapi-mkl@2023.2.0 intel-oneapi-mpi cmake@3.26.3%gcc@12.2.0 ninja
export PYTHONPATH=${PYTHONPATH}:/scratch/cdspani/Kartikey/lammps-2Aug2023/python
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/scratch/cdspani/Kartikey/lammps-2Aug2023/build
# export LD_LIBRARY_PATH="$HOME/.local/lib64:/apps/anaconda/lib:$LD_LIBRARY_PATH"
# LD_LIBRARY_PATH="/shared/LAMMPS/lammps/build_shared:$LD_LIBRARY_PATH"
# export PYTHONPATH="/shared/LAMMPS/lammps/python:$PYTHONPATH"
mpirun -n $SLURM_NTASKS python3 Normal_Stress_Parallel_0.py
# mpirun -n $SLURM_NTASKS python3 MC_good_cluster.py > /dev/null
# mpirun -n $SLURM_NTASKS python3 bond_length_0.py > /dev/null
# mpirun -n $SLURM_NTASKS python3 yield_bs.py > /dev/null






