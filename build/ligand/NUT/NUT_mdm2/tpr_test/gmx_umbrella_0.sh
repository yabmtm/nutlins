#PBS -N mdm2_nut.pf
#PBS -o mdm2_nut.out
#PBS -q normal
#PBS -l nodes=1:ppn=12
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS

cd $PBS_O_WORKDIR

module use -a /home/tuf10875/pkg/modulefiles/
module load gromacs/5.1.2-mine openmpi

mpirun -np 12 mdrun_mpi -deffnm buildit3

