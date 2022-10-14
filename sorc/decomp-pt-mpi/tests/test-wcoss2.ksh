#!/bin/ksh -x

PDY=2019080512
CDATE=${PDY:0:8}
CYC=${PDY:8:2}
UTIL_BASE=/lfs/h2/emc/physics/noscrub/Youhua.Tang/UFS/AQM-utils/
DATA_BASE=/lfs/h2/emc/physics/noscrub/Youhua.Tang/expt_dirs/aqm_cold_aqmna13_1day

. $DATA_BASE/var_defns.sh

export NX=$ESGgrid_NX
export NY=$ESGgrid_NY
export LAYOUT_X
export LAYOUT_Y
export FCST_LEN_HRS

let NSTEP=$FCST_LEN_HRS+1

if [ $RUN_ENVIR = 'nco' ]; then
# RUN_DIR=$STMP/tmpnwprd/$RUN/$PDY
 RUN_DIR=${COMIN_BASEDIR}/${RUN}.${CDATE}/${CYC}
elif [ $RUN_ENVIR = 'community' ]; then
 RUN_DIR=$DATA_BASE/$PDY
else
 echo " unknown run environment $RUN_ENVIR "
 exit 1
fi 

cd $RUN_DIR
if [ ! -s PT/pt-$PDY.nc ]; then 
  mkdir PT
  cd PT
# need python with netCDF4
#module load hpc-miniconda3/4.6.14
#module load ufswm/1.0.0
  $UTIL_BASE/python_utils/stack-pt-merge-wcoss2.py -s $PDY -n $NSTEP
fi
if [ ! -s $RUN_DIR/PT/pt-$PDY.nc ]; then 
  echo " can not find $DATA_BASE/$PDY/PT/pt-$PDY.nc"
  exit 1
fi

cd $RUN_DIR
if [ ! -s PT/pt-0000.nc ]; then 
let NPE=$LAYOUT_X*$LAYOUT_Y
export NPE
cat>run-decomp-ptemis-mpi.bash<<!
#!/bin/bash -x
#PBS -N run_decomp_ptemis_mpi
#PBS -l place=vscatter,select=27:ncpus=64:mem=80GB
#PBS -l walltime=30:00
#PBS -A AQM-DEV 
#PBS -q debug
#PBS -V
source $UTIL_BASE/build.ver.wcoss2
. /etc/profile
module use $UTIL_BASE/modulefiles
module load wcoss2.intel

export TOPO=$NEXUS_FIX_DIR/$NEXUS_GRID_FN
export PT_IN=$RUN_DIR/PT/pt-$PDY.nc

export SLURM_ACCOUNT=naqfc
cd $RUN_DIR
mpiexec -n $NPE -cpu-bind core $UTIL_BASE/exec/decomp-ptemis-mpi
!
qsub run-decomp-ptemis-mpi.bash
fi
if [ ! -s PT/pt-0000.nc ]; then
   echo "PT emission decomposition failed"
   exit 1
else
 exit 0
fi    
