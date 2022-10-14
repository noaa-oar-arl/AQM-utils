#!/bin/ksh -x

PDY=2019071512
UTIL_BASE=/scratch1/RDARCH/rda-arl-gpu/YouHua.Tang/UFS/AQM-utils
DATA_BASE=/scratch1/RDARCH/rda-arl-gpu/YouHua.Tang/expt_dirs/aqm_cold_aqmna13_1day

. $DATA_BASE/var_defns.sh

export NX=$ESGgrid_NX
export NY=$ESGgrid_NY
export LAYOUT_X
export LAYOUT_Y
export FCST_LEN_HRS

let NSTEP=$FCST_LEN_HRS+1

if [ $RUN_ENVIR = 'nco' ]; then
 RUN_DIR=$STMP/tmpnwprd/$RUN/$PDY
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
  $UTIL_BASE/python_utils/stack-pt-merge.py -s $PDY -n $NSTEP
fi
if [ ! -s $DATA_BASE/$PDY/PT/pt-$PDY.nc ]; then 
  echo " can not find $DATA_BASE/$PDY/PT/pt-$PDY.nc"
  exit 1
fi

cd $RUN_DIR
if [ ! -s PT/pt-0000.nc ]; then 
 let NPE=$LAYOUT_X*$LAYOUT_Y

 export TOPO=$NEXUS_FIX_DIR/$NEXUS_GRID_FN
 export PT_IN=$RUN_DIR/PT/pt-$PDY.nc

 export SLURM_ACCOUNT=naqfc
 time srun -n $NPE --time=30:00 -q debug $UTIL_BASE/exec/decomp-ptemis-mpi
fi
if [ ! -s PT/pt-0000.nc ]; then
   echo "PT emission decomposition failed"
   exit 1
else
 exit 0
fi    
