#!/bin/ksh -x

PDY=2019071512
UTIL_BASE=/scratch1/RDARCH/rda-arl-gpu/YouHua.Tang/UFS/AQM-utils
DATA_BASE=/scratch1/RDARCH/rda-arl-gpu/YouHua.Tang/expt_dirs/uwm2_cold_1day

. $DATA_BASE/var_defns.sh
. $DATA_BASE/config.sh

export NX
export NY
export LAYOUT_X
export LAYOUT_Y
export FCST_LEN_HRS

let NSTEP=$FCST_LEN_HRS+1

cd $DATA_BASE/$PDY
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

cd $DATA_BASE/$PDY
if [ ! -s PT/pt-0000.nc ]; then 
 let NPE=$LAYOUT_X*$LAYOUT_Y

 export TOPO=$DATA_BASE/orog/C775_oro_data.tile7.halo0.nc
 export PT_IN=$DATA_BASE/$PDY/PT/pt-$PDY.nc

 export SLURM_ACCOUNT=naqfc
 time srun -n $NPE --time=30:00 -q debug $UTIL_BASE/exec/decomp-ptemis-mpi
fi
if [ ! -s PT/pt-0000.nc ]; then
   echo "PT emission decomposition failed"
   exit 1
else
 exit 0
fi    
