# Interpolate GRIB2 GEFS-Aerosols fields to FV3 cold start files
In this directory, `ics_from_grib2.py` can be called from the command line when provided arguments to interpolate GRIB2 GEFS-Aerosols fields to FV3 cold start files for GCAFS

## Instructions

### Environment
First, please ensure that you have an appropriate python environment loaded with all of the required dependencies. See `requirements.txt` in this directory for details. You may have supported environments already set up on HPC platforms, it may be worth asking the developer of this code if one already exists.

The key components needed are:
- numpy
- grib2io (developed by NOAA-MDL; it has several dependencies)
- scipy
- xarray

### Example Usage
`python ics_from_grib2.py --grib-file /path/to/gefs.chem.t00z.a3d_0p50.f000.grib2 --fv3-prefix /path/to/coldstarts/input/gfs_data.`

The above will interpolate the fields from the specified GRIB file and place the fields in a modified set of FV3 cold start ICs in gfs_data.tile{tile}.nc
