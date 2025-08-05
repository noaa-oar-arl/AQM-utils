# AQM-utils

AQM-utils is a collection of utilities and source code for Air Composition Modeling data processing, analysis, and workflow automation. It supports building, running, and post-processing AQM outputs on multiple platforms (NOAA Tier 1 and Tier 2 machines) and includes both Python and Fortran tools for scientific and operational use.

### Available Utilities

**Python Utilities** (located in `python_utils/`):
- `add_aqm_ics.py`: Add initial conditions to AQM files.
- `append-chem-lbc.py`: Append chemical lateral boundary conditions.
- `intpl_ESMPy_QC_allspecies.py`: Interpolate fields using ESMpy for all species.
- `RAVE_remake.allspecies.aqmna13km.g793.py`: RAVE post-processing for 13km grid.
- `RAVE_remake.allspecies.aqmna9km.g1144.py`: RAVE post-processing for 9km grid.
- `stack-pt-merge.py`: Merge point source stacks.
- Subdirectories for specialized utilities:
  - `gefsaero_grib_to_cube/`: Interpolate GRIB2 GEFS-Aerosols fields to FV3 cold start files
  - `gdas-cmaqprep`: A Python tool for processing GDAS (Global Data Assimilation System) data for use with the CMAQ (Community Multiscale Air Quality) modeling system.

**Fortran Source Modules** (located in `sorc/`):
- `aqm_bias_correct.fd/`: Bias correction routines.
- `aqm_bias_interpolate.fd/`: Interpolation routines for bias correction.
- `aqm_post_bias_cor_grib2.fd/`: Post-processing for bias-corrected GRIB2 files.
- `aqm_post_grib2.fd/`: General post-processing for GRIB2 files.
- `aqm_post_maxi_bias_cor_grib2.fd/`: Maxi bias correction post-processing.
- `aqm_post_maxi_grib2.fd/`: Maxi post-processing for GRIB2 files.
- `convert_airnow_csv.fd/`: Convert AirNow CSV data for model use.
- `gefs2lbcs_para.fd/`: GEFS to LBCs conversion routines.

Additional configuration files, modulefiles, and parameter files are provided to support building and running the utilities on supported platforms.

## Build

### Use the `build.sh`

AQM-utils provides a `build.sh` script to automate the build process for supported platforms.

#### Basic Usage

1. Make sure you have the required modules and environment variables set for your platform (see below).
2. Run the build script from the repository root:
   ```bash
   ./build.sh
   ```

What the Script Does
- Loads the appropriate modulefile for your platform
- Creates a build directory if needed
- Runs CMake and Make with recommended options
- Installs binaries to the correct location

#### `build.sh` options

```bash
./build.sh -h
Start ... Thu Jul 24 12:09:37 EDT 2025
WARNING: UNKNOWN PLATFORM

Usage: ./build.sh -p <prefix> | -t <target> -h

  -p  installation prefix <prefix>    DEFAULT: <none>
  -t  target to build for <target>    DEFAULT: UNKNOWN
  -c  additional CMake options        DEFAULT: <none>
  -v  build with verbose output       DEFAULT: NO
  -f  force a clean build             DEFAULT: NO
  -h  display this message and quit
```

Manual Build Steps (if not using build.sh)
You can still build manually using the steps below for each platform.


#### Manual Building
```
cd AQM-utils
source versions/build.ver.wcoss2
module use modulefiles
module load build_[machine].[compiler]
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=.. -DCMAKE_INSTALL_BINDIR=exec
make -j2
where `[machine]` is `hera` for example.
```

## Run pre-commit
## Pre-commit Setup and Usage

This repository uses [pre-commit](https://pre-commit.com/) to automatically check and format code before each commit.

### Installation

1. Install pre-commit (if not already installed):
   ```bash
   pip install pre-commit
   ```
2. Install the hooks:
   ```bash
   pre-commit install
   ```

### Usage

- To run pre-commit checks on all files:
  ```bash
  pre-commit run --all-files
  ```
- Hooks will run automatically on staged files during each commit.
- If hooks modify files, re-stage and commit again until no changes are made.
