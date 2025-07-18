#!/usr/bin/env python3
"""
Script to generate GCAFS initial conditions from GRIB2 files and FV3 cold start files.

This script reads GRIB2 chemical species data and FV3 cubed-sphere cold start files
to generate initial conditions for the GCAFS (Global Chemistry Aerosol Forecasting System).
"""

import os
import sys
import argparse
import numpy as np
from pathlib import Path
from read_aerosols_grib2 import read_aerosol_species_from_grib2

def generate_gcafs_ics(grib_file, fv3_prefix, output_dir=None):
    """
    Generate GCAFS initial conditions from GRIB2 file and FV3 cold start files.
    
    Parameters:
    -----------
    grib_file : str
        Path to input GRIB2 file containing chemical species data
    fv3_prefix : str
        Path prefix for FV3 cold start files (should allow reading all 6 tiles)
        Example: "/path/to/fv3_data/gfs_data.tile" (will read gfs_data.tile1.nc through gfs_data.tile6.nc)
    output_dir : str, optional
        Directory to write output files. If None, uses current directory.
    
    Returns:
    --------
    bool
        True if successful, False otherwise
    """
    print(f"Starting GCAFS initial conditions generation...")
    print(f"Input GRIB file: {grib_file}")
    print(f"FV3 file prefix: {fv3_prefix}")
    
    # Validate input files exist
    if not os.path.exists(grib_file):
        print(f"ERROR: GRIB file not found: {grib_file}")
        return False
    
    # Check for all 6 FV3 tile files
    fv3_files = []
    for tile in range(1, 7):
        tile_file = f"{fv3_prefix}tile{tile}.nc"
        if not os.path.exists(tile_file):
            print(f"ERROR: FV3 tile file not found: {tile_file}")
            return False
        fv3_files.append(tile_file)
    
    print(f"Found all 6 FV3 tile files: {fv3_files}")
    
    # Set output directory
    if output_dir is None:
        output_dir = os.getcwd()
    
    # TODO: Implement the actual processing logic
    # This function will eventually:
    # 1. Read GRIB2 chemical species data
    # 2. Read FV3 cubed-sphere grid information from tile files
    # 3. Interpolate chemical species from regular grid to cubed-sphere
    # 4. Write out GCAFS initial condition files
    
    print("Processing logic to be implemented...")

    print(f"Extracting aerosol species data from GRIB2 file: {grib_file}")
    aerosol_data = read_aerosol_species_from_grib2(grib_file)
    if aerosol_data is None:
        print("ERROR: Failed to read aerosol species data from GRIB2 file.")
        return False
    print(f"Successfully extracted aerosol species data: {aerosol_data.keys()}")
    grid_info_grib = aerosol_data['_grid_info']
    print(grid_info_grib)
    levels_grib = aerosol_data['_levels']
    print(f"Levels in GRIB file: {levels_grib}")

    return True


def main():
    """Main function to parse arguments and run the GCAFS IC generation."""
    parser = argparse.ArgumentParser(
        description="Generate GCAFS initial conditions from GRIB2 and FV3 files"
    )
    
    parser.add_argument(
        "--grib-file",
        type=str,
        required=True,
        help="Path to input GRIB2 file containing chemical species data"
    )
    
    parser.add_argument(
        "--fv3-prefix",
        type=str,
        required=True,
        help="Path prefix for FV3 cold start files (e.g., /path/to/gfs_data.tile)"
    )
    
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Output directory for GCAFS initial condition files (default: current directory)"
    )
    
    args = parser.parse_args()
    
    # Run the IC generation
    success = generate_gcafs_ics(
        grib_file=args.grib_file,
        fv3_prefix=args.fv3_prefix,
        output_dir=args.output_dir
    )
    
    if success:
        print("GCAFS initial conditions generation completed successfully!")
        sys.exit(0)
    else:
        print("GCAFS initial conditions generation failed!")
        sys.exit(1)


if __name__ == "__main__":
    # Example usage (uncomment to test)
    # grib_file = "/scratch3/NCEPDEV/da/Cory.R.Martin/sample_files_for_ics/gefs.chem.t00z.a3d_0p50.f000.grib2"
    # fv3_prefix = "/scratch3/NCEPDEV/da/Cory.R.Martin/20241214.030000.sfc_data.tile"
    # generate_gcafs_ics(grib_file, fv3_prefix)
    
    main()