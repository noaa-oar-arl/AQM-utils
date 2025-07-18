#!/usr/bin/env python3
"""
Debug script to test the FV3 tile modification function
"""

import numpy as np
import sys
import os
import xarray as xr
from pathlib import Path

# Add the current directory to Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def debug_fv3_tile():
    """Debug the FV3 tile file modification process"""
    
    tile_file = '/scratch3/NCEPDEV/da/Cory.R.Martin/july2025/sample_files_for_ics/gdas.20250701/00/model/atmos/input/gfs_data.tile2.nc'
    
    print("=== Debug: FV3 Tile File Analysis ===")
    
    # Check if file exists
    if not Path(tile_file).exists():
        print(f"❌ File does not exist: {tile_file}")
        return
    
    print(f"✅ File exists: {tile_file}")
    
    # Examine the file structure
    try:
        with xr.open_dataset(tile_file) as ds:
            print(f"\n=== File Structure ===")
            print(f"Dimensions: {dict(ds.dims)}")
            print(f"Number of data variables: {len(ds.data_vars)}")
            print(f"Number of coordinates: {len(ds.coords)}")
            
            print(f"\n=== Data Variables (first 10) ===")
            for i, var in enumerate(list(ds.data_vars.keys())[:10]):
                shape = ds[var].shape
                dtype = ds[var].dtype
                print(f"  {var}: {shape} ({dtype})")
            
            print(f"\n=== Coordinates ===")
            for coord in ds.coords:
                shape = ds[coord].shape
                dtype = ds[coord].dtype
                print(f"  {coord}: {shape} ({dtype})")
            
            # Check for required coordinates
            print(f"\n=== Coordinate Check ===")
            required_coords = ['geolon', 'geolat', 'lev', 'lat', 'lon']
            for coord in required_coords:
                if coord in ds.coords:
                    print(f"  ✅ {coord}: present")
                else:
                    print(f"  ❌ {coord}: missing")
            
            # Test creating a simple field
            print(f"\n=== Testing Field Addition ===")
            nlev = ds.dims['lev']
            nlat = ds.dims['lat'] 
            nlon = ds.dims['lon']
            
            print(f"Target dimensions: {nlev} levels, {nlat}x{nlon} grid")
            
            # Create test data
            test_data = np.random.random((nlev, nlat, nlon)) * 1e-6
            print(f"Test data shape: {test_data.shape}")
            
            # Try to create a new variable manually
            print(f"\n=== Manual Variable Creation Test ===")
            
            # Make a copy of the dataset
            ds_out = ds.copy(deep=True)
            
            # Create the DataArray for the test field
            test_field = xr.DataArray(
                test_data,
                dims=['lev', 'lat', 'lon'],
                coords={
                    'lev': ds['lev'] if 'lev' in ds.coords else ds.coords['lev'],
                    'lat': ds['lat'] if 'lat' in ds.coords else ds.coords['lat'],
                    'lon': ds['lon'] if 'lon' in ds.coords else ds.coords['lon']
                },
                name='test_dust1'
            )
            
            # Add metadata
            test_field.attrs = {
                'units': 'kg kg-1',
                'long_name': 'test dust mixing ratio',
                'coordinates': 'geolon geolat'
            }
            
            # Add to dataset
            ds_out['test_dust1'] = test_field
            
            print(f"✅ Successfully created test field")
            print(f"Original variables: {len(ds.data_vars)}")
            print(f"Modified variables: {len(ds_out.data_vars)}")
            
            # Check if we can write to a test file
            test_output = 'test_tile_output.nc'
            print(f"\n=== Testing File Write ===")
            
            try:
                # Create backup first
                backup_file = f"{tile_file}.backup_debug"
                if not Path(backup_file).exists():
                    import shutil
                    shutil.copy2(tile_file, backup_file)
                    print(f"✅ Created backup: {backup_file}")
                
                # Write test file
                ds_out.to_netcdf(test_output)
                print(f"✅ Successfully wrote test file: {test_output}")
                
                # Verify the test file
                with xr.open_dataset(test_output) as test_ds:
                    if 'test_dust1' in test_ds.data_vars:
                        print(f"✅ Test field present in output file")
                        print(f"✅ Output file has {len(test_ds.data_vars)} variables (original: {len(ds.data_vars)})")
                        
                        # Check if original variables are preserved
                        missing_vars = []
                        for var in ds.data_vars:
                            if var not in test_ds.data_vars:
                                missing_vars.append(var)
                        
                        if missing_vars:
                            print(f"❌ Missing variables in output: {missing_vars}")
                        else:
                            print(f"✅ All original variables preserved")
                    else:
                        print(f"❌ Test field missing from output file")
                
                # Clean up
                if Path(test_output).exists():
                    Path(test_output).unlink()
                    print(f"🧹 Cleaned up test file")
                    
            except Exception as write_error:
                print(f"❌ Write test failed: {write_error}")
                import traceback
                traceback.print_exc()
                
    except Exception as e:
        print(f"❌ Error reading file: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    debug_fv3_tile()
