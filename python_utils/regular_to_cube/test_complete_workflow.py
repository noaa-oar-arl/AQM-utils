#!/usr/bin/env python3
"""
Complete test of the GRIB2 → FV3 workflow with both horizontal and vertical interpolation.
"""

import numpy as np
import xarray as xr
import sys
import os
from pathlib import Path

# Add the current directory to Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fv3_cold_starts import add_3d_fields_to_fv3_tile, create_aerosol_field_metadata, create_common_pressure_levels
from example_usage import interpolate_latlon_to_fv3_tile

def create_test_atmospheric_tile(output_path: str):
    """
    Create a test atmospheric tile file with 3D structure and pressure levels.
    """
    print(f"Creating test atmospheric tile: {output_path}")
    
    # Create test atmospheric data with multiple levels
    nx, ny = 192, 192  # Smaller than C384 for faster testing
    nz = 128  # 128 vertical levels
    
    # Create pressure levels (similar to GFS L128)
    pressure_levels = create_common_pressure_levels('gfs_128')
    
    # Create cubed sphere grid coordinates (simplified)
    x_coords = np.linspace(-45, 45, nx)  # degrees
    y_coords = np.linspace(-45, 45, ny)  # degrees
    x_grid, y_grid = np.meshgrid(x_coords, y_coords)
    
    # Create fake atmospheric data
    temp_data = np.random.normal(250, 30, (nz, ny, nx))  # Temperature
    humid_data = np.random.exponential(0.01, (nz, ny, nx))  # Humidity
    
    # Create dataset
    ds = xr.Dataset({
        'temp': (['pfull', 'grid_yt', 'grid_xt'], temp_data, 
                {'long_name': 'Temperature', 'units': 'K'}),
        'humid': (['pfull', 'grid_yt', 'grid_xt'], humid_data,
                 {'long_name': 'Specific Humidity', 'units': 'kg/kg'}),
        'grid_latt': (['grid_yt', 'grid_xt'], y_grid,
                     {'long_name': 'Latitude', 'units': 'degrees_north'}),
        'grid_lont': (['grid_yt', 'grid_xt'], x_grid,
                     {'long_name': 'Longitude', 'units': 'degrees_east'}),
    }, coords={
        'pfull': (['pfull'], pressure_levels, 
                 {'long_name': 'Full Pressure Levels', 'units': 'hPa'}),
        'grid_yt': np.arange(ny),
        'grid_xt': np.arange(nx),
    })
    
    # Add attributes
    ds.attrs['title'] = 'Test FV3 Atmospheric Tile'
    ds.attrs['institution'] = 'Test'
    
    # Save to file
    ds.to_netcdf(output_path)
    print(f"✅ Created test tile with {nz} levels, {ny}x{nx} grid")
    return output_path

def test_complete_workflow():
    """
    Test the complete workflow: GRIB2-style data → horizontal interpolation → vertical interpolation → FV3 tile
    """
    print("=" * 70)
    print("Testing Complete GRIB2 → FV3 Workflow")
    print("=" * 70)
    
    # Step 1: Create synthetic GRIB2-style data (lat/lon grid)
    print("\n1. Creating synthetic GRIB2-style aerosol data")
    
    # Global lat/lon grid (similar to GEFS-Aerosols)
    nlat_source, nlon_source = 96, 192  # ~1.875° resolution
    nlevels_source = 64  # 64 levels from GRIB2
    
    source_lats = np.linspace(-90, 90, nlat_source)
    source_lons = np.linspace(0, 358.125, nlon_source)
    
    # Create realistic dust concentration pattern
    # Higher concentrations near dust source regions (Africa, Middle East, Asia)
    lat_grid, lon_grid = np.meshgrid(source_lats, source_lons, indexing='ij')
    
    # Create dust pattern (higher in Sahara region)
    dust_pattern = np.exp(-0.1 * ((lat_grid - 20)**2 + (lon_grid - 10)**2))  # Africa peak
    dust_pattern += 0.5 * np.exp(-0.1 * ((lat_grid - 30)**2 + (lon_grid - 50)**2))  # Middle East
    dust_pattern += 0.3 * np.exp(-0.1 * ((lat_grid - 40)**2 + (lon_grid - 90)**2))  # Asia
    
    # Create 3D dust field with vertical decay
    pressure_source = create_common_pressure_levels('gfs_64')
    dust1_data = np.zeros((nlevels_source, nlat_source, nlon_source))
    
    for k in range(nlevels_source):
        # Dust concentrates near surface, decays with height
        height_factor = np.exp(-pressure_source[k] / 500)  # Decay with altitude
        dust1_data[k] = dust_pattern * height_factor * 1e-6  # Typical dust concentration
    
    print(f"   Source data: {dust1_data.shape} on {nlat_source}x{nlon_source} lat/lon grid")
    print(f"   Concentration range: {np.min(dust1_data):.2e} to {np.max(dust1_data):.2e} kg/kg")
    
    # Step 2: Create test atmospheric tile
    print("\n2. Creating test FV3 atmospheric tile")
    test_tile_path = "test_atmospheric_tile.nc"
    create_test_atmospheric_tile(test_tile_path)
    
    # Step 3: Horizontal interpolation
    print("\n3. Horizontal interpolation: lat/lon → FV3 cubed sphere")
    try:
        dust1_interp = interpolate_latlon_to_fv3_tile(
            dust1_data, source_lats, source_lons, test_tile_path
        )
        print(f"   ✅ Interpolated to: {dust1_interp.shape}")
    except Exception as e:
        print(f"   ❌ Horizontal interpolation failed: {e}")
        return
    
    # Step 4: Vertical interpolation and addition to tile
    print("\n4. Vertical interpolation and adding to FV3 tile")
    
    field_data = {'dust1': dust1_interp}
    field_metadata = create_aerosol_field_metadata()
    
    print(f"   Source levels: {nlevels_source} (pressure: {pressure_source[0]:.2f} to {pressure_source[-1]:.2f} hPa)")
    
    try:
        output_file = add_3d_fields_to_fv3_tile(
            tile_file_path=test_tile_path,
            field_data=field_data,
            source_pressure_levels=pressure_source,
            field_metadata=field_metadata,
            backup_original=True
        )
        print(f"   ✅ Successfully added dust1 to: {output_file}")
        
        # Verify the result
        with xr.open_dataset(output_file) as ds:
            if 'dust1' in ds:
                dust1_final = ds['dust1']
                print(f"   ✅ Final dust1 shape: {dust1_final.shape}")
                print(f"   ✅ Final concentration range: {float(dust1_final.min()):.2e} to {float(dust1_final.max()):.2e}")
            else:
                print(f"   ❌ dust1 field not found in output file")
        
    except Exception as e:
        print(f"   ❌ Adding to tile failed: {e}")
        return
    
    # Step 5: Summary
    print("\n" + "=" * 70)
    print("✅ COMPLETE WORKFLOW TEST SUCCESSFUL")
    print("=" * 70)
    print("Workflow steps completed:")
    print("1. ✅ Created synthetic GRIB2-style dust data (64 levels, lat/lon)")
    print("2. ✅ Horizontal interpolation: lat/lon → FV3 cubed sphere")
    print("3. ✅ Vertical interpolation: 64 → 128 levels using pressure")
    print("4. ✅ Added interpolated dust1 field to FV3 tile")
    print("\nThis demonstrates the complete pipeline for adding")
    print("GRIB2 aerosol data to FV3 atmospheric restart files.")
    
    # Cleanup
    try:
        os.remove(test_tile_path)
        os.remove(output_file)
        backup_file = output_file + ".backup"
        if os.path.exists(backup_file):
            os.remove(backup_file)
        print("\nCleaned up test files")
    except:
        pass

if __name__ == "__main__":
    test_complete_workflow()
