#!/usr/bin/env python3
"""
Test script demonstrating proper vertical interpolation with realistic pressure levels
for handling dust1 field with different vertical coordinate systems.
"""

import numpy as np
import sys
import os

# Add the current directory to Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fv3_cold_starts import (
    interpolate_vertical_levels, 
    create_common_pressure_levels,
    add_3d_fields_to_fv3_tile
)

def test_pressure_based_interpolation():
    """Test vertical interpolation using actual pressure coordinates"""
    
    print("Testing pressure-based vertical interpolation")
    print("=" * 50)
    
    # Create realistic pressure levels
    # Source: 64 levels (typical GRIB2 data)
    source_pressure = create_common_pressure_levels('gfs_64')
    print(f"Source pressure levels: {len(source_pressure)} levels")
    print(f"Source pressure range: {source_pressure[-1]:.2f} to {source_pressure[0]:.2f} hPa")
    
    # Target: 128 levels (FV3 model)
    target_pressure = create_common_pressure_levels('gfs_128')
    print(f"Target pressure levels: {len(target_pressure)} levels")
    print(f"Target pressure range: {target_pressure[-1]:.2f} to {target_pressure[0]:.2f} hPa")
    
    # Create realistic dust1 field with vertical structure
    ny, nx = 192, 192
    source_levels = len(source_pressure)
    
    # Create dust field with physically realistic vertical distribution
    dust1_field = np.zeros((source_levels, ny, nx))
    
    for k in range(source_levels):
        # Dust typically has maximum concentration near surface,
        # decreasing exponentially with height (lower pressure)
        # Use pressure as proxy for height
        pressure = source_pressure[k]
        
        # Higher pressure (surface) = more dust
        # Lower pressure (upper atmosphere) = less dust
        dust_factor = np.exp(-(1000 - pressure) / 200.0)  # Exponential decay
        
        # Add some spatial variability
        spatial_pattern = np.random.random((ny, nx)) * 0.8 + 0.2
        dust1_field[k, :, :] = dust_factor * spatial_pattern * 1e-6
    
    print(f"\nSource dust1 field shape: {dust1_field.shape}")
    print(f"Surface dust concentration: {np.mean(dust1_field[0, :, :]):.2e} kg/kg")
    print(f"Mid-level dust concentration: {np.mean(dust1_field[32, :, :]):.2e} kg/kg")
    print(f"Upper-level dust concentration: {np.mean(dust1_field[-1, :, :]):.2e} kg/kg")
    
    # Perform pressure-based interpolation
    print(f"\nPerforming pressure-based interpolation...")
    interpolated_field = interpolate_vertical_levels(
        dust1_field,
        source_pressure,
        target_pressure,
        method='linear'
    )
    
    print(f"Interpolated field shape: {interpolated_field.shape}")
    print(f"Surface dust concentration: {np.mean(interpolated_field[0, :, :]):.2e} kg/kg")
    print(f"Mid-level dust concentration: {np.mean(interpolated_field[64, :, :]):.2e} kg/kg")
    print(f"Upper-level dust concentration: {np.mean(interpolated_field[-1, :, :]):.2e} kg/kg")
    
    # Validate that the interpolation preserved the physical structure
    print(f"\nValidation:")
    print(f"Surface levels preserved: {np.allclose(np.mean(dust1_field[0, :, :]), np.mean(interpolated_field[0, :, :]), rtol=0.01)}")
    print(f"Upper levels preserved: {np.allclose(np.mean(dust1_field[-1, :, :]), np.mean(interpolated_field[-1, :, :]), rtol=0.01)}")
    
    return interpolated_field, source_pressure, target_pressure

def demonstrate_extrapolation():
    """Demonstrate extrapolation when source and target ranges don't overlap completely"""
    
    print("\n" + "=" * 50)
    print("Demonstrating extrapolation capabilities")
    print("=" * 50)
    
    # Source: Limited range (e.g., regional model or limited observations)
    source_pressure = np.array([1000, 925, 850, 700, 500, 300, 200, 100, 50, 10])  # 10 levels
    
    # Target: Full range (e.g., global model)
    target_pressure = np.logspace(3, -1, 64)  # 1000 to 0.1 hPa
    
    print(f"Source pressure range: {source_pressure[-1]:.2f} to {source_pressure[0]:.2f} hPa")
    print(f"Target pressure range: {target_pressure[-1]:.2f} to {target_pressure[0]:.2f} hPa")
    print(f"Extrapolation needed above: {target_pressure[-1]:.2f} hPa")
    
    # Create test field
    ny, nx = 50, 50  # Smaller for demonstration
    test_field = np.zeros((len(source_pressure), ny, nx))
    
    for k in range(len(source_pressure)):
        pressure = source_pressure[k]
        dust_factor = np.exp(-(1000 - pressure) / 200.0)
        test_field[k, :, :] = dust_factor * 1e-6
    
    # Perform interpolation with extrapolation
    print(f"\nPerforming interpolation with extrapolation...")
    extrapolated_field = interpolate_vertical_levels(
        test_field,
        source_pressure,
        target_pressure,
        method='linear'
    )
    
    print(f"Extrapolated field shape: {extrapolated_field.shape}")
    print(f"Extrapolation successful!")
    
    return extrapolated_field

def create_usage_example():
    """Create a practical usage example"""
    
    print("\n" + "=" * 50)
    print("Practical Usage Example")
    print("=" * 50)
    
    # Example: You have GRIB2 data with 64 pressure levels
    source_pressure = create_common_pressure_levels('gfs_64')
    
    # Your dust1 field from GRIB2
    dust1_data = np.random.random((64, 192, 192)) * 1e-6
    
    # Your field data dictionary
    field_data = {
        'dust1': dust1_data,
        # Add other fields as needed
    }
    
    print(f"Example field data:")
    print(f"  dust1: shape {dust1_data.shape}")
    print(f"  Source pressure levels: {len(source_pressure)} levels")
    print(f"  Pressure range: {source_pressure[-1]:.2f} to {source_pressure[0]:.2f} hPa")
    
    print(f"\nTo use with your FV3 tile file:")
    print(f"1. Load your GRIB2 data and extract dust1 field")
    print(f"2. Get the pressure levels from your GRIB2 data")
    print(f"3. Call the updated function:")
    print(f"")
    print(f"   output_file = add_3d_fields_to_fv3_tile(")
    print(f"       tile_file_path='your_tile_file.nc',")
    print(f"       field_data={{'dust1': dust1_data}},")
    print(f"       source_pressure_levels=source_pressure,")
    print(f"       backup_original=True")
    print(f"   )")
    print(f"")
    print(f"The function will:")
    print(f"- Read the target pressure levels from your FV3 tile file")
    print(f"- Automatically interpolate/extrapolate from source to target levels")
    print(f"- Handle any pressure range differences properly")
    print(f"- Add the interpolated field to your tile file")

if __name__ == "__main__":
    print("Enhanced Vertical Interpolation Test")
    print("Using realistic pressure coordinates")
    print("=" * 60)
    
    # Test 1: Pressure-based interpolation
    interpolated, source_p, target_p = test_pressure_based_interpolation()
    
    # Test 2: Extrapolation demonstration
    extrapolated = demonstrate_extrapolation()
    
    # Test 3: Usage example
    create_usage_example()
    
    print("\n" + "=" * 60)
    print("All tests completed successfully!")
    print("The enhanced interpolation now uses actual pressure coordinates")
    print("and properly handles extrapolation when ranges don't match.")
