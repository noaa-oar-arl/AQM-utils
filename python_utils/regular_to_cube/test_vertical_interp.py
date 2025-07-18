#!/usr/bin/env python3
"""
Test script demonstrating the vertical interpolation functionality
for handling dust1 field with 64 levels instead of expected 128 levels.
"""

import numpy as np
import sys
import os

# Add the current directory to Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fv3_cold_starts import interpolate_vertical_levels

def test_vertical_interpolation():
    """Test the vertical interpolation function"""
    
    # Create test data similar to your dust1 field
    # Original shape: (64, 192, 192) -> Target shape: (128, 192, 192)
    source_levels = 64
    target_levels = 128
    ny, nx = 192, 192
    
    print(f"Testing vertical interpolation from {source_levels} to {target_levels} levels")
    print(f"Horizontal grid: {ny} x {nx}")
    
    # Create synthetic dust1 data with realistic values
    # Dust typically decreases with height, so we'll create that pattern
    test_data = np.zeros((source_levels, ny, nx))
    
    for k in range(source_levels):
        # Exponential decay with height (higher k = higher altitude = less dust)
        height_factor = np.exp(-k / 20.0)  # Decay constant
        # Add some random spatial variability
        spatial_pattern = np.random.random((ny, nx)) * 1e-6
        test_data[k, :, :] = height_factor * spatial_pattern
    
    print(f"Original data shape: {test_data.shape}")
    print(f"Original data range: {np.min(test_data):.2e} to {np.max(test_data):.2e}")
    
    # Perform vertical interpolation
    interpolated_data = interpolate_vertical_levels(
        test_data, 
        source_levels, 
        target_levels, 
        method='linear'
    )
    
    print(f"Interpolated data shape: {interpolated_data.shape}")
    print(f"Interpolated data range: {np.min(interpolated_data):.2e} to {np.max(interpolated_data):.2e}")
    
    # Verify the interpolation preserved the general structure
    print("\nVerification:")
    print(f"Surface level (k=0) original: {np.mean(test_data[0, :, :]):.2e}")
    print(f"Surface level (k=0) interpolated: {np.mean(interpolated_data[0, :, :]):.2e}")
    
    print(f"Mid-level original (k={source_levels//2}): {np.mean(test_data[source_levels//2, :, :]):.2e}")
    print(f"Mid-level interpolated (k={target_levels//2}): {np.mean(interpolated_data[target_levels//2, :, :]):.2e}")
    
    print(f"Top level original (k={source_levels-1}): {np.mean(test_data[-1, :, :]):.2e}")
    print(f"Top level interpolated (k={target_levels-1}): {np.mean(interpolated_data[-1, :, :]):.2e}")
    
    return interpolated_data

def create_example_field_data():
    """Create example field data that would cause the original ValueError"""
    
    # This simulates your dust1 field with 64 levels instead of 128
    dust1_64_levels = np.random.random((64, 192, 192)) * 1e-6
    
    # Other dust fields that already have the correct dimensions
    dust2_128_levels = np.random.random((128, 192, 192)) * 1e-6
    
    field_data = {
        'dust1': dust1_64_levels,  # This will trigger vertical interpolation
        'dust2': dust2_128_levels  # This will pass through normally
    }
    
    print("\nExample field data created:")
    for field_name, field_array in field_data.items():
        print(f"{field_name}: shape {field_array.shape}")
    
    return field_data

if __name__ == "__main__":
    print("Testing vertical interpolation for dust1 field")
    print("=" * 50)
    
    # Test the interpolation function
    interpolated = test_vertical_interpolation()
    
    print("\n" + "=" * 50)
    print("Example of field data that would cause the original error:")
    
    # Show example of what your field data might look like
    example_fields = create_example_field_data()
    
    print("\nWith the updated fv3_cold_starts.py, the dust1 field with 64 levels")
    print("will be automatically interpolated to 128 levels, resolving the ValueError.")
