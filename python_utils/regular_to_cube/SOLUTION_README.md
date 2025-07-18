# Solution for dust1 Vertical Interpolation Error

## Problem
You encountered this error:
```
ValueError: Field dust1 has shape (64, 192, 192), expected (128, 192, 192) we need to interpolate in the vertical
```

This error occurs because your `dust1` field has 64 vertical levels, but the target FV3 tile file expects 128 levels.

## Solution
I've updated the `fv3_cold_starts.py` file to automatically handle vertical interpolation when there's a level mismatch. The solution includes:

### 1. New Vertical Interpolation Function
Added `interpolate_vertical_levels()` function that:
- Uses scipy.interpolate.interp1d for vertical interpolation
- Supports linear, cubic, and nearest neighbor interpolation methods
- Interpolates each horizontal grid point independently
- Preserves the vertical structure of the data

### 2. Updated Main Function
Modified `add_3d_fields_to_fv3_tile()` to:
- Detect when field dimensions don't match (specifically vertical levels)
- Automatically perform vertical interpolation from source to target levels
- Provide informative logging about the interpolation process
- Only error if horizontal dimensions don't match (which would be a real error)

### 3. Key Features
- **Automatic Detection**: Automatically detects level mismatches
- **Preserves Data Structure**: Maintains the physical relationships in the data
- **Flexible**: Supports different interpolation methods
- **Robust**: Handles edge cases and provides clear error messages
- **Backwards Compatible**: Still works for fields with correct dimensions

## Usage

### Before (would cause error):
```python
field_data = {
    'dust1': dust1_array_64_levels,  # Shape: (64, 192, 192) - ERROR!
}
```

### After (automatically handles interpolation):
```python
field_data = {
    'dust1': dust1_array_64_levels,  # Shape: (64, 192, 192) - Will be interpolated to (128, 192, 192)
}

# This will now work without errors
output_file = add_3d_fields_to_fv3_tile(
    tile_file_path="your_file.nc",
    field_data=field_data
)
```

## What Happens Now
1. Function detects that dust1 has 64 levels instead of expected 128
2. Logs a warning message about performing vertical interpolation
3. Calls `interpolate_vertical_levels()` to interpolate from 64 to 128 levels
4. Continues processing with the interpolated 128-level field
5. Successfully adds the field to your FV3 tile file

## Files Modified
- `fv3_cold_starts.py`: Added vertical interpolation functionality
- Added `scipy` import for interpolation functions

## Test Files Created
- `test_vertical_interp.py`: Demonstrates the interpolation functionality
- `example_usage.py`: Shows how to use the updated function

The solution maintains the physical integrity of your dust data while ensuring compatibility with the FV3 tile file structure.
