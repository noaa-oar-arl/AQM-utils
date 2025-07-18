#!/usr/bin/env python3
# fv3_cold_starts.py
# Functions to add 3D aerosol fields to FV3 tile files for cold starts

import numpy as np
import xarray as xr
import logging
from pathlib import Path
from typing import Dict, List, Optional, Union
import shutil

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def add_3d_fields_to_fv3_tile(tile_file_path: str, 
                              field_data: Dict[str, np.ndarray],
                              field_metadata: Optional[Dict[str, Dict]] = None,
                              output_path: Optional[str] = None,
                              backup_original: bool = True) -> str:
    """
    Add 3D fields to an FV3 tile file (e.g., gfs_data.tile1.nc).
    
    Args:
        tile_file_path (str): Path to the FV3 tile file
        field_data (Dict[str, np.ndarray]): Dictionary of field names and their 3D data arrays
                                          Expected shape: (nlev, ny, nx) for each field
        field_metadata (Dict[str, Dict], optional): Metadata for each field including:
                                                  - units: str
                                                  - long_name: str
                                                  - standard_name: str (optional)
        output_path (str, optional): Path for output file. If None, overwrites original
        backup_original (bool): Whether to create a backup of the original file
    
    Returns:
        str: Path to the output file
    
    Raises:
        FileNotFoundError: If the tile file doesn't exist
        ValueError: If field dimensions don't match tile dimensions
    """
    
    tile_path = Path(tile_file_path)
    if not tile_path.exists():
        raise FileNotFoundError(f"Tile file not found: {tile_file_path}")
    
    logger.info(f"Opening FV3 tile file: {tile_file_path}")
    
    # Open the tile file
    with xr.open_dataset(tile_file_path) as ds:
        # Get tile dimensions
        if 'pfull' in ds.dims:
            nlev = ds.dims['pfull']
        elif 'lev' in ds.dims:
            nlev = ds.dims['lev']
        else:
            # Try to infer from existing 3D variables
            for var in ds.data_vars:
                if len(ds[var].dims) == 3:
                    nlev = ds[var].shape[0]
                    break
            else:
                raise ValueError("Cannot determine number of levels from tile file")
        
        # Get horizontal dimensions
        if 'grid_yt' in ds.dims and 'grid_xt' in ds.dims:
            ny, nx = ds.dims['grid_yt'], ds.dims['grid_xt']
        elif 'lat' in ds.dims and 'lon' in ds.dims:
            ny, nx = ds.dims['lat'], ds.dims['lon']
        else:
            # Try to infer from existing 2D/3D variables
            for var in ds.data_vars:
                if len(ds[var].dims) >= 2:
                    var_dims = ds[var].dims
                    if len(var_dims) == 3:
                        ny, nx = ds[var].shape[1], ds[var].shape[2]
                    else:
                        ny, nx = ds[var].shape[0], ds[var].shape[1]
                    break
            else:
                raise ValueError("Cannot determine horizontal dimensions from tile file")
        
        logger.info(f"Tile dimensions: {nlev} levels, {ny}x{nx} horizontal grid")
        
        # Create a copy of the dataset to modify
        ds_out = ds.copy(deep=True)
        
        # Add each field to the dataset
        for field_name, field_array in field_data.items():
            logger.info(f"Adding field: {field_name}")
            
            # Validate field dimensions
            if field_array.shape != (nlev, ny, nx):
                raise ValueError(f"Field {field_name} has shape {field_array.shape}, "
                               f"expected ({nlev}, {ny}, {nx})")
            
            # Determine coordinate names based on what's available in the dataset
            if 'pfull' in ds.coords:
                lev_coord = 'pfull'
            elif 'lev' in ds.coords:
                lev_coord = 'lev'
            else:
                # Create a level coordinate if it doesn't exist
                lev_coord = 'pfull'
                ds_out[lev_coord] = (lev_coord, np.arange(1, nlev + 1))
                ds_out[lev_coord].attrs = {'units': 'level', 'long_name': 'pressure level'}
            
            if 'grid_yt' in ds.coords and 'grid_xt' in ds.coords:
                lat_coord, lon_coord = 'grid_yt', 'grid_xt'
            elif 'lat' in ds.coords and 'lon' in ds.coords:
                lat_coord, lon_coord = 'lat', 'lon'
            else:
                # Use existing coordinate names from a 3D variable
                for var in ds.data_vars:
                    if len(ds[var].dims) == 3:
                        lat_coord, lon_coord = ds[var].dims[1], ds[var].dims[2]
                        break
                else:
                    raise ValueError("Cannot determine coordinate names")
            
            # Create the DataArray for the field
            field_da = xr.DataArray(
                field_array,
                dims=[lev_coord, lat_coord, lon_coord],
                coords={
                    lev_coord: ds_out[lev_coord],
                    lat_coord: ds_out[lat_coord],
                    lon_coord: ds_out[lon_coord]
                },
                name=field_name
            )
            
            # Add metadata if provided
            if field_metadata and field_name in field_metadata:
                field_da.attrs.update(field_metadata[field_name])
            else:
                # Default metadata
                field_da.attrs = {
                    'units': 'kg kg-1',
                    'long_name': f'{field_name} mass mixing ratio'
                }
            
            # Add the field to the dataset
            ds_out[field_name] = field_da
            
            logger.info(f"Added {field_name}: shape {field_array.shape}, "
                       f"min={np.min(field_array):.2e}, max={np.max(field_array):.2e}")
    
    # Determine output path
    if output_path is None:
        output_path = tile_file_path
        # Create backup if requested
        if backup_original:
            backup_path = f"{tile_file_path}.backup"
            if not Path(backup_path).exists():
                logger.info(f"Creating backup: {backup_path}")
                shutil.copy2(tile_file_path, backup_path)
    
    # Write the modified dataset
    logger.info(f"Writing modified dataset to: {output_path}")
    ds_out.to_netcdf(output_path)
    
    logger.info(f"Successfully added {len(field_data)} fields to {output_path}")
    
    return output_path


def add_aerosol_fields_to_all_tiles(tile_directory: str,
                                   aerosol_data: Dict[str, np.ndarray],
                                   field_metadata: Optional[Dict[str, Dict]] = None,
                                   tile_pattern: str = "gfs_data.tile*.nc",
                                   output_directory: Optional[str] = None,
                                   backup_original: bool = True) -> List[str]:
    """
    Add aerosol fields to all FV3 tile files in a directory.
    
    Args:
        tile_directory (str): Directory containing FV3 tile files
        aerosol_data (Dict[str, np.ndarray]): Dictionary of aerosol field names and their 3D data
        field_metadata (Dict[str, Dict], optional): Metadata for each field
        tile_pattern (str): Glob pattern to match tile files
        output_directory (str, optional): Directory for output files. If None, overwrites originals
        backup_original (bool): Whether to create backups of original files
    
    Returns:
        List[str]: List of paths to modified tile files
    
    Raises:
        FileNotFoundError: If tile directory doesn't exist
        ValueError: If no tile files found matching pattern
    """
    
    tile_dir = Path(tile_directory)
    if not tile_dir.exists():
        raise FileNotFoundError(f"Tile directory not found: {tile_directory}")
    
    # Find all tile files matching the pattern
    tile_files = list(tile_dir.glob(tile_pattern))
    if not tile_files:
        raise ValueError(f"No tile files found matching pattern '{tile_pattern}' in {tile_directory}")
    
    logger.info(f"Found {len(tile_files)} tile files to process")
    
    output_files = []
    
    for tile_file in sorted(tile_files):
        logger.info(f"Processing tile file: {tile_file.name}")
        
        # Determine output path
        if output_directory:
            output_dir = Path(output_directory)
            output_dir.mkdir(parents=True, exist_ok=True)
            output_path = output_dir / tile_file.name
        else:
            output_path = None
        
        # Add fields to this tile
        try:
            output_file = add_3d_fields_to_fv3_tile(
                str(tile_file),
                aerosol_data,
                field_metadata,
                str(output_path) if output_path else None,
                backup_original
            )
            output_files.append(output_file)
            logger.info(f"Successfully processed {tile_file.name}")
            
        except Exception as e:
            logger.error(f"Failed to process {tile_file.name}: {e}")
            continue
    
    logger.info(f"Successfully processed {len(output_files)} out of {len(tile_files)} tile files")
    
    return output_files


def create_aerosol_field_metadata() -> Dict[str, Dict]:
    """
    Create metadata dictionary for aerosol fields matching FV3 tracer names.
    
    Returns:
        Dict[str, Dict]: Metadata for each aerosol field using FV3 naming convention
    """
    
    metadata = {
        # Dust species (dust1-dust5)
        'dust1': {
            'units': 'kg kg-1',
            'long_name': 'dust mixing ratio bin 1 (0.2-2 μm)',
            'standard_name': 'mass_fraction_of_dust_dry_aerosol_particles_in_air'
        },
        'dust2': {
            'units': 'kg kg-1',
            'long_name': 'dust mixing ratio bin 2 (2-3.6 μm)',
            'standard_name': 'mass_fraction_of_dust_dry_aerosol_particles_in_air'
        },
        'dust3': {
            'units': 'kg kg-1',
            'long_name': 'dust mixing ratio bin 3 (3.6-6 μm)',
            'standard_name': 'mass_fraction_of_dust_dry_aerosol_particles_in_air'
        },
        'dust4': {
            'units': 'kg kg-1',
            'long_name': 'dust mixing ratio bin 4 (6-12 μm)',
            'standard_name': 'mass_fraction_of_dust_dry_aerosol_particles_in_air'
        },
        'dust5': {
            'units': 'kg kg-1',
            'long_name': 'dust mixing ratio bin 5 (12-20 μm)',
            'standard_name': 'mass_fraction_of_dust_dry_aerosol_particles_in_air'
        },
        # Sea salt species (seas1-seas5)
        'seas1': {
            'units': 'kg kg-1',
            'long_name': 'sea salt mixing ratio bin 1 (0.06-0.2 μm)',
            'standard_name': 'mass_fraction_of_sea_salt_dry_aerosol_particles_in_air'
        },
        'seas2': {
            'units': 'kg kg-1',
            'long_name': 'sea salt mixing ratio bin 2 (0.2-1 μm)',
            'standard_name': 'mass_fraction_of_sea_salt_dry_aerosol_particles_in_air'
        },
        'seas3': {
            'units': 'kg kg-1',
            'long_name': 'sea salt mixing ratio bin 3 (1-3 μm)',
            'standard_name': 'mass_fraction_of_sea_salt_dry_aerosol_particles_in_air'
        },
        'seas4': {
            'units': 'kg kg-1',
            'long_name': 'sea salt mixing ratio bin 4 (3-10 μm)',
            'standard_name': 'mass_fraction_of_sea_salt_dry_aerosol_particles_in_air'
        },
        'seas5': {
            'units': 'kg kg-1',
            'long_name': 'sea salt mixing ratio bin 5 (10-20 μm)',
            'standard_name': 'mass_fraction_of_sea_salt_dry_aerosol_particles_in_air'
        },
        # Sulfate and other species
        'so4': {
            'units': 'kg kg-1',
            'long_name': 'sulfate aerosol mixing ratio',
            'standard_name': 'mass_fraction_of_sulfate_dry_aerosol_particles_in_air'
        },
        'so2': {
            'units': 'kg kg-1',
            'long_name': 'sulfur dioxide mixing ratio',
            'standard_name': 'mass_fraction_of_sulfur_dioxide_in_air'
        },
        'dms': {
            'units': 'kg kg-1',
            'long_name': 'dimethyl sulfide mixing ratio',
            'standard_name': 'mass_fraction_of_dimethyl_sulfide_in_air'
        },
        'msa': {
            'units': 'kg kg-1',
            'long_name': 'methane sulfonic acid mixing ratio',
            'standard_name': 'mass_fraction_of_methane_sulfonic_acid_in_air'
        },
        # Black carbon species
        'bc1': {
            'units': 'kg kg-1',
            'long_name': 'hydrophobic black carbon aerosol mixing ratio',
            'standard_name': 'mass_fraction_of_black_carbon_dry_aerosol_particles_in_air'
        },
        'bc2': {
            'units': 'kg kg-1',
            'long_name': 'hydrophilic black carbon aerosol mixing ratio',
            'standard_name': 'mass_fraction_of_black_carbon_dry_aerosol_particles_in_air'
        },
        # Organic carbon species
        'oc1': {
            'units': 'kg kg-1',
            'long_name': 'hydrophobic organic carbon aerosol mixing ratio',
            'standard_name': 'mass_fraction_of_particulate_organic_matter_dry_aerosol_particles_in_air'
        },
        'oc2': {
            'units': 'kg kg-1',
            'long_name': 'hydrophilic organic carbon aerosol mixing ratio',
            'standard_name': 'mass_fraction_of_particulate_organic_matter_dry_aerosol_particles_in_air'
        },
        # Particulate matter (diagnostic fields)
        'pm25': {
            'units': 'kg kg-1',
            'long_name': 'particulate matter with diameter <= 2.5 μm',
            'standard_name': 'mass_fraction_of_pm2p5_ambient_aerosol_particles_in_air'
        },
        'pm10': {
            'units': 'kg kg-1',
            'long_name': 'particulate matter with diameter <= 10 μm',
            'standard_name': 'mass_fraction_of_pm10_ambient_aerosol_particles_in_air'
        }
    }
    
    return metadata


def example_usage():
    """
    Example of how to use the FV3 tile modification functions.
    """
    
    # Example paths
    tile_file = "/gpfs/f6/ira-sti/proj-shared/Cory.R.Martin/july2025/gcafs/gdas.init/output/gdas.20250701/00/model/atmos/input/gfs_data.tile1.nc"
    tile_directory = "/gpfs/f6/ira-sti/proj-shared/Cory.R.Martin/july2025/gcafs/gdas.init/output/gdas.20250701/00/model/atmos/input/"
    
    # Example aerosol data (replace with actual data from grib2_to_cube.py)
    # This would come from the GRIB2 reader output
    aerosol_data = {
        'dust1': np.random.random((64, 192, 192)) * 1e-6,  # Example shape and values
        'dust2': np.random.random((64, 192, 192)) * 1e-6,
        'dust3': np.random.random((64, 192, 192)) * 1e-6,
        'dust4': np.random.random((64, 192, 192)) * 1e-6,
        'dust5': np.random.random((64, 192, 192)) * 1e-6,
        'seas1': np.random.random((64, 192, 192)) * 1e-6,
        'seas2': np.random.random((64, 192, 192)) * 1e-6,
        'seas3': np.random.random((64, 192, 192)) * 1e-6,
        'seas4': np.random.random((64, 192, 192)) * 1e-6,
        'seas5': np.random.random((64, 192, 192)) * 1e-6,
        'so4': np.random.random((64, 192, 192)) * 1e-7,
        'so2': np.random.random((64, 192, 192)) * 1e-8,
        'bc1': np.random.random((64, 192, 192)) * 1e-7,
        'bc2': np.random.random((64, 192, 192)) * 1e-7,
        'oc1': np.random.random((64, 192, 192)) * 1e-7,
        'oc2': np.random.random((64, 192, 192)) * 1e-7,
        # ... etc for other species
    }
    
    # Get metadata for aerosol fields
    metadata = create_aerosol_field_metadata()
    
    try:
        # Example 1: Add fields to a single tile file
        logger.info("Adding aerosol fields to single tile file...")
        output_file = add_3d_fields_to_fv3_tile(
            tile_file,
            aerosol_data,
            metadata,
            backup_original=True
        )
        logger.info(f"Single tile processed: {output_file}")
        
        # Example 2: Add fields to all tile files in a directory
        logger.info("Adding aerosol fields to all tile files...")
        output_files = add_aerosol_fields_to_all_tiles(
            tile_directory,
            aerosol_data,
            metadata,
            backup_original=True
        )
        logger.info(f"Processed {len(output_files)} tile files")
        
    except Exception as e:
        logger.error(f"Example failed: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    example_usage()
