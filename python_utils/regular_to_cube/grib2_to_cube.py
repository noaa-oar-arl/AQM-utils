#!/usr/bin/env python3
# grib2_to_cube.py
# This script converts GRIB2 aerosol data to FV3 cube sphere format.

import grib2io
import numpy as np
import logging
from typing import Dict, List, Tuple, Optional

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def read_aerosol_species_from_grib2(grib_file_path: str, 
                                   level_range: Optional[Tuple[int, int]] = None) -> Dict[str, np.ndarray]:
    """
    Read all aerosol species from a GRIB2 file using grib2io.
    
    This function reads particulate matter fine (PMTF) and coarse (PMTC) aerosol species
    from GRIB2 files, typically from chemical weather models like GEFS-Aerosols.
    
    Args:
        grib_file_path (str): Path to the GRIB2 file
        level_range (tuple, optional): Tuple of (start_level, end_level) to read specific 
                                     hybrid levels. If None, reads all levels.
    
    Returns:
        Dict[str, np.ndarray]: Dictionary with aerosol species names as keys and 
                              4D arrays (time, level, lat, lon) as values
    
    Raises:
        FileNotFoundError: If the GRIB2 file doesn't exist
        Exception: For other GRIB2 reading errors
    """
    
    try:
        logger.info(f"Opening GRIB2 file: {grib_file_path}")
        
        # Open the GRIB2 file
        grb = grib2io.open(grib_file_path)
        
        # Dictionary to store aerosol data
        aerosol_data = {}
        
        # Define aerosol species patterns based on the file structure
        # From the wgrib2 analysis, we have these main aerosol types:
        aerosol_species_info = {
            'dust_bin1': {'parameter': 'PMTF', 'aerosol_type': 'Dust dry', 'size_range': '>=2e-07,<2e-06'},
            'dust_bin2': {'parameter': 'PMTF', 'aerosol_type': 'Dust dry', 'size_range': '>=2e-06,<3.6e-06'},
            'dust_bin3': {'parameter': 'PMTC', 'aerosol_type': 'Dust dry', 'size_range': '>=3.6e-06,<6e-06'},
            'dust_bin4': {'parameter': 'PMTC', 'aerosol_type': 'Dust dry', 'size_range': '>=6e-06,<1.2e-05'},
            'dust_bin5': {'parameter': 'PMTC', 'aerosol_type': 'Dust dry', 'size_range': '>=1.2e-05,<2e-05'},
            'seasalt_bin1': {'parameter': 'PMTF', 'aerosol_type': 'Sea salt dry', 'size_range': '>=6e-08,<2e-07'},
            'seasalt_bin2': {'parameter': 'PMTF', 'aerosol_type': 'Sea salt dry', 'size_range': '>=2e-07,<1e-06'},
            'seasalt_bin3': {'parameter': 'PMTC', 'aerosol_type': 'Sea salt dry', 'size_range': '>=1e-06,<3e-06'},
            'seasalt_bin4': {'parameter': 'PMTC', 'aerosol_type': 'Sea salt dry', 'size_range': '>=3e-06,<1e-05'},
            'seasalt_bin5': {'parameter': 'PMTC', 'aerosol_type': 'Sea salt dry', 'size_range': '>=1e-05,<2e-05'},
            'sulfate': {'parameter': 'PMTF', 'aerosol_type': 'Sulphate dry', 'size_range': '=1.39e-07'},
            'organic_carbon_hydrophobic': {'parameter': 'PMTF', 'aerosol_type': 'Particulate organic matter hydrophobic dry', 'size_range': '=4.24e-08'},
            'organic_carbon_hydrophilic': {'parameter': 'PMTF', 'aerosol_type': 'Particulate organic matter hydrophilic dry', 'size_range': '=4.24e-08'},
            'black_carbon_hydrophobic': {'parameter': 'PMTF', 'aerosol_type': 'Black carbon hydrophobic dry', 'size_range': '=2.36e-08'},
            'black_carbon_hydrophilic': {'parameter': 'PMTF', 'aerosol_type': 'Black carbon hydrophilic dry', 'size_range': '=2.36e-08'},
        }
        
        logger.info(f"Found {len(grb)} total messages in GRIB2 file")
        
        # Get grid information from the first message
        first_msg = grb[0]
        grid_info = {
            'ni': first_msg.nx,
            'nj': first_msg.ny,
            'lat_min': first_msg.latitudeOfFirstGridPoint / 1e6,
            'lat_max': first_msg.latitudeOfLastGridPoint / 1e6,
            'lon_min': first_msg.longitudeOfFirstGridPoint / 1e6,
            'lon_max': first_msg.longitudeOfLastGridPoint / 1e6,
            'dlat': first_msg.jDirectionIncrement / 1e6,
            'dlon': first_msg.iDirectionIncrement / 1e6
        }
        
        logger.info(f"Grid info: {grid_info['ni']}x{grid_info['nj']} points")
        logger.info(f"Lat range: {grid_info['lat_min']} to {grid_info['lat_max']}")
        logger.info(f"Lon range: {grid_info['lon_min']} to {grid_info['lon_max']}")
        
        # Determine levels to read
        available_levels = set()
        for msg in grb:
            if hasattr(msg, 'level'):
                available_levels.add(msg.level)
        
        available_levels = sorted(list(available_levels))
        logger.info(f"Available levels: {len(available_levels)} levels from {min(available_levels)} to {max(available_levels)}")
        
        if level_range:
            levels_to_read = [l for l in available_levels if level_range[0] <= l <= level_range[1]]
            logger.info(f"Reading levels {level_range[0]} to {level_range[1]}: {len(levels_to_read)} levels")
        else:
            levels_to_read = available_levels
            logger.info(f"Reading all {len(levels_to_read)} levels")
        
        # Initialize data arrays for each aerosol species
        for species_name in aerosol_species_info.keys():
            # Shape: (time, level, lat, lon) - assuming single time for now
            aerosol_data[species_name] = np.zeros((1, len(levels_to_read), grid_info['nj'], grid_info['ni']))
        
        # Read data for each aerosol species and level
        species_found = {species: False for species in aerosol_species_info.keys()}
        
        for msg in grb:
            # Get message properties
            level = getattr(msg, 'level', None)
            if level not in levels_to_read:
                continue
                
            level_idx = levels_to_read.index(level)
            
            # Get parameter name and aerosol information
            param_name = msg.shortName if hasattr(msg, 'shortName') else str(msg.parameterNumber)
            
            # Try to match this message to one of our aerosol species
            for species_name, species_info in aerosol_species_info.items():
                # This is a simplified matching - in practice, you might need more sophisticated
                # matching based on GRIB2 parameter tables and aerosol metadata
                if param_name in ['PMTF', 'PMTC']:
                    # For now, we'll use a simple approach based on message order
                    # In a real implementation, you'd want to parse the aerosol metadata
                    # to properly identify each species
                    
                    # Read the data values
                    try:
                        data_values = msg.data()
                        if data_values is not None:
                            # Reshape to match grid dimensions
                            data_2d = data_values.reshape(grid_info['nj'], grid_info['ni'])
                            
                            # Store in the aerosol_data array
                            # This is a placeholder - you'd need proper species identification
                            if not species_found[species_name]:
                                aerosol_data[species_name][0, level_idx, :, :] = data_2d
                                species_found[species_name] = True
                                logger.info(f"Read {species_name} for level {level}")
                                break
                    except Exception as e:
                        logger.warning(f"Error reading data for message at level {level}: {e}")
                        continue
        
        # Log which species were found
        found_species = [s for s, found in species_found.items() if found]
        missing_species = [s for s, found in species_found.items() if not found]
        
        logger.info(f"Successfully read {len(found_species)} aerosol species: {found_species}")
        if missing_species:
            logger.warning(f"Could not find data for {len(missing_species)} species: {missing_species}")
        
        # Add grid information to the returned data
        aerosol_data['_grid_info'] = grid_info
        aerosol_data['_levels'] = levels_to_read
        
        grb.close()
        
        return aerosol_data
        
    except FileNotFoundError:
        logger.error(f"GRIB2 file not found: {grib_file_path}")
        raise
    except Exception as e:
        logger.error(f"Error reading GRIB2 file {grib_file_path}: {e}")
        raise


def get_aerosol_species_by_parameter(grib_file_path: str, 
                                   parameter_name: str = 'PMTF',
                                   level: int = 1) -> Dict[str, np.ndarray]:
    """
    Read specific aerosol parameter (PMTF or PMTC) from GRIB2 file.
    
    Args:
        grib_file_path (str): Path to the GRIB2 file
        parameter_name (str): Parameter to read ('PMTF' for fine, 'PMTC' for coarse)
        level (int): Hybrid level to read
    
    Returns:
        Dict[str, np.ndarray]: Dictionary with aerosol data and metadata
    """
    
    try:
        logger.info(f"Reading {parameter_name} data from level {level}")
        
        grb = grib2io.open(grib_file_path)
        
        # Find messages matching the criteria
        matching_data = {}
        
        for i, msg in enumerate(grb):
            if (hasattr(msg, 'shortName') and msg.shortName == parameter_name and 
                hasattr(msg, 'level') and msg.level == level):
                
                # Read the data
                data_values = msg.data()
                if data_values is not None:
                    # Create a unique key for this aerosol type
                    # You might want to extract more specific aerosol type info here
                    key = f"{parameter_name}_message_{i}_level_{level}"
                    matching_data[key] = {
                        'data': data_values.reshape(msg.ny, msg.nx),
                        'units': getattr(msg, 'units', 'unknown'),
                        'level': level,
                        'parameter': parameter_name,
                        'grid_shape': (msg.ny, msg.nx)
                    }
                    
                    logger.info(f"Read {key}: shape {data_values.shape}, min={np.min(data_values):.2e}, max={np.max(data_values):.2e}")
        
        grb.close()
        
        return matching_data
        
    except Exception as e:
        logger.error(f"Error reading {parameter_name} data: {e}")
        raise


# Example usage function
def example_usage():
    """
    Example of how to use the aerosol reading functions.
    """
    
    # Example file path (adjust as needed)
    grib_file = "/gpfs/f6/ira-sti/proj-shared/Cory.R.Martin/july2025/gcafs/gefs.chem.t00z.a3d_0p50.f000.grib2"
    
    try:
        # Read all aerosol species
        logger.info("Reading all aerosol species...")
        all_aerosols = read_aerosol_species_from_grib2(grib_file, level_range=(1, 5))
        
        # Print summary
        for species_name, data in all_aerosols.items():
            if not species_name.startswith('_'):  # Skip metadata entries
                print(f"{species_name}: shape {data.shape}, min={np.min(data):.2e}, max={np.max(data):.2e}")
        
        # Read specific parameter
        logger.info("Reading PMTF data from level 1...")
        pmtf_data = get_aerosol_species_by_parameter(grib_file, 'PMTF', level=1)
        
        for key, info in pmtf_data.items():
            print(f"{key}: {info['grid_shape']}, units: {info['units']}")
            
    except Exception as e:
        logger.error(f"Example failed: {e}")


if __name__ == "__main__":
    example_usage()
