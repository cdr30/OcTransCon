"""
Module containing main routines to execute OceanDiv

"""

import parse_args
import namelist
import load
import process


def main():
    """
    Parse command line arguments and options and run OceanDiv.
    
    """
    args = parse_args.get_args()
    config = namelist.get_namelist(args)    
    
    # Load data
    ohcs = load.load_data(config, dtype='ohc')
    flxs = load.load_data(config, dtype='flx')
    basins = load.load_geodata(config, geotype='basins')
    areas = load.load_geodata(config, geotype='areas')

    # Process data 
    process.process_by_basin(config, ohcs, flxs, basins, areas)
    
        
        
        
    
