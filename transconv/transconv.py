"""
Module containing main routines to execute TransConv.

"""

import parse_args
import namelist
import load
import process


def main():
    """
    Parse command line arguments and options and run TransConv.
    
    """
    args = parse_args.get_args()
    config = namelist.get_namelist(args)    
    
    # Load data
    zsums = load.load_data(config, dtype='zsum')
    flxs = load.load_data(config, dtype='flux')
    basins = load.load_geodata(config, geotype='basins')
    areas = load.load_geodata(config, geotype='areas')
    
    print flxs[0][0]
    print flxs[0][-1]
    
    
    # Process data 
    zsums, flxs, basins, areas = process.unify_masks(zsums, flxs, basins, areas)
    process.process_by_basin(config, zsums, flxs, basins, areas)
    
   
        
        
        
    