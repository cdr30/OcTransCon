"""
Module containing routines to save data to file.

"""

import iris

def save_as_netcdf(config, out_cubes):
    """ Write output cubes to netcdf """
    
    for outvar in out_cubes.keys():
        outdir = config.get('output', 'dir')
        outname = config.get('output', 'name')
        out_cube = out_cubes[outvar]
        
        fout = '%s%s_%s.nc' % (outdir, outname, outvar)
        print 'SAVING: ' + fout
        iris.fileformats.netcdf.save(out_cube, fout) 
    