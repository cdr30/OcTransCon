"""
Module containing routines to save data to file.

"""

import numpy as np
import cPickle as pickle
import iris

def save_temporary_file(config, kout, nbasin):
    """ Save basin output as a temporary pickle file """
    
    tmpdir  = config.get('output', 'tmpdir')
    tmpname = config.get('output', 'label')
    tmpf = '%s%s_basin%i.tmp' % (tmpdir, tmpname, nbasin)
    
    print 'SAVING TEMPORARY FILE: ' + tmpf

    pickle.dump(kout, open(tmpf,'wb'))


# def create_output_cube(cube, name, units):
#     """ Update name and units within output cube. """
#     cube.rename(name)
#     cube.units = units
#   
#     return cube
#     
# 
# def create_output_cubes(config, save_template):
#     """ Create <iris.cube.Cube> objects to store output """
#     
#     
#     
#     
#     fwd_ohc = save_template.
#     
#     
#     fwd_ohc, fwd_flx, fwd_ohts
#     bwd_ohc, bwd_flx, bwd_ohts
#     
#     ohc_ob, flx_ob
#     ohc_ob_err
#     flx_ob_err
#     
#     

def save_as_netcdf(config, basins, save_template):
    """ Combine data from temporary files and save as netcdf files """
    
    savevars = [('ohc_ob', 'J'),
                ('ohc_ob_err', 'J'),
                ('ohc_kfwd', 'J'),
                ('ohc_kfwd_err', 'J'),
                ('ohc_ksmooth', 'J'),
                ('ohc_ksmooth_err', 'J'),
                ('flx_ob', 'W'),
                ('flx_ob_err', 'W'),
                ('flx_kfwd', 'W'),
                ('flx_kfwd_err', 'W'),
                ('flx_ksmooth', 'W'),                
                ('flx_ksmooth_err', 'W'),
                ('oht_kfwd', 'W'),
                ('oht_kfwd_err', 'W'),
                ('oht_ksmooth', 'W'),
                ('oht_ksmooth_err' 'W')]

# FINISH CHANGING VAR NAMES AND NAME OF REPOSITORY AND CHECK THAT CODE STILL RUNS FOR NATL.    
#   
    
    nbasins = np.unique(basins.data.filled())
    
    
    for savevar, units in savevars:
        print savevar, units
        
    
    ohc_cube = save_template.copy()
    ohc_cube.rename(config.get('output', 'ohc'))
    ohc_cube.units = config.get('output', 'ohcunits')
    
    for nbasin in nbasins:
        if nbasin != basins.data.fill_value:
            tmpdir  = config.get('output', 'tmpdir')
            tmpname = config.get('output', 'name')
            tmpf = '%s%s_basin%i.tmp' % (tmpdir, tmpname, nbasin)
            
            
            
            kout = pickle.load(open(tmpf, 'rb'))
            print kout.keys()
            print 'Loading temporary file  %i' % nbasin
            bind = np.where(basins.data.filled() == nbasin)
            
            #for j, i in zip(bind[0], bind[1]):
            print ohc_cube.data[:, bind[0], bind[1]].shape
            print kout['ohc_ksmooth'][:, None].shape
            ohc_cube.data[:, bind[0], bind[1]] = kout['ohc_ksmooth'][:, None]
        
    fout = '/Users/chris_roberts/test.nc'
    print 'SAVING: ' + fout
    iris.fileformats.netcdf.save(ohc_cube, fout)
        
        
        
        
        
        
        
        
        
        
