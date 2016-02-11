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


def load_temporary_file(config, nbasin):
    """ Load basin output from temporary pickle file """
    
    tmpdir  = config.get('output', 'tmpdir')
    tmpname = config.get('output', 'label')
    tmpf = '%s%s_basin%i.tmp' % (tmpdir, tmpname, nbasin)           
    kout = pickle.load(open(tmpf, 'rb'))
    
    return kout

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

def create_save_template(config, cube):
    """ Create template cube for data output """
    
    if config.getboolean('kfilter', 'smooth'):
        cutoff = config.getint('kfilter', 'cutoff')
        save_template = cube.copy()[cutoff/2:-cutoff/2]
    else:
        save_template = cube.copy()
        
    return cube



def create_save_cubes():
    """ Create cubes for data output """
    
    
        savevars = [('ohc_ob', 'J', 3),
                ('ohc_ob_err', 'J', 2),
                ('ohc_kfwd', 'J', 3),
                ('ohc_kfwd_err', 'J', 3),
                ('ohc_ksmooth', 'J', 3),
                ('ohc_ksmooth_err', 'J', 3),
                ('flx_ob', 'W', 3),
                ('flx_ob_err', 'W', 2),
                ('flx_kfwd', 'W', 3),
                ('flx_kfwd_err', 'W', 3),
                ('flx_ksmooth', 'W', 3),                
                ('flx_ksmooth_err', 'W', 3),
                ('oht_kfwd', 'W', 3),
                ('oht_kfwd_err', 'W', 3),
                ('oht_ksmooth', 'W', 3),
                ('oht_ksmooth_err', 'W', 3)]
    

def save_as_netcdf(config, basins, save_template):
    """ Combine data from temporary files and save as netcdf files """
    
    savevars = [('ohc_ob', 'J', 3),
                ('ohc_ob_err', 'J', 2),
                ('ohc_kfwd', 'J', 3),
                ('ohc_kfwd_err', 'J', 3),
                ('ohc_ksmooth', 'J', 3),
                ('ohc_ksmooth_err', 'J', 3),
                ('flx_ob', 'W', 3),
                ('flx_ob_err', 'W', 2),
                ('flx_kfwd', 'W', 3),
                ('flx_kfwd_err', 'W', 3),
                ('flx_ksmooth', 'W', 3),                
                ('flx_ksmooth_err', 'W', 3),
                ('oht_kfwd', 'W', 3),
                ('oht_kfwd_err', 'W', 3),
                ('oht_ksmooth', 'W', 3),
                ('oht_ksmooth_err', 'W', 3)]

# FINISH CHANGING VAR NAMES AND NAME OF REPOSITORY AND CHECK THAT CODE STILL RUNS FOR NATL.    
#   
    nbasins = np.unique(basins.data.filled())
    
    
    for savevar, saveunit, savedim in savevars:
        
        if savedim == 2:
            save_cube = save_template.copy()[0]
        elif savedim == 3:
            save_cube = save_template.copy()
    
        if not config.getboolean('areas', 'calc_basin_totals'):
            saveunit += '/m2'
    
        save_cube.rename(savevar)
        save_cube.units = saveunit
    
        for nbasin in nbasins:
            if nbasin != basins.data.fill_value:

                kout = load_temporary_file(config, nbasin)
                bind = np.where(basins.data.filled() == nbasin)
                
                if savedim == 3:
                    save_cube.data[:, bind[0], bind[1]] = kout[savevar][:, None]
                elif savedim == 2:
                    save_cube.data[bind[0], bind[1]] = kout[savevar]
                    
                #print 'Loading temporary file  %i' % nbasin
        
        outdir = config.get('output', 'outdir')
        label = config.get('output', 'label')
        fout = outdir + label + '_' + savevar + '.nc'
        print 'SAVING: ' + fout
        iris.fileformats.netcdf.save(save_cube, fout)
        
            
            #for j, i in zip(bind[0], bind[1]):
            #print ohc_cube.data[:, bind[0], bind[1]].shape
            #print kout['ohc_ksmooth'][:, None].shape
            
        
    #fout = '/Users/chris_roberts/test.nc'
    
    
        
        
        
        
        
        
        
        
        
        
