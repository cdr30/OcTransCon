"""
Module containing routines to save data to file.

"""

import numpy as np
import cPickle as pickle
import iris

def save_temporary_file(config, kout, nbasin):
    """ Save basin output as a temporary pickle file """
    
    tmpdir  = config.get('output', 'tmpdir')
    tmpname = config.get('output', 'name')
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
#     fwd_zsum = save_template.
#     
#     
#     fwd_ohc, fwd_flx, fwd_trans
#     bwd_ohc, bwd_flx, bwd_trans
#     
#     ohc_ob, flx_ob
#     ohc_ob_err
#     flx_ob_err
#     
#     

def save_as_netcdf(config, basins, save_template):
    """ Combine data from temporary files and save as netcdf files """
    
    savevars = ['flx_kf', 'flx_kb', 'tran_kf', 'flx_ob',
                 'zsum_kb_err', 'flx_ob_err', 'tran_kb', 'tran_kb_err',
                 'flx_kb_err', 'zsum_kf', 'zsum_kf_err', 'zsum_kb',
                 'flx_kf_err', 'zsum_ob', 'zsum_ob_err', 'tran_kf_err']
    
    nbasins = np.unique(basins.data.filled())
    
    
    for savevar in savevars:
        varname = config.get('output', savevar.split('_')[0])
        vartype = '_'.join(savevar.split('_')[1:])
        print savevar, varname + '_' + vartype
        
    
    zsum_cube = save_template.copy()
    zsum_cube.rename(config.get('output', 'zsum'))
    zsum_cube.units = config.get('output', 'zsumunits')
    
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
            print zsum_cube.data[:, bind[0], bind[1]].shape
            print kout['zsum_kb'][:, None].shape
            zsum_cube.data[:, bind[0], bind[1]] = kout['zsum_kb'][:, None]
        
    fout = '/Users/chris_roberts/test.nc'
    print 'SAVING: ' + fout
    iris.fileformats.netcdf.save(zsum_cube, fout)
        
        
        
        
        
        
        
        
        
        