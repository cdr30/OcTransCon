"""
Module containing routines to pre-process data and invoke Kalman filter.

"""

import numpy as np
import scipy.signal
import kalman
import tools


class ShapeError(Exception):
    pass


def return_dates(cube):
    """ Return array of <datetime.datetime> objects. """
    
    tcoord = cube.coord('time')
    dts = tcoord.units.num2date(tcoord.points)

    return np.array(dts)
    

def try_append_mask(data, masks):
    """ Attempt to append data.mask to list of masks """
    
    try:
        masks.append(data.mask)
    except AttributeError:
        pass
    
    return masks


def update_mask(cube, mask):
    """ Return <iris.cube.Cube> with updated mask """
    
    if len(cube.shape) == 2:
        cube.data = np.ma.MaskedArray(cube.data, mask=mask)
    elif len(cube.shape) == 3:
        mask = np.array(np.ones_like(cube.data) * mask[np.newaxis, :, :], dtype=np.bool)
        cube.data = np.ma.MaskedArray(cube.data, mask=mask)
    else:
        raise ShapeError('Cannot mask data that is not 2D or 3D.')
        
    return cube
    
    
def unify_masks(ohcs, flxs, basins, areas):
    """ Unify masked areas in all data sets. """
    
    masks = []
    masks = try_append_mask(basins.data, masks)
    masks = try_append_mask(areas.data, masks)

    for ohc in ohcs:
        masks = try_append_mask(ohc[0].data, masks)
        
    for flx in flxs:
        masks = try_append_mask(flx[0].data, masks)

    new_mask = np.zeros_like(masks[0],dtype=np.int)
    for mask in masks:
        new_mask += np.array(mask, dtype=np.int)

    new_mask = new_mask > 0
    basins = update_mask(basins, new_mask)
    areas = update_mask(areas, new_mask)
    flxs = [update_mask(flx, new_mask) for flx in flxs]
    ohcs = [update_mask(ohc, new_mask) for ohc in ohcs]
    
    return ohcs, flxs, basins, areas
    
    
def low_pass(input_signal, cutoff):
    """ Low-pass data with Butterworth filter. """

    cutoff_freq = 1. / cutoff
    Wn = 2. * cutoff_freq 
    b, a = scipy.signal.butter(2, Wn, 'low')
    output_signal = scipy.signal.filtfilt(b, a, input_signal)
    
    return output_signal[cutoff/2:-cutoff/2]


def calc_basin_avg(cube, basins, areas, nbasin, total=False):
    """ Calculate area-weighted average for the specified basin. """
    
    bind = np.where(basins.data.filled() == nbasin)
    bareas = areas.data[bind[0], bind[1]]
    bdat = cube.data[:,bind[0], bind[1]]
    
    if total:
        bwts = bareas
    else:
        bwts = bareas/bareas.sum()
        
    bavg = (bwts[np.newaxis,:] * bdat).sum(axis=1)

    return bavg


def calc_ensemble_mean(datalist):
    """ Return ensemble mean of data in datalist """
    
    mn = np.zeros_like(datalist[0])
    
    for dat in datalist:
        mn += dat
        
    return mn/len(datalist)

 
def calc_ermsd(datalist, err_tol=1e-6):
    """ Return ensemble root-mean-square deviation for data in datalist"""
    
    mn = calc_ensemble_mean(datalist)
    
    ss = 0
    for dat in datalist:
        ss += np.sum((dat - mn) * (dat - mn))
    
    nk = len(datalist) * len(datalist[0])
    
    ermsd = np.sqrt(ss/nk)
    
    if ermsd < err_tol:
        raise ValueError('RMSD between obs must be > %f for uncertainty estimates' % err_tol)
    
    return ermsd


def process_basin(config, ohcs, flxs, basins, areas, nbasin):
    """ Pre-process data and invoke Kalman filter for the specified basin """
    
    # Extract dates
    dates = return_dates(ohcs[0])
    
    #Calc basin averages/totals
    use_total = config.getboolean('areas', 'calc_basin_totals')
    ohcs_bavg = [calc_basin_avg(cube, basins, areas, nbasin, total=use_total) for cube in ohcs]
    flxs_bavg = [calc_basin_avg(cube, basins, areas, nbasin, total=use_total) for cube in flxs]

    # Low pass data
    if config.getboolean('kfilter', 'smooth'):
        cutoff = config.getint('kfilter', 'cutoff')
        dates = dates[cutoff/2:-cutoff/2]
        ohcs_bavg = [low_pass(dat, cutoff) for dat in ohcs_bavg]
        flxs_bavg = [low_pass(dat, cutoff) for dat in flxs_bavg]
    
    #Reference ohc data to zero at t0.
    ohcs_bavg = [dat - dat[0] for dat in ohcs_bavg]
    
    # Calculate ensemble means
    flx_ob = calc_ensemble_mean(flxs_bavg)
    ohc_ob = calc_ensemble_mean(ohcs_bavg)
    
    # Calculate obs errors
    flx_ob_err = calc_ermsd(flxs_bavg)
    ohc_ob_err = calc_ermsd(ohcs_bavg)
    
    # Apply Kalman smoother
    kout = kalman.apply_ksmooth(config, flx_ob, ohc_ob, flx_ob_err, ohc_ob_err)
    
    # Append data to output dictionary
    kout['dates'] = dates
    kout['flx_ob'] = flx_ob
    kout['ohc_ob'] = ohc_ob
    kout['flx_ob_err'] = flx_ob_err
    kout['ohc_ob_err'] = ohc_ob_err
    
    return kout


def create_output_cubes(config, output_template):
    """ Create cubes for data output """
    
    outvars = [('ohc_ob', 'J', 3),
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
    
    out_cubes = {}
    
    for outvar, outunit, outdim in outvars:
        if outdim == 2:
            out_cube = output_template.copy()[0]
        elif outdim == 3:
            out_cube = output_template.copy()
    
        if not config.getboolean('areas', 'calc_basin_totals'):
            outunit += '/m2'
    
        out_cube.rename(outvar)
        out_cube.units = outunit
        out_cubes[outvar] = out_cube
    
    return out_cubes


def update_output_cubes(out_cubes, kout, basins, nbasin):
    """ 
    Update cubes with Kalman filter output for the specified basin in cubes.
    
    """

    for key in out_cubes.keys():
        cube = out_cubes[key]  
        bind = np.where(basins.data.filled() == nbasin)
        
        if len(cube.shape) == 3:
            cube.data[:, bind[0], bind[1]] = kout[key][:, None]
        elif len(cube.shape) == 2:
            cube.data[bind[0], bind[1]] = kout[key]
            
        out_cubes[key] = cube

    return out_cubes


def create_output_template(config, cube):
    """ 
    Create output cube and, if necessary, trim time axis 
    to account for low pass filtering. 
    
    """
    if config.getboolean('kfilter', 'smooth'):
        cutoff = config.getint('kfilter', 'cutoff')
        output_template = cube[cutoff/2:-cutoff/2].copy()
    else:
        output_template = cube.copy()
        
    return output_template


def process_by_basin(config, ohcs, flxs, basins, areas):
    """  For each basin, process data and invoke Kalman filter """
        
    ### Create cubes for data output using ohc as template
    output_template = create_output_template(config, ohcs[0])
    out_cubes = create_output_cubes(config, output_template)

    ### Process basins
    nbasins = np.unique(basins.data.filled())    
    for ntask, nbasin in enumerate(nbasins):
        if config.getboolean('output', 'print_stdout'): 
            tools.print_progress('Processing basins', len(nbasins), ntask+1, nbar=20)   
        if nbasin != basins.data.fill_value:
            kout = process_basin(config, ohcs, flxs, basins, areas, nbasin)
            out_cubes = update_output_cubes(out_cubes, kout, basins, nbasin)
            
    return out_cubes


