"""
Module containing main routines to process data 

"""


import numpy as np
import scipy.signal
import math



import matplotlib.pyplot as plt


class ShapeError(Exception):
    pass


def try_append_mask(data, masks):
    """ Attempt to append data.mask to masks """
    
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
    
    
def unify_masks(zsums, flxs, basins, areas):
    """ Unify masked areas in all data sets. """
    
    masks = []
    masks = try_append_mask(basins.data, masks)
    masks = try_append_mask(areas.data, masks)

    for zsum in zsums:
        masks = try_append_mask(zsum[0].data, masks)
        
    for flx in flxs:
        masks = try_append_mask(flx[0].data, masks)

    new_mask = np.zeros_like(masks[0],dtype=np.int)
    for mask in masks:
        new_mask += np.array(mask, dtype=np.int)

    new_mask = new_mask > 0
    basins = update_mask(basins, new_mask)
    areas = update_mask(areas, new_mask)
    flxs = [update_mask(flx, new_mask) for flx in flxs]
    zsums = [update_mask(zsum, new_mask) for zsum in zsums]
    
    return zsums, flxs, basins, areas
    
    
def low_pass(input_signal, cutoff):
    """ Low-pass data with Butterworth filter. """

    cutoff_freq = 1. / cutoff
    Wn = 2. * cutoff_freq 
    b, a = scipy.signal.butter(2, Wn, 'low')
    output_signal = scipy.signal.filtfilt(b, a, input_signal)
    
    return output_signal[cutoff/2:-cutoff/2]


def calc_basin_avg(cube, basins, areas, nbasin):
    """ Calculate area-weighted average for the specified basin. """
    
    bind = np.where(basins.data.filled() == nbasin)
    bareas = areas.data[bind[0], bind[1]]
    bdat = cube.data[:,bind[0], bind[1]]
    bwts = bareas/bareas.sum()
    bavg = (bwts[np.newaxis,:] * bdat).sum(axis=1)

    return bavg


def calc_ensemble_mean(datalist):
    """ Return ensemble mean of data in datalist """
    
    mn = np.zeros_like(datalist[0])
    
    for dat in datalist:
        mn += dat
        
    return mn/len(datalist)

 
def calc_ermsd(datalist):
    """ Return ensemble root-mean-square deviation for data in datalist"""
    
    mn = calc_ensemble_mean(datalist)
    
    ss = 0
    for dat in datalist:
        ss += np.sum((dat - mn) * (dat - mn))
    
    nk = len(datalist) * len(datalist[0])
    
    return np.sqrt(ss/nk)


def process_basin(config, zsums, flxs, basins, areas, nbasin):
    """ Pre-process data and invoke Kalman filter for the specified basin """
    
    cutoff = config.getint('kfilter', 'cutoff')
    
    #Calc basin averages
    zsums_bavg = [calc_basin_avg(cube, basins, areas, nbasin) for cube in zsums]
    flxs_bavg = [calc_basin_avg(cube, basins, areas, nbasin) for cube in flxs]
    
    # Low pass data
    zsums_lp = [low_pass(dat, cutoff) for dat in zsums_bavg]
    flxs_lp = [low_pass(dat, cutoff) for dat in flxs_bavg]
    
    #Reference OHC data to zero at t0.
    zsums_lp = [dat - dat[0] for dat in zsums_lp]
    
    # Calculate ensemble means
    flx_ob = calc_ensemble_mean(flxs_lp)
    zsum_ob = calc_ensemble_mean(zsums_lp)
    
    # Calculate obs errors
    flx_ob_err = calc_ermsd(flxs_lp)
    zsum_ob_err = calc_ermsd(zsums_lp)
    
    # Apply Kalman smoother
    kalman.ksmooth(flx_ob, zsum_ob, flx_ob_err, zsum_ob_err)
        
    plt.plot(zsums_lp[0], 'r')
    plt.plot(zsums_lp[1], 'r')
    plt.plot(zsum_ob, 'k')
    plt.plot(zsum_ob + zsum_ob_err * 1, 'k:')
    plt.plot(zsum_ob - zsum_ob_err * 1, 'k:'); plt.show()

    

def process_by_basin(config, zsums, flxs, basins, areas):
    """  For each basin, process data and invoke Kalman filter """
    
    nbasins = np.unique(basins.data.filled())
    
    for nbasin in nbasins:
        if nbasin != basins.data.fill_value:
            process_basin(config, zsums, flxs, basins, areas, nbasin)
            
            

#     
#      
#      Copy cubes to store ksmoothed OHC, Flx and OHT and associated errors.
# 
#      Find list of unique basins
#     nbasins = find_nbasins(config, basins)
#     
#     for nbasin in nbasins:
#         
#          Extract area-weighted basin data
# 
#          Smooth data
#         
#          Calculate obs errors
#         
#          Calculate ensemble  means
#         
#          Pass observed flx and zsum to Kalman filter.
#         
#          Apply Kalman filter and extract values
#         
#          Insert data into cubes
#         pass
#     
#      Save cubes.
        