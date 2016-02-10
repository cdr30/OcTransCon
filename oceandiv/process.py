"""
Module containing routines to pre-process data and invoke Kalman filter.

"""

import numpy as np
import scipy.signal
import math

import kalman
import save

import matplotlib.pyplot as plt


class ShapeError(Exception):
    pass


def return_dates(cube):
    """ Return array of <datetime.datetime> objects. """
    
    tcoord = cube.coord('time')
    dts = tcoord.units.num2date(tcoord.points)

    return np.array(dts)
    

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
    
    # Extract dates
    dates = return_dates(zsums[0])
    
    #Calc basin averages/totals
    use_total = config.getboolean('areas', 'calc_basin_totals')
    zsums_bavg = [calc_basin_avg(cube, basins, areas, nbasin, total=use_total) for cube in zsums]
    flxs_bavg = [calc_basin_avg(cube, basins, areas, nbasin, total=use_total) for cube in flxs]

    # Low pass data
    if config.getboolean('kfilter', 'smooth'):
        cutoff = config.getint('kfilter', 'cutoff')
        dates = dates[cutoff/2:-cutoff/2]
        zsums_bavg = [low_pass(dat, cutoff) for dat in zsums_bavg]
        flxs_bavg = [low_pass(dat, cutoff) for dat in flxs_bavg]
    
    #Reference zsum data to zero at t0.
    zsums_bavg = [dat - dat[0] for dat in zsums_bavg]
    
    # Calculate ensemble means
    flx_ob = calc_ensemble_mean(flxs_bavg)
    zsum_ob = calc_ensemble_mean(zsums_bavg)
    
    # Calculate obs errors
    flx_ob_err = calc_ermsd(flxs_bavg)
    zsum_ob_err = calc_ermsd(zsums_bavg)
    
    # Apply Kalman smoother
    kout = kalman.apply_ksmooth(config, flx_ob, zsum_ob, flx_ob_err, zsum_ob_err)
    
    # Write data to temporary file
    print 'Processing basin %i' % nbasin
    kout['dates'] = dates
    kout['flx_ob'] = flx_ob
    kout['zsum_ob'] = zsum_ob
    kout['flx_ob_err'] = flx_ob_err
    kout['zsum_ob_err'] = zsum_ob_err
    
    save.save_temporary_file(config, kout, nbasin)
    
    # Plot stuff - debugging
    
#     plt.plot(dates, zsum_ob, '0.5', linewidth=3)
#     plt.plot(dates, kout['zsum_kf'], 'k')
#     plt.plot(dates, kout['zsum_kf'] + kout['zsum_kf_err'], 'k:')
#     plt.plot(dates, kout['zsum_kf'] - kout['zsum_kf_err'], 'k:')
#     
#     plt.plot(dates, kout['zsum_kb'], 'r')
#     plt.plot(dates, kout['zsum_kb'] + kout['zsum_kb_err'], 'r:')
#     plt.plot(dates, kout['zsum_kb'] - kout['zsum_kb_err'], 'r:')
#     plt.show()
#     
#     
#     plt.plot(dates, flx_ob, '0.5', linewidth=3)
#     plt.plot(dates, kout['flx_kf'], 'k')
#     plt.plot(dates, kout['flx_kf'] + kout['flx_kf_err'], 'k:')
#     plt.plot(dates, kout['flx_kf'] - kout['flx_kf_err'], 'k:')
#     
#     plt.plot(dates, kout['flx_kb'], 'r')
#     plt.plot(dates, kout['flx_kb'] + kout['flx_kb_err'], 'r:')
#     plt.plot(dates, kout['flx_kb'] - kout['flx_kb_err'], 'r:')
#     plt.show()
#     
#     plt.plot(dates, kout['tran_kf'], 'k')
#     plt.plot(dates, kout['tran_kf'] + kout['tran_kf_err'], 'k:')
#     plt.plot(dates, kout['tran_kf'] - kout['tran_kf_err'], 'k:')
#     
#     plt.plot(dates, kout['tran_kb'], 'r')
#     plt.plot(dates, kout['tran_kb'] + kout['tran_kb_err'], 'r:')
#     plt.plot(dates, kout['tran_kb'] - kout['tran_kb_err'], 'r:')
#     
#     plt.show()



def process_by_basin(config, zsums, flxs, basins, areas):
    """  For each basin, process data and invoke Kalman filter """
    
    zsums, flxs, basins, areas = unify_masks(zsums, flxs, basins, areas)
    nbasins = np.unique(basins.data.filled())
    
    for nbasin in nbasins:
        if nbasin != basins.data.fill_value:
            process_basin(config, zsums, flxs, basins, areas, nbasin)
            
    if config.getboolean('kfilter', 'smooth'):
        cutoff = config.getint('kfilter', 'cutoff')
        save_template = zsums[0].copy()[cutoff/2:-cutoff/2]
    else:
        save_template = zsums[0].copy()
    
    save.save_as_netcdf(config, basins, save_template)        
                
    
    
    

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
        