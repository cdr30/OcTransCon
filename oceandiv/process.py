"""
Module containing routines to pre-process data and invoke Kalman filter.

"""

import numpy as np
import scipy.signal
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

 
def calc_ermsd(datalist):
    """ Return ensemble root-mean-square deviation for data in datalist"""
    
    mn = calc_ensemble_mean(datalist)
    
    ss = 0
    for dat in datalist:
        ss += np.sum((dat - mn) * (dat - mn))
    
    nk = len(datalist) * len(datalist[0])
    
    return np.sqrt(ss/nk)


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
    
    # Write data to temporary file
    print 'Processing basin %i' % nbasin
    kout['dates'] = dates
    kout['flx_ob'] = flx_ob
    kout['ohc_ob'] = ohc_ob
    kout['flx_ob_err'] = flx_ob_err
    kout['ohc_ob_err'] = ohc_ob_err
    
    return kout
    
    #save.save_temporary_file(config, kout, nbasin)
    
    # Plot stuff - debugging
    
#     plt.plot(dates, ohc_ob, '0.5', linewidth=3)
#     plt.plot(dates, kout['ohc_kfwd'], 'k')
#     plt.plot(dates, kout['ohc_kfwd'] + kout['ohc_kfwd_err'], 'k:')
#     plt.plot(dates, kout['ohc_kfwd'] - kout['ohc_kfwd_err'], 'k:')
#      
#     plt.plot(dates, kout['ohc_ksmooth'], 'r')
#     plt.plot(dates, kout['ohc_ksmooth'] + kout['ohc_ksmooth_err'], 'r:')
#     plt.plot(dates, kout['ohc_ksmooth'] - kout['ohc_ksmooth_err'], 'r:')
#     plt.show()
#      
#      
#     plt.plot(dates, flx_ob, '0.5', linewidth=3)
#     plt.plot(dates, kout['flx_kfwd'], 'k')
#     plt.plot(dates, kout['flx_kfwd'] + kout['flx_kfwd_err'], 'k:')
#     plt.plot(dates, kout['flx_kfwd'] - kout['flx_kfwd_err'], 'k:')
#     
#     plt.plot(dates, kout['flx_ksmooth'], 'r')
#     plt.plot(dates, kout['flx_ksmooth'] + kout['flx_ksmooth_err'], 'r:')
#     plt.plot(dates, kout['flx_ksmooth'] - kout['flx_ksmooth_err'], 'r:')
#     plt.show()
#     
#     plt.plot(dates, kout['oht_kfwd'], 'k')
#     plt.plot(dates, kout['oht_kfwd'] + kout['oht_kfwd_err'], 'k:')
#     plt.plot(dates, kout['oht_kfwd'] - kout['oht_kfwd_err'], 'k:')
#     
#     plt.plot(dates, kout['oht_ksmooth'], 'r')
#     plt.plot(dates, kout['oht_ksmooth'] + kout['oht_ksmooth_err'], 'r:')
#     plt.plot(dates, kout['oht_ksmooth'] - kout['oht_ksmooth_err'], 'r:')
#     
#     plt.show()



def process_by_basin(config, ohcs, flxs, basins, areas):
    """  For each basin, process data and invoke Kalman filter """

    ### Create cubes for data output using ohc as template

    
    save_cubes = save.create_save_cubes(save_template)
    
    
    ### 
    nbasins = np.unique(basins.data.filled())
    
    for nbasin in nbasins:
        if nbasin != basins.data.fill_value:
            kout = process_basin(config, ohcs, flxs, basins, areas, nbasin)
    
            save_cubes = save.add_
            

    
    #save.save_as_netcdf(config, basins, save_template)        
                
    
    
    

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
#          Pass observed flx and ohc to Kalman filter.
#         
#          Apply Kalman filter and extract values
#         
#          Insert data into cubes
#         pass
#     
#      Save cubes.
        
