
import numpy as np
from scipy.interpolate import interp1d
from pylab import *

MS2KTS = 1.94
MISSING = np.nan

def _removedMasked(ary):
    if np.ma.is_masked(ary):
        ary_mask = ary.mask
        ary = ary[np.where(~ary_mask)]
    return ary

def _RH(temp, dwpt, pres):
    return 100*(np.exp((17.625*dwpt)/(243.04+dwpt))/np.exp((17.625*temp)/(243.04+temp))) 

def _groupLayers(cand_layers):
    """
    Function to group candidate layers into contiguous regions.
    """
    group_layers = np.where(np.diff(cand_layers) == 1, 1, 0)
    try:
        breaks = np.where(group_layers == 0)[0] + 1
    except IndexError:
        breaks = []
    breaks = [ 0 ] + list(breaks) + [ -1 ]

    cand_idxs = []
    for idx in xrange(len(breaks) - 1):
        cand_idxs.append(cand_layers[breaks[idx]:breaks[idx + 1]])
    
    return cand_idxs

def _splitProfile(prof, lb_idx, tol, _depth=0):
    """
    _splitProfile()
    Purpose:    Thin a profile in such a way that none of the intermediate 
                    values deviate from a linear interplation of the thinned
                    values by more than a specified tolerance.
    Arguments:  prof [type=np.ndarray]
                    The profile to thin
                lb_idx [type=int]
                    The position (in the full profile) of the first index in 
                    this segment of the profile.
                tol [type=float]
                    Tolerance on the linear interpolation
    Returns:    A list of indices representing those in the thinned profile
                    
    """
    interp = interp1d([0., len(prof) - 1.], [prof[0], prof[-1]])
    interp_prof = interp(np.arange(float(len(prof))))

    if np.abs(interp_prof - prof).max() < tol:
        ret_val = [ lb_idx ]
    else:
        split_point = np.argmax(np.abs(prof - interp_prof))
        if split_point > 1:
            first_split = _splitProfile(prof[:split_point], lb_idx, tol, _depth=_depth + 1)
        else:
            first_split = [ lb_idx ]

        if split_point < len(prof) - 2:
            second_split = _splitProfile(prof[split_point:], lb_idx + split_point, tol, _depth=_depth + 1)
        else:
            second_split = [ lb_idx + split_point ]

        ret_val = first_split + second_split

    if _depth == 0:
        ret_val += [ len(prof) - 1 ]
    return ret_val

def _mandatoryPresLevels(pres, tol=0.1):
    '''
        Searches for the surface, 1000, 925, 850, 700,
        500, 400, 300, 250, 200, 150, 100, 70, 50, 30, 20, and 10 mb.
    '''
    mandatory_pres = np.asarray([1000, 925, 850, 700, 500, 400, 300, 250, 200,\
                      150, 100, 70, 50, 30, 20, 10])
    mandatory_pres_idx = np.ones(mandatory_pres.shape, dtype=int)
    pres = np.round(pres,0)
    for i in xrange(len(mandatory_pres_idx)):
        ind = np.where(np.fabs(pres - mandatory_pres[i]) < tol)[0]
        if len(ind) == 0:
            ind = [-9999]
        mandatory_pres_idx[i] = ind[0]

    return mandatory_pres_idx

def _findTropopause(**snd):
    '''
        Finding the tropopause level using the WMO definition of a 
        tropopause as being the lowest level where the 2 km layer
        aloft has a lapse rate greater than 2 C/km.
    '''

    def mean(array):
        return (array[1:] + array[:-1]) * 0.5

    above_500 = np.where((snd['pres'] > 75) & (snd['pres'] < 550))[0]
    temp = snd['temp'][above_500]
    hght = snd['hght'][above_500]
    dT = np.diff(temp)
    dz = np.diff(hght)
    mean_z = mean(hght)

    gamma = dT / dz

    keep_gamma = np.isfinite(gamma) & (gamma > 0.002)
    if type(gamma) == np.ma.core.MaskedArray:
        keep_gamma &= (~gamma.mask)
    idx = np.where(keep_gamma)[0]

    for i in idx:
        hbot = mean_z[i]
        htop = hbot + 2000.
        htop_idx = np.argmin(np.abs(mean_z - htop))
        mean_lapse = (temp[htop_idx] - temp[i]) / (hght[htop_idx] - hght[i])
        if mean_lapse > 0.002:
            return above_500[i]
    return MISSING

def _unfoldWindDir(wdir):
    fold_idxs = np.where((wdir[:-1] > 350) & (wdir[1:] < 10) | (wdir[:-1] < 10) & (wdir[1:] > 350))[0]
    wdir_unfold = wdir.copy()

    for idx in fold_idxs:
        if wdir[idx] > wdir[idx + 1]:
            wdir_unfold[(idx + 1):] += 360
        else:
            wdir_unfold[(idx + 1):] -= 360

    return wdir_unfold

def _findMaxWind(wspd, min_wspd=(30 * MS2KTS), min_diff=(10 * MS2KTS), smooth_pts=11):
    '''
    _maxWind()
    Purpose:    Find the level of the maximum winds in the profile
    Arguments:  the snd dictionary.
    Returns:    The index from the sounding of the maximum wind speed.

    Assumes that the winds are in knots.
    '''
    
    # Smooth the wind profile
    wgts = np.ones((smooth_pts,)) / smooth_pts
    wspd_smooth = np.convolve(wspd, wgts, mode='same')

    # Find the indices of the maxima in the smoothed profile
    max_idxs = np.where(
        (wspd_smooth[1:-1] > wspd_smooth[:-2]) & 
        (wspd_smooth[1:-1] >= wspd_smooth[2:])
    )[0] + 1

    # Sort them by wind speed
    sorted_max_idxs = np.argsort(wspd[max_idxs])[::-1]
    max_idxs = max_idxs[sorted_max_idxs]

    # Cut off the wind speeds at 30 m/s, leaving only the max wind speed
    max_idxs = np.concatenate((
        max_idxs[:1],
        max_idxs[np.where(wspd[max_idxs[1:]] > min_wspd)]
    ))

    if len(max_idxs) > 1:
        # If we still have any other wind maxima left, cut off the ones that 
        #  are less than 10 m/s less than the previous
        ws_diff = np.diff(wspd[max_idxs])
        diff_cutoff = np.where(ws_diff < min_diff)[0][0] + 1
        max_idxs = max_idxs[:diff_cutoff]

    return max_idxs

def _windSigLevels(standard='RWS', **snd):
    wdir = _removedMasked(snd['wdir'])
    wspd = _removedMasked(snd['wspd'])

    wdir_unfold = _unfoldWindDir(wdir)

    wspd_tol = 0.5 * MS2KTS
    wdir_tol = 10

    max_wind_sl = _findMaxWind(wspd)
    wind_spd_sl = _splitProfile(wspd, 0, wspd_tol)
    wind_dir_sl = _splitProfile(wdir_unfold, 0, wdir_tol)

    wind_sl = np.unique(np.concatenate((max_wind_sl, wind_spd_sl, wind_dir_sl)))
    return np.sort(wind_sl)

def _findInversions(temp, dewp, pres, trop_idx):
    lr_layers = np.diff(temp)
    cand_layers = np.where(lr_layers > 0)[0]
    groupedLayers = _groupLayers(cand_layers)
    rh = _RH(temp, dewp, pres)

    inversion_layers = []
    for layer in groupedLayers:
        lbot = layer[0]
        ltop = layer[-1]
        if pres[lbot] < 300 or lbot > trop_idx:
            break
        elif (temp[ltop] - temp[lbot] >= 2.5) or \
             (np.abs(rh[ltop] - rh[lbot]) > 20) or \
             (pres[lbot] - pres[ltop] >= 20):
            inversion_layers.append(lbot)
            inversion_layers.append(ltop)
    return np.asarray(inversion_layers, dtype=int)

def _findIsothermals(temp, pres, tol=0.5):
    """
    _findIsothermals()
    Purpose:    Find the significant levels associated with the tops and 
                    bottoms of isothermal layers.
    Arguments:  temp [type=np.ndarray]
                    Temperature array in degrees C
                pres [type=np.ndarray]
                    Pressure array in hPa
    Keywords:   tol [type=float]
                    Temperature tolerance on "isothermal."  Default is 0.5 C.
    Returns:    A list of indices of the tops and bottoms of isothermal layers.
    """
    def pareLayer(cand_idxs):
        """
        Function to take a layer and pare it down to only the truly isothermal
        portion (not just the small lapse rate).  The process may split the
        layer, so if it does, pare all the sublayers and return the result.
        """
        temp_layer = temp[cand_idxs]

        while (len(temp_layer) > 1 and 
            np.max(temp_layer - temp_layer.mean()) > tol):

            worst = np.argmax(temp_layer - temp_layer.mean())
            cand_idxs = np.delete(cand_idxs, worst)

            cand_idxs = _groupLayers(cand_idxs)
            if len(cand_idxs) == 1:
                cand_idxs = cand_idxs[0]
            else:
                split_idxs = []
                for cidx in cand_idxs:
                    cidx = pareLayer(cidx)
                    for sidx in cidx:
                        if len(sidx) > 0:
                            split_idxs.append(sidx)
                return split_idxs

            temp_layer = temp[cand_idxs]

        ret_val = []
        pres_cand = pres[cand_idxs]
        if len(pres_cand) > 0 and pres_cand.max() - pres_cand.min() >= 20:
            ret_val = [ cand_idxs ]
        return ret_val

    isotherm_idxs = []

    lr_layers = np.abs(np.diff(temp))
    cand_layers = np.where(lr_layers < 2 * tol)[0]
    cand_idxs = _groupLayers(cand_layers)

    for cand_idx in cand_idxs:
        cand_idx = pareLayer(cand_idx)

        for cidx in cand_idx:
            isotherm_idxs.extend([cidx[0], cidx[-1]])

    return isotherm_idxs

def _findAddSigLevels(temp, pres, trop_idx, tol=[1.0, 2.0], max_trop=300):

    cutoff_idx = min(np.argmin(np.abs(pres - max_trop)), trop_idx)

    idxs_below = _splitProfile(temp[:cutoff_idx], 0, tol[0])[:-1]
    idxs_above = _splitProfile(temp[cutoff_idx:], cutoff_idx, tol[1])
    additional_idxs = np.sort(np.unique(np.concatenate((idxs_below, idxs_above))))

    return additional_idxs

def _findSigRHLevels(temp, dwpt, pres, tol=[15]):
    """
    _findSigRHLevels()
    Purpose:    Convert dwpt to relative humidity and find the
                relative humidity significant levels.
    Keywords:   temp [type=np.ndarray]
                    Temperature array in degrees C
                dewp [type=np.ndarray]
                    Dewpoint array in degrees C
                pres [type=np.ndarray]
                    Pressure array in hPa
                tol [type=float]
                    Tolerance on the linear interpolation                        
    Returns:    A list of indices represending the indices of the significant RH
                    levels.
    """
    rh = _RH(temp, dwpt, pres)
    idxs_rhSig = _splitProfile(rh, 0, tol[0])
    additional_idxs = np.sort(idxs_rhSig)

    return additional_idxs

def _thermSigLevels(standard='RWS', **snd):
    """
    _thermSigLevels()
    Purpose:    Find the temperature significant levels
    Keywords:   [ same as for the soundingFilter() function ]
    Returns:    A list of indices represending the indices of the significant
                    levels.
    """

    if standard == 'RWS':
        temp_tol = [0.5, 1.0]
        rh_tol = [5.]
        trop_idx = len(snd['pres']) 
        max_trop = 100.

    elif standard == 'WMO':
        temp_tol = [1.0, 2.0]
        rh_tol = [10.]
        trop_idx = _findTropopause(**snd)
        max_trop = 300.

    inv_sl = _findInversions(snd['temp'], snd['dewp'], snd['pres'], trop_idx)
    iso_sl = _findIsothermals(snd['temp'], snd['pres'])
    add_sl = _findAddSigLevels(snd['temp'], snd['pres'], trop_idx, tol=temp_tol, max_trop=max_trop)
    rh_sl = _findSigRHLevels(snd['temp'], snd['dewp'], snd['pres'], tol=rh_tol)

    therm_sl = np.unique(np.concatenate((inv_sl, iso_sl, add_sl, rh_sl)))
    return np.sort(therm_sl)

def soundingFilter(missing=MISSING, standard='RWS', **snd):
    """
    soundingFilter()
    Purpose:    Filter high-resolution sounding observations by finding the 
                    standard and significant levels.
    Keywords:   temp [type=np.ndarray]
                    Temperature array in degrees C
                dewp [type=np.ndarray]
                    Dewpoint array in degrees C
                pres [type=np.ndarray]
                    Pressure array in hPa
                wspd [type=np.ndarray]
                    Wind speed array in kts
                wdir [type=np.ndarray]
                    Wind direction array in degrees
    Returns:    A dictionary contaning filtered sounding observations.  Should
                    levels not line up between temperature and wind significant
                    levels, missing observations will be filled with np.nan.
    """

    pres_ml = _mandatoryPresLevels(snd['pres'])
    trop_idx = _findTropopause(**snd)
    therm_sl = _thermSigLevels(standard=standard, **snd)
    wind_sl = _windSigLevels(standard=standard, **snd)

    pres_ml = np.asarray(pres_ml[pres_ml >= 0], dtype=int)

    therm_sl = np.concatenate((therm_sl, pres_ml))
    wind_sl = np.concatenate((wind_sl, pres_ml))

    if np.isfinite(trop_idx):
        therm_sl = np.concatenate((therm_sl, [ trop_idx ]))
        wind_sl = np.concatenate((wind_sl, [ trop_idx ]))

    if therm_sl.dtype == np.float64:
        therm_sl = np.asarray(therm_sl[np.isfinite(therm_sl)], dtype=np.int64)

    if wind_sl.dtype == np.float64:
        wind_sl = np.asarray(wind_sl[np.isfinite(wind_sl)], dtype=np.int64)

    pres_sl = np.union1d(therm_sl, wind_sl)
    missing_therm = np.searchsorted(pres_sl, therm_sl)
    missing_wind = np.searchsorted(pres_sl, wind_sl)

    therm_sl = -np.ones(pres_sl.shape, dtype=pres_sl.dtype)
    therm_sl[missing_therm] = pres_sl[missing_therm]

    wind_sl = -np.ones(pres_sl.shape, dtype=pres_sl.dtype)
    wind_sl[missing_wind] = pres_sl[missing_wind]

    snd_filtered = {
        'temp':snd['temp'][pres_sl], #np.where(therm_sl == -1, missing, snd['temp'][therm_sl]), 
        'dewp':snd['dewp'][pres_sl], #np.where(therm_sl == -1, missing, snd['dewp'][therm_sl]),
        'pres':snd['pres'][pres_sl],
        'hght':snd['hght'][pres_sl],
        'wspd':np.where(wind_sl == -1, missing, snd['wspd'][wind_sl]),
        'wdir':np.where(wind_sl == -1, missing, snd['wdir'][wind_sl])
    }

    return snd_filtered
