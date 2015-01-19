
import numpy as np
from scipy.interpolate import interp1d
from pylab import *

def _centered2nd(prof):
    return prof[2:] - 2 * prof[1:-1] + prof[:-2]

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

def _windSigLevels(**snd):
    max_wind_sl = _maxWind(snd['wspd'], snd['wdir'], snd['pres'])
    wind_spd_sl = _splitProfile(snd['wspd']*0.514444, 0, 5)
    return np.concatenate((max_wind_sl, wind_spd_sl))


def groupLayers(cand_layers):
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

def _findInversions(temp, dewp, pres, trop_idx):
    lr_layers = np.diff(temp)
    cand_layers = np.where(lr_layers > 0)[0]
    groupedLayers = groupLayers(cand_layers)
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
#   print inversion_layers
    return np.asarray(inversion_layers, dtype=int)

def _maxWind(wspd, wdir, pres, min_wspd=(30 * 1.94), min_diff=(10 * 1.94), smooth_pts=11):
    '''
    _maxWind()
    Purpose:    Find the level of the maximum winds in the profile
    Arguments:  the snd dictionary.
    Returns:    The index from the sounding of the maximum wind speed.

    Assumes that the winds are in knots.
    '''

    if np.ma.is_masked(wspd):
        wspd_mask = wspd.mask
        wspd = wspd[np.where(~wspd_mask)]

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

    def groupLayers(cand_layers):
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

            cand_idxs = groupLayers(cand_idxs)
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
    cand_idxs = groupLayers(cand_layers)

    for cand_idx in cand_idxs:
        cand_idx = pareLayer(cand_idx)

        for cidx in cand_idx:
            isotherm_idxs.extend([cidx[0], cidx[-1]])

    return isotherm_idxs

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

def _findAddSigLevels(temp, pres, trop_idx, tol=[1.0, 2.0]):

    cutoff_idx = min(np.argmin(np.abs(pres - 300)), trop_idx)

    idxs_below = _splitProfile(temp[:cutoff_idx], 0, tol[0])[:-1]
    idxs_above = _splitProfile(temp[cutoff_idx:], cutoff_idx, tol[1])
    additional_idxs = np.sort(np.unique(np.concatenate((idxs_below, idxs_above))))

    return additional_idxs

def _RH(temp, dwpt, pres):
    return 100*(np.exp((17.625*dwpt)/(243.04+dwpt))/np.exp((17.625*temp)/(243.04+temp))) 

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

def _thermSigLevels(**snd):
    """
    _thermSigLevels()
    Purpose:    Find the temperature significant levels
    Keywords:   [ same as for the soundingFilter() function ]
    Returns:    A list of indices represending the indices of the significant
                    levels.
    """

    trop_idx = _findTropopause(**snd)

    inv_sl = _findInversions(snd['temp'], snd['dewp'], snd['pres'], trop_idx)
    iso_sl = _findIsothermals(snd['temp'], snd['pres'])
    add_sl = _findAddSigLevels(snd['temp'], snd['pres'], trop_idx)

    therm_sl = np.concatenate((inv_sl, iso_sl, add_sl))
    return np.sort(therm_sl)

def reverseDiff(array):
    return array[1:] + array[:-1]

def _highestLevel(**snd):
    return np.ma.argmin(snd['pres'])

def _findTropopause(**snd):
    '''
        Finding the tropopause level using the WMO definition of a 
        tropopause as being the lowest level where the 2 km layer
        aloft has a lapse rate greater than 2 C/km.
        
        The algorithm below was taken from:
        http://www.inscc.utah.edu/~reichler/publications/papers/2003_tropo.pdf
    '''
    above_500 = np.where((snd['pres'] > 75) & (snd['pres'] < 550))[0]
    dT = np.diff(snd['temp'][above_500])
    dp = np.diff(snd['pres'][above_500])
    at = reverseDiff(snd['temp'][above_500])
    ap = reverseDiff(snd['pres'][above_500])
    ah = reverseDiff(snd['hght'][above_500])
    R = 287.
    g = 9.81
    Cp = 1004.
    kappa = R / Cp

    halfp = ap / 2.
    halfh = ah / 2.
    halfgamma = (dT/dp) * (ap/at) * (kappa * g/R)
    idx = np.where(halfgamma > 0.002)[0]
    for i in idx:
        hbot = halfh[i]
        htop = hbot + 2000.
        htop_idx = np.argmin(np.abs(halfh - htop))
        mean_lapse = np.mean(halfgamma[i:htop_idx])
        if mean_lapse > 0.002:
            return np.ma.argmin(np.fabs(snd['pres']-halfp[i]))    
    return np.nan

def soundingFilter(**snd):
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

    from datetime import datetime

    pres_ml = _mandatoryPresLevels(snd['pres'])
    trop_idx = _findTropopause(**snd)
    therm_sl = _thermSigLevels(**snd)
    rh_sl = _findSigRHLevels(snd['temp'], snd['dewp'], snd['pres'])
    wind_sl = _windSigLevels(**snd)

    all_idxs = np.concatenate((pres_ml, [trop_idx], therm_sl, rh_sl, wind_sl))
    all_idxs = np.unique(all_idxs[np.isfinite(all_idxs)])
    all_idxs = np.asarray(all_idxs[all_idxs >= 0], dtype=int)

    """
    ticks = np.concatenate((np.arange(1000,0,-100), [50,30, 20, 10]))
    subplot(121)
    gca().set_yscale('log')
    ylabel("Pressure [mb]")
    xlabel("Temperature [C]")
    plot(snd['temp'], snd['pres'], 'r-')
    plot(snd['dewp'], snd['pres'], 'g-')
    plot(snd['temp'][all_idxs], snd['pres'][all_idxs], 'ko')
    plot(snd['dewp'][all_idxs], snd['pres'][all_idxs], 'ko') 
    gca().invert_yaxis()
    yticks(ticks, np.asarray(ticks, dtype=str))
    subplot(122)
    gca().set_yscale('log')
    xlabel("Wind Speed [kts]")
    plot(snd['wspd'], snd['pres'], 'b-')
    plot(snd['wspd'][all_idxs], snd['pres'][all_idxs], 'ko')
    yticks(ticks, np.asarray(ticks, dtype=str))
    gca().invert_yaxis()
    show()
    """

    pres_sl = np.union1d(therm_sl, wind_sl)
    missing_therm = np.searchsorted(pres_sl, therm_sl)
    missing_wind = np.searchsorted(pres_sl, wind_sl)

    therm_sl = -np.ones(pres_sl.shape, dtype=pres_sl.dtype)
    therm_sl[missing_therm] = pres_sl[missing_therm]

    wind_sl = -np.ones(pres_sl.shape, dtype=pres_sl.dtype)
    wind_sl[missing_wind] = pres_sl[missing_wind]

    snd_filtered = {
        'temp':np.where(therm_sl == -1, np.nan, snd['temp'][therm_sl]), 
        'dewp':np.where(therm_sl == -1, np.nan, snd['dewp'][therm_sl]),
        'pres':snd['pres'][pres_sl],
        'wspd':np.where(wind_sl == -1, np.nan, snd['wspd'][wind_sl]),
        'wdir':np.where(wind_sl == -1, np.nan, snd['wdir'][wind_sl])
    }

    return snd_filtered
