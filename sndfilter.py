
import numpy as np

def _windSigLevels(**snd):
    idxs = np.arange(snd['wdir'].shape[0])
    np.random.shuffle(idxs)
    return np.sort(idxs[:20])

def _findInversions(temp, pres):
    return

def _findIsothermals(temp, pres, tol=0.5):
    lr_layers = np.abs(np.diff(temp))
    cand_layers = np.where(lr_layers < 2 * tol)
    return

def _mandatoryPresLevels(pres, tol=0.1):
    '''
        Searches for the surface, 1000, 925, 850, 700,
        500, 400, 300, 250, 200, 150, 100, 70, 50, and 10 mb.
    '''
    mandatory_pres = np.asarray([1000, 925, 850, 700, 500, 400, 300, 250, 200,\
                      150, 100, 70, 50, 10])
    mandatory_pres_idx = np.ones(mandatory_pres.shape, dtype=int)
    pres = np.round(pres,0)
    for i in xrange(len(mandatory_pres_idx)):
        ind = np.where(np.fabs(pres - mandatory_pres[i]) < tol)[0]
        if len(ind) == 0:
            ind = [-9999]
        mandatory_pres_idx[i] = ind[0]

    return mandatory_pres_idx

def _findAddSigLevels(temp, pres):
    return

def _thermSigLevels(**snd):
    inv_sl = _findInversions(snd['temp'], snd['pres'])
    iso_sl = _findIsothermals(snd['temp'], snd['pres'])
    add_sl = _findAddSigLevels(snd['temp'], snd['pres'])

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
    pres_ml = _mandatoryPresLevels(snd['pres'])
    _findTropopause(**snd)
    therm_sl = _thermSigLevels(**snd)
    
    wind_sl = _windSigLevels(**snd)

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

if __name__ == "__main__":
    import Nio as nio
    from pylab import *
    name_translate = {
        'tdry':'temp',
        'dp':'dewp',
        'pres':'pres',
        'wspd':'wspd',
        'deg':'wdir',
        'alt':'hght'
    }

    snd_file = nio.open_file('sgpnwsupa1sX11.b1.20130520.165157.cdf', mode='r')
    snd_file = nio.open_file('sgpnwsupa1sX11.b1.20130520.110839.cdf', mode='r')
    snd = dict( (name_translate[k], snd_file.variables[k][:]) for k in ['tdry', 'dp', 'pres', 'wspd', 'deg', 'alt'])
    snd_filtered = soundingFilter(**snd)
