
import numpy as np

import Nio as nio

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

def _findAddSigLevels(temp, pres):
    return

def _thermSigLevels(**snd):
    inv_sl = _findInversions(snd['temp'], snd['pres'])
    iso_sl = _findIsothermals(snd['temp'], snd['pres'])
    add_sl = _findAddSigLevels(snd['temp'], snd['pres'])

    therm_sl = np.concatenate((inv_sl, iso_sl, add_sl))
    return np.sort(therm_sl)

def soundingFilter(**snd):
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
    name_translate = {
        'tdry':'temp',
        'dp':'dewp',
        'pres':'pres',
        'wspd':'wspd',
        'deg':'wdir'
    }

    base_path = "/data6/tsupinie/20110524/arm/"
    snd_file = nio.open_file("%s/sgpsondewnpnS05.b1.20110524.202900.cdf" % base_path, mode='r')
    snd = dict( (name_translate[k], snd_file.variables[k][:]) for k in ['tdry', 'dp', 'pres', 'wspd', 'deg'])
    snd_filtered = soundingFilter(**snd)
