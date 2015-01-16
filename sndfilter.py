
import numpy as np

def _windSigLevels(**snd):
    idxs = np.arange(snd['wdir'].shape[0])
    np.random.shuffle(idxs)
    return np.sort(idxs[:20])

def _findInversions(temp, pres):
    return

def _findIsothermals(temp, pres, tol=0.5, min_depth=5):
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
                min_depth [type=float]
                    Minimum number of points to consider a "layer."  Default is
                    5 points.
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
        if len(cand_idxs) >= min_depth:
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

def _findAddSigLevels(temp, pres):
    return 

def _thermSigLevels(**snd):
    """
    _thermSigLevels()
    Purpose:    Find the temperature significant levels
    Keywords:   [ same as for the soundingFilter() function ]
    Returns:    A list of indices represending the indices of the significant
                    levels.
    """
    inv_sl = _findInversions(snd['temp'], snd['pres'])
    iso_sl = _findIsothermals(snd['temp'], snd['pres'])
    add_sl = _findAddSigLevels(snd['temp'], snd['pres'])

    therm_sl = np.concatenate((inv_sl, iso_sl, add_sl))
    return np.sort(therm_sl)

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
