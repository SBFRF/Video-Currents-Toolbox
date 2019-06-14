import numpy as np
from scipy.io import loadmat
import warnings, sys, shutil
sys.path.append('/home/spike/repos')
from testbedutils import sblib as sb

def loadRSEX(obsFileName, idxG):
    """
    loads RSEX data from pats mat file
    Args:
       obsFileName: file name
       idxG:  which gauge to pull from

    Returns:
        epochTime: time in
        xFRF: frf cross-shore position
        yFRF: frf alonghsore position
        measuredU: cross-shore velocity
        measuredV: alongshore velocity
        Hs : wave height
        depthBelowSurface: how far below water surface is gauge
    """
    obsData = loadmat(obsFileName)
    epochTime = sb.mtime2epoch(obsData['uv']['dn'].squeeze()[idxG].squeeze())
    xFRF = obsData['uv']['x'].squeeze()[idxG].squeeze()
    yFRF = obsData['uv']['y'].squeeze()[idxG].squeeze()
    measuredU = obsData['uv']['u'].squeeze()[idxG].squeeze()
    measuredV = obsData['uv']['v'].squeeze()[idxG].squeeze()
    Hs = obsData['uv']['Hs'].squeeze()[idxG].squeeze()
    offset = (obsData['uv']['pmab'].squeeze()[idxG] - obsData['uv'].squeeze()['mab'][idxG]).squeeze()  # positive values mean gauge is pointed downwards
    depthBelowSurface = obsData['uv']['p'].squeeze()[idxG].squeeze() - obsData['uv']['mab'].squeeze()[idxG].squeeze() - offset
    warnings.warn('check Depth below Surface value before using')

    return epochTime, xFRF, yFRF, measuredU, measuredV, Hs, depthBelowSurface

def timeMatchOCM(obsTime, opticalTime):
    """
    function looks though optical time for the closest value that is below the time step of the optical time
    the obs sample period.  It does this so multiple OCM measurements can be compared to the observations
    Args:
        obsTime:  epochtime (or some other numeric)
        opticalTime: epoch time (or some other numeric) -- this is master time to be matched to

    :return:
    """
    idxOptical, idxObs = [], []
    obsSamplePeriod = np.median(np.diff(obsTime))
    # opticalSamplePeriod = np.median(np.diff(opticalTime))
    for idxOT, oTime in enumerate(opticalTime):
        idxMaybe = np.argmin(np.abs(oTime - obsTime))
        if np.abs(oTime - obsTime[idxMaybe]) < obsSamplePeriod:
            idxObs.append(idxMaybe)
            idxOptical.append(idxOT)
           #print('obs[idx]{} optical {} diff {}'.format(obsTime[idxMaybe], oTime, np.abs(obsTime[idxMaybe] - oTime)))
    return idxObs, idxOptical


def combineOCMmats(fileList, yFRF, var1='prob', var2 = 'QCspan', var3 = 'SNR', var4 = 'cispan'):
    """
    loads multiple OCM files as processed by matlab wrapper
    Args
        fileList:
        yFRF: pulls/combines locations of only y location of interest
        var1: any variable in the data file (as produced from matlab scripts)
        var2: any variable in the data file (as produced from matlab scripts)
        var3: any variable in the data file (as produced from matlab scripts)
        var4: any variable in the data file (as produced from matlab scripts)

    Returns
        concatenated values across multiple files

    """
    assert np.size(fileList) > 0, 'No files in file list'
    for fname in fileList:  # pre-process to one array
        data = loadmat(fname)
        idxY = np.argmin(np.abs(data['dataSave']['y'][
                                    0, 0].squeeze() - yFRF)).squeeze()  # find which column we're interested in for comparing to observations
        if fname == fileList[0]:
            meanV = data['dataSave']['meanV'][0, 0].squeeze()[:, idxY]
            ocmT = data['dataSave']['t'][0,0].squeeze()[:,idxY]
            color1 = data['dataSave'][var1][0, 0].squeeze()[:, idxY]
            if var2 == 'ci':
                try:
                    color2 = data['dataSave'][var2][0, 0].squeeze()[:, :, idxY]
                except IndexError:
                    print('found Bad Shaped CI in file {}'.format(fname))
                    shutil.move(fname, '/home/spike/repos/myOCM/data/processed/badCIs')
            else:
                color2 = data['dataSave'][var2][0, 0].squeeze()[:, idxY]
            color3 = data['dataSave'][var3][0, 0].squeeze()[:, idxY]
            color4 = data['dataSave'][var4][0, 0].squeeze()[:, idxY]
        else: # append
            meanV = np.append(meanV, data['dataSave']['meanV'][0, 0].squeeze()[:, idxY])
            ocmT = np.append(ocmT, data['dataSave']['t'][0,0].squeeze()[:,idxY] )
            color1 = np.append(color1, data['dataSave'][var1][0, 0].squeeze()[:, idxY])
            if var2 == 'ci':
                try:
                    color2 = np.append(color2, data['dataSave'][var2][0, 0].squeeze()[:, :, idxY], axis=0)
                except IndexError:
                    print('found Bad Shaped CI in file {}'.format(fname))
                    shutil.move(fname, '/home/spike/repos/myOCM/data/processed/badCIs')
            else:
                color2 = data['dataSave'][var2][0, 0].squeeze()[:, idxY]
            color3 = np.append(color3, data['dataSave'][var3][0, 0].squeeze()[:, idxY])
            color4 = np.append(color4, data['dataSave'][var4][0, 0].squeeze()[:, idxY])

    return ocmT, meanV, color1, color2, color3, color4


def combineStatsOCM(stats, **kwargs):
    """
    combines stats from nested statistics dictionary (made by ocm_compare.py), pulls only RMSE and bias right now

    Args:
        stats: made by ocm_compare.py
        kwargs:
         'stepWindow'[step]:  Tstep - time length to step the window (in points)
         'winWindow'[win]:   Twin - the time length of the FFT window (in points)
         'dyWindow'[dy]:    output resolution in alongshore (how frequently to output in alongshore)
         'dyWinWindow'[dyWin]: how wide alongshore window is
         'resWindow'[res]:   resolution of interpolation scheme in preprocessing
         'probFitWindow'[probFitThreshold]:  filter applied to data

    Returns
        bias, RMSE

    """
    stepWindow= kwargs.get('stepWindow', stats.keys())
    winWindow =  kwargs.get('winWindow', stats[list(stepWindow)[0]].keys())
    dyWindow =  kwargs.get('dyWindow', stats[list(stepWindow)[0]][list(winWindow)[0]].keys())
    dyWinWindow = kwargs.get('dyWinWindow', stats[list(stepWindow)[0]][list(winWindow)[0]][list(dyWindow)[0]].keys())
    resWindow = kwargs.get('resWindow', stats[list(stepWindow)[0]][list(winWindow)[0]][list(dyWindow)[0]][list(dyWinWindow)[0]].keys())
    probFitWindow = kwargs.get('probFitWindow', stats[list(stepWindow)[0]][list(winWindow)[0]][list(dyWindow)[0]][list(dyWinWindow)[0]][list(resWindow)[0]].keys())
    bias, RMSE = [], []
    for step in stepWindow:
        for win in winWindow:
            for dy in dyWindow:
                for dyWin in dyWinWindow:
                    for res in resWindow:
                        for probFit in probFitWindow:
                            bias.append(stats[step][win][dy][dyWin][res][probFit]['bias'])
                            RMSE.append(stats[step][win][dy][dyWin][res][probFit]['RMSE'])
    return bias, RMSE