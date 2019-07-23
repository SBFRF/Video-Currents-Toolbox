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

def timeMatchOCM(obsTime, opticalTime, allData=True):
    """
    function looks though optical time for the closest value that is below the median time step of the optical time (tWin)
    the obs sample period.  It does this so multiple OCM measurements can be compared to the observations
    Args:
        obsTime:  epochtime (or some other numeric)
        opticalTime: epoch time (or some other numeric) -- this is master time to be matched to
        allData (bool): if True will return as many opticalTimes as possible, if False, will return
            data that are within half of the observation sample window. eg if obs samples are ever 8 minutes (tWin = 8*60)
            function will return data that are within a 4 minute window (default =True)

    Returns:
        idxObs: indices of observations
        idxOptical: indices of optical

    """
    idxOptical, idxObs = [], []
    if allData == True:
        obsSamplePeriod = np.median(np.diff(obsTime))
    else:
        obsSamplePeriod = np.median(np.diff(obsTime))/2
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
                                    0, 0].squeeze() - yFRF)).squeeze()     # find which yLocation of vBar output we're interested in for comparing to observations
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
        stats: made by ocm_compare.py Assumed structure of stats[Tstep][Twin][yStep][yWin][ciSpan][probFit][QCspan][gauge]
        kwargs:
         'stepWindow'[step]:  Tstep - time length to step the window (in points)
         'winWindow'[win]:   Twin - the time length of the FFT window (in points)
         'dyWindow'[dy]:    output resolution in alongshore (how frequently to output in alongshore)
         'dyWinWindow'[dyWin]: how wide alongshore window is
         'resWindow'[res]:   resolution of interpolation scheme in prepossessing
         'probFitWindow'[probFitThreshold]:  filter applied to data
         'qcSpanWindow': the span of the confidence interval
         'stat1' (str): some key from sblib.statsBryant corresponding to some statistic of interest
         'stat2' (str): some key from sblib.statsBryant corresponding to some statistic of interest

    Returns
        stat1, stat2

    """
    stepWindow= kwargs.get('stepWindow', stats.keys())
    winWindow =  kwargs.get('winWindow', stats[list(stepWindow)[0]].keys())
    dyWindow =  kwargs.get('dyWindow', stats[list(stepWindow)[0]][list(winWindow)[0]].keys())
    dyWinWindow = kwargs.get('dyWinWindow', stats[list(stepWindow)[0]][list(winWindow)[0]][list(dyWindow)[0]].keys())
    ciSpanWindow = kwargs.get('resWindow', stats[list(stepWindow)[0]][list(winWindow)[0]][list(dyWindow)[0]][list(dyWinWindow)[0]].keys())
    probFitWindow = kwargs.get('probFitWindow', stats[list(stepWindow)[0]][list(winWindow)[0]][list(dyWindow)[0]][
                                        list(dyWinWindow)[0]][list(ciSpanWindow)[0]].keys())
    QCspanWindow = kwargs.get('qcSpanWindow',
                               stats[list(stepWindow)[0]][list(winWindow)[0]][list(dyWindow)[0]][list(dyWinWindow)[0]][
                                   list(ciSpanWindow)[0]][list(probFitWindow)[0]].keys())
    gaugeList = kwargs.get('gaugeList', stats[list(stepWindow)[0]][list(winWindow)[0]][list(dyWindow)[0]][list(dyWinWindow)[0]][
                                   list(ciSpanWindow)[0]][list(probFitWindow)[0]][list(QCspanWindow)[0]].keys())
    if not isinstance(gaugeList, list):
        gaugeList = [gaugeList]

    stat1 = kwargs.get('stat1', 'bias')
    stat2 = kwargs.get('stat2', 'RMSE')
    stat1Out, stat2Out = [], []
    for step in stepWindow:
        for win in winWindow:
            for dy in dyWindow:
                for dyWin in dyWinWindow:
                    for ciSpan in ciSpanWindow:
                        for probFit in probFitWindow:
                            for QCspan in QCspanWindow:
                                for gauge in gaugeList:
                                    if gauge is not 'all':
                                        gauge= int(gauge)
                                    if gauge == 'all' and 'all' in stats[int(step)][int(win)][int(dy)][int(dyWin)][ciSpan][probFit][int(QCspan)].keys():
                                        statDict = stats[int(step)][int(win)][int(dy)][int(dyWin)][ciSpan][probFit][int(QCspan)][gauge]
                                        stat1Out.append(statDict[stat1])
                                        stat2Out.append(statDict[stat2])
                                    elif gauge == 'all' and 'all' not in stats[int(step)][int(win)][int(dy)][int(dyWin)][ciSpan][probFit][int(QCspan)].keys():
                                        for key2 in stats[int(step)][int(win)][int(dy)][int(dyWin)][ciSpan][probFit][int(QCspan)].keys():
                                            if len(stats[int(step)][int(win)][int(dy)][int(dyWin)][ciSpan][probFit][int(QCspan)][key2]) > 0:
                                                gauge = key2
                                            if gauge is not 'all':  # if its still all, there is no data
                                                statDict = stats[int(step)][int(win)][int(dy)][int(dyWin)][ciSpan][probFit][int(QCspan)][
                                                        gauge]
                                                stat1Out.append(statDict[stat1])
                                                stat2Out.append(statDict[stat2])
                                    else: # normal gauges
                                        try:
                                            statDict = stats[int(step)][int(win)][int(dy)][int(dyWin)][ciSpan][probFit][int(QCspan)][gauge]
                                        # np.size(statDict.keys()) > 1:
                                            stat1Out.append(statDict[stat1])
                                            stat2Out.append(statDict[stat2])
                                        except KeyError:
                                            continue
    return stat1Out, stat2Out