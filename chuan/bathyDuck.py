# Snapshots of alongshore currents from vbar and adv measurements 

import glob, sys
sys.path.append('/home/chuan/Documents/Code/cmtb')
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib.backends.backend_pdf import PdfPages
from testbedutils import sblib as sb

advDir = 'https://chlthredds.erdc.dren.mil/thredds/dodsC/frf/projects/bathyduck/data/'
advFn = [
          'BathyDuck-ocean_currents_p11_201510.nc', # (124.8, 813.0)
        #   'BathyDuck-ocean_currents_p12_201510.nc', # (149.5, 814.0)
        #   'BathyDuck-ocean_currents_p13_201510.nc', # (200.1, 813.4)
        #   'BathyDuck-ocean_currents_p14_201510.nc', # (249.4, 813.0)
          'BathyDuck-ocean_currents_p21_201510.nc', # (124.9, 741.8)
        #   'BathyDuck-ocean_currents_p22_201510.nc', # (150.5, 740.2)
        #   'BathyDuck-ocean_currents_p23_201510.nc', # (199.9, 740.7)
        #   'BathyDuck-ocean_currents_p24_201510.nc', # (249.9, 741.3)
        #   'BathyDuck-ocean_currents_p83_201510.nc', # (199.2, 521.4)
        #   'BathyDuck-ocean_currents_p84_201510.nc', # (182.0, 521.8)
         ]
# advC = ['orange', 'green']
vbarDir = '/mnt/argus/argus_vbar2nc_transfer_test/vbar125/2015/10-Oct/'
windFn = 'http://134.164.129.55/thredds/dodsC/FRF/meteorology/wind/derived/derived.ncml'

nStations = len(advFn)
time = [None] * nStations
y = [None] * nStations
v = [None] * nStations

for i in range(nStations):
    with Dataset(advDir + advFn[i]) as rootgrp:
        time[i] = rootgrp['time'][:].astype('datetime64[s]')
        y[i] = rootgrp['alongVelYloc'][0]
        v[i] = rootgrp['alongVel'][:]

row = 0
col = 0

with Dataset(windFn) as rootgrp:
    timeWind = rootgrp['time'][:]
    wdspd = rootgrp['windSpeed'][:]

with PdfPages('bathyDuck.pdf') as pdf:
    for fn in sorted(glob.glob(vbarDir + '*.nc')):
        Y, M, D, h, m = fn[-18:-14], fn[-14:-12], fn[-12:-10], fn[-9:-7], fn[-7:-5]
        vbarDateTime = np.datetime64('{}-{}-{}T{}:{}'.format(Y, M, D, h, m), 's')
        plotted = False
        for i in range(nStations):
            ind = vbarDateTime == time[i]
            if not np.any(ind):
                continue
            v_ = v[i][ind][0]           
            if not plotted:
                with Dataset(fn) as rootgrp:
                    vbary = rootgrp['yFRF'][:]
                    vbarV = rootgrp['vMean'][:][0]
                    lowerCI = rootgrp['vLowerCI'][:][0]
                    upperCI = rootgrp['vUpperCI'][:][0]
                    QCspan = rootgrp['vQCspan'][:][0]
                if row == 0 and col == 0:
                    fig, ax = plt.subplots(nrows=3, ncols=2, figsize=(8.5, 11),
                                tight_layout=True, sharex=True, sharey=True)
                    ax[0, 0].set_ylim(-0.5, 0.5)
                p = ax[row, col].plot(vbary, vbarV)
                color = p[0].get_color()
                ax[row, col].plot(vbary, lowerCI, '--', alpha=0.5, c=color)
                ax[row, col].plot(vbary, upperCI, '--', alpha=0.5, c=color)
                _, wdspd_, _ = sb.timeMatch(timeWind, wdspd, 
                    [vbarDateTime.astype('float64')], vbarV)
                ax[row, col].text(0.05, 0.9, fn[-18:-5], 
                                  transform=ax[row, col].transAxes)
                if len(wdspd_) != 0:
                    ax[row, col].text(0.7, 0.9, 'wind = {:.3} m/s'.format(wdspd_[0]), 
                                    transform=ax[row, col].transAxes)
                ax2 = ax[row, col].twinx()
                ax2.plot(vbary, QCspan, c='y', alpha=0.5)
                plotted = True
            ax[row, col].plot(y[i], v_[0] / 100, '.', c='k')
        if plotted == True:
            if col == 1:
                if row == 2:
                    pdf.savefig()
                    row = 0
                else:
                    row += 1
                col = 0
            else: 
                col += 1
    if not (row == 0 and col == 0):
        pdf.savefig() 

# plt.show()