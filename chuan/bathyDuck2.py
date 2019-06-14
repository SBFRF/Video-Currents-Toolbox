import glob, sys
import matplotlib.pyplot as plt
import numpy as np
sys.path.append('/home/chuan/Documents/Code/cmtb')
from datetime import datetime
from getdatatestbed.getDataFRF import getObs
from netCDF4 import Dataset
from testbedutils.sblib import timeMatch

plt.rc('font', size=12)
plt.rc('date.autoformatter', hour='%b %d %Hh')
plt.rc('date.autoformatter', day='%b %d')

startTime = '2015-10-13T00:00:00Z'
endTime = '2015-10-18T00:00:00Z'
waveStation = 'waverider-26m'

advDir = 'https://chlthredds.erdc.dren.mil/thredds/dodsC/frf/projects/bathyduck/data/'

advFn = ['BathyDuck-ocean_currents_p11_201510.nc', # (124.8, 813.0)
         'BathyDuck-ocean_currents_p21_201510.nc'] # (124.9, 741.8)
vbarDir = '/mnt/argus/argus_vbar2nc_transfer_test/vbar125/2015/10-Oct/'
loc = 125

# advFn = ['BathyDuck-ocean_currents_p12_201510.nc', # (149.5, 814.0)
#          'BathyDuck-ocean_currents_p22_201510.nc'] # (150.5, 740.2)
# vbarDir = '/mnt/argus/argus_vbar2nc_transfer_test/vbar150/2015/10-Oct/'
# loc = 150

# advFn = ['BathyDuck-ocean_currents_p13_201510.nc', # (200.1, 813.4)
#          'BathyDuck-ocean_currents_p23_201510.nc', # (199.9, 740.7)
#          'BathyDuck-ocean_currents_p83_201510.nc'] # (199.2, 521.4)
# vbarDir = '/mnt/argus/argus_vbar2nc_transfer_test/vbar200/2015/10-Oct/'
# loc = 200

startTime = datetime.strptime(startTime, '%Y-%m-%dT%H:%M:%SZ')
endTime = datetime.strptime(endTime, '%Y-%m-%dT%H:%M:%SZ')

go = getObs(startTime, endTime)
wo = go.getWaveSpec(waveStation)
wdo = go.getWind()

figTS, axTS = plt.subplots(3, sharex=True, tight_layout=True, figsize=(7.5, 7.5),
                           gridspec_kw={'height_ratios': [1, 1, 2]})
axTS[0].plot(wo['time'], wo['Hs'], label='Wave height')
axTS[0].tick_params('y', colors='#1f77b4')
axTS[0].set_ylabel('Wave height (m)', color='#1f77b4')
axTS[0].set_title('2015 Oct - {} m away from shore'.format(loc))
h, l = axTS[0].get_legend_handles_labels()
ax_ = axTS[0].twinx()
ax_.plot(wo['time'], wo['waveDm'], c='#ff7f0e', label='Mean wave direction', ls='--')
ax_.tick_params('y', colors='#ff7f0e')
ax_.set_ylabel('Wave direction (deg)', color='#ff7f0e')
ax_.set_ylim(0, 360)
h_, l_ = ax_.get_legend_handles_labels()
ax_.legend(h + h_, l + l_, ncol=2, fontsize=12)
axTS[1].plot(wdo['time'], wdo['windspeed'], label='Wind speed')
axTS[1].tick_params('y', colors='#1f77b4')
axTS[1].set_ylabel('Wind speed (m/s)', color='#1f77b4')
h, l = axTS[1].get_legend_handles_labels()
ax_ = axTS[1].twinx()
ax_.plot(wdo['time'], wdo['winddir'], c='#ff7f0e', label='Wind direction', ls='--')
ax_.set_ylabel('Wind direction (deg)', color='#ff7f0e')
ax_.set_ylim(0, 360)
ax_.tick_params('y', colors='#ff7f0e')
h_, l_ = ax_.get_legend_handles_labels()
ax_.legend(h + h_, l + l_, ncol=2, fontsize=12)
axTS[2].set_ylim(-0.5, 0.5)
axTS[2].set_ylabel('Alongshore current (m/s)')
figTS.autofmt_xdate()

figOne2one, axOne2one = plt.subplots(tight_layout=True, figsize=(7.5, 7.5))

for fn in advFn:
    with Dataset(advDir + fn) as rootgrp:
        tAdv = rootgrp['time'][:].astype('datetime64[s]').astype(datetime)
        yAdv = rootgrp['alongVelYloc'][0]
        vAdv = rootgrp['alongVel'][:, 0] / 100

    vAdv = vAdv[np.logical_and(tAdv >= startTime, tAdv <= endTime)]
    tAdv = tAdv[np.logical_and(tAdv >= startTime, tAdv <= endTime)]
    
    vVbar = [None] * 1000
    tVbar = [None] * 1000
    sVbar = [None] * 1000
    i = 0
    for fn in sorted(glob.glob(vbarDir + '*.nc')):
        with Dataset(fn) as rootgrp:
            yVbar = rootgrp['yFRF'][:]
            ind = np.argmin(np.abs(yVbar - yAdv))
            valid = rootgrp['vValid'][:][0][ind]
            vVbar[i] = rootgrp['vMean'][:][0][ind]
            tVbar[i] = datetime.utcfromtimestamp(rootgrp['time'][:][0])
            sVbar[i] = rootgrp['vSkill'][:][0][ind]
            if valid == 0:
                lInval, = axTS[2].plot(tVbar[i], vVbar[i], '.b', fillstyle='none', 
                                       zorder=10)
            else:
                lVal, = axTS[2].plot(tVbar[i], vVbar[i], 'xk', zorder=10)
        i += 1
    tVbar = tVbar[:i]
    vVbar = vVbar[:i]
    axTS[2].plot(tAdv, vAdv, label='In situ at y = {} m'.format(yAdv))

    tAdv, vVbarTM, vAdvTM = timeMatch(tVbar, vVbar, tAdv, vAdv)
    tAdv, sVbarTM, vAdvTM = timeMatch(tVbar, sVbar, tAdv, vAdv)
    cax = axOne2one.scatter(vAdvTM, vVbarTM, c=sVbarTM)
    axOne2one.plot([-0.1, 0.2], [-0.1, 0.2], c='k')

lInval.set_label('Optical (invalid)')
lVal.set_label('Optical (valid)')
axTS[2].legend(ncol=2, fontsize=12)

axOne2one.set_xlim(-0.1, 0.2)
axOne2one.set_ylim(-0.1, 0.2)
axOne2one.set_aspect('equal')
axOne2one.set_title('2015 Oct - {} m away from shore'.format(loc), fontsize=12)
axOne2one.set_ylabel('Optical alongshore current (m)')
axOne2one.set_xlabel('In situ alongshore current (m)')

cbar = figOne2one.colorbar(cax, shrink=0.5)
cbar.ax.set_ylabel('Skill')

figTS.savefig('bathyDuck_ts_{}.png'.format(loc))
figOne2one.savefig('bathyDuck_one2one_{}.png'.format(loc))
plt.show()