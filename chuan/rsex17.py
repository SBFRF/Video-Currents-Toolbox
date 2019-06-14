import sys
import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
sys.path.append('/home/chuan/Documents/Code/cmtb')
from datetime import datetime
from getdatatestbed.getDataFRF import getObs
from glob import glob
from netCDF4 import Dataset
from testbedutils.sblib import timeMatch

plt.rc('date.autoformatter', day='%b %d')

startTime = '2017-9-22T00:00:00Z'
endTime = '2017-10-20T00:00:00Z'
waveStation = 'waverider-26m'

advDir = '/home/chuan/Documents/Misc/RSEX17/'
advFn = 'waveCurStats_RSEX17.mat'

# advStations = 'hm22', 'hm32'
# vbarDir = ['/mnt/argus/argus_vbar2nc_transfer_test/vbar200/2017/10-Oct/',]
# loc = 200
advStations = 'hm31',
vbarDir = ['/mnt/argus/argus_vbar2nc_transfer_test/vbar150/2017/09-Sep/',
           '/mnt/argus/argus_vbar2nc_transfer_test/vbar150/2017/10-Oct/',]
loc = 150

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
axTS[0].set_title('2017 Oct - {} m away from shore'.format(loc), fontsize=12)
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
axTS[2].set_ylabel('Alongshore current (m/s)')

figOne2one, axOne2one = plt.subplots(tight_layout=True, figsize=(7.5, 7.5))

m = sio.loadmat(advDir + advFn)

vbarFiles = []
for dr in vbarDir:
    vbarFiles += glob(dr + '*.nc')

for station in advStations:
    f = m[station][0, 0]
    dn = f['dn'][0]
    tAdv = ((dn - datetime(1971,1,2).toordinal()) * 60 * 60 * 24
            ).astype('datetime64[s]').astype(datetime)
    vAdv = f['measuredV_plot'][0]
    yAdv = f['y'][0][0]

    vVbar = [None] * 1000
    tVbar = [None] * 1000
    sVbar = [None] * 1000
    i = 0

    for fn in vbarFiles:
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
axOne2one.set_title('2017 Oct - {} m away from shore'.format(loc), fontsize=12)
axOne2one.set_ylabel('Optical alongshore current (m)')
axOne2one.set_xlabel('In situ alongshore current (m)')

cbar = figOne2one.colorbar(cax, shrink=0.5)
cbar.ax.set_ylabel('Skill')

figTS.autofmt_xdate()

figTS.savefig('rsex17_ts_{}.png'.format(loc))
figOne2one.savefig('rsex17_one2one_{}.png'.format(loc))
plt.show()