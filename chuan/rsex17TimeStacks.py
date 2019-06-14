import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from glob import glob
from netCDF4 import Dataset
from matplotlib.backends.backend_pdf import PdfPages
from datetime import datetime

advDir = '/home/chuan/Documents/Misc/RSEX17/'
advFn = 'waveCurStats_RSEX17.mat'

# vbarDir = ['/mnt/argus/argus_vbar2nc_transfer_test/vbar150/2017/09-Sep/',
#            '/mnt/argus/argus_vbar2nc_transfer_test/vbar150/2017/10-Oct/',]
# advStations = 'hm31',

vbarDir = ['/mnt/argus/argus_vbar2nc_transfer_test/vbar200/2017/10-Oct/',]
advStations = 'hm22', 'hm32'


advN = len(advStations)
tAdv = [None] * advN
yAdv = [None] * advN
vAdv = [None] * advN

m = sio.loadmat(advDir + advFn)
for i in range(advN):
    f = m[advStations[i]][0, 0]
    dn = f['dn'][0]
    tAdv[i] = (dn - datetime(1971,1,2).toordinal()) * 60 * 60 * 24
    vAdv[i] = f['measuredV_plot'][0]
    yAdv[i] = f['y'][0][0]

vbarFiles = []
for dr in vbarDir:
    vbarFiles += glob(dr + '*.nc')

with PdfPages('timeStack.pdf') as pdf:
    for fn in vbarFiles:
        print(fn.split('/')[-1])
        with Dataset(fn) as rootgrp:
            tVbar = rootgrp['time'][:][0]
            tVbarRaw = rootgrp['rawTime'][:]
            iVbarRaw = rootgrp['rawIntensity'][:]
            yVbarRaw = rootgrp['rawY'][:]
            vVbar = rootgrp['vMean'][:][0]
            yVbar = rootgrp['yFRF'][:]
            validVbar = rootgrp['vValid'][:][0]

        fig, ax = plt.subplots(2, tight_layout=True, sharex=True, figsize=(11, 8.5))
        ax[0].imshow(iVbarRaw.T, aspect='auto', extent=(yVbarRaw[0], yVbarRaw[-1],
                     tVbarRaw[-1], tVbarRaw[0]), cmap='gray')
        ax[0].set_title(fn.split('/')[-1])
        ax[0].set_ylabel('Time (s)')
        ax[1].scatter(yVbar[validVbar == 1], vVbar[validVbar == 1], c='b', s=2,
                      label='Optical (valid)')
        ax[1].scatter(yVbar[validVbar == 0], vVbar[validVbar == 0], c='r', s=2,
                      label='Optical (invalid)')
        ax[1].set_ylim(-0.5, 0.5)
        ax[1].set_xlabel('Alongshore position (m)')
        ax[1].set_ylabel('Alongshore velocity (m/s)')

        exist = False
        for i in range(advN):
            ind = np.argmin(np.abs(tAdv[i] - tVbar))
            if np.abs(tAdv[i][ind] - tVbar) < 3600:
                exist = True
                l1, = ax[1].plot(yAdv[i], vAdv[i][ind], 'xk')
                l2 = ax[0].axvline(yAdv[i], c='lime', ls='--')
        if exist:
            l1.set_label('In situ')
            l2.set_label('In situ location')
            ax[0].legend()
        ax[1].legend()

        pdf.savefig(fig)
        plt.close(fig)