import matplotlib.pyplot as plt
import numpy as np
from glob import glob
from netCDF4 import Dataset
from matplotlib.backends.backend_pdf import PdfPages
from datetime import datetime

advDir = 'https://chlthredds.erdc.dren.mil/thredds/dodsC/frf/projects/bathyduck/data/'

# vbarDir = '/mnt/argus/argus_vbar2nc_transfer_test/vbar125/2015/10-Oct/'
# advFn = ['BathyDuck-ocean_currents_p11_201510.nc', # (124.8, 813.0)
#          'BathyDuck-ocean_currents_p21_201510.nc'] # (124.9, 741.8)

# vbarDir = '/mnt/argus/argus_vbar2nc_transfer_test/vbar150/2015/10-Oct/'
# advFn = ['BathyDuck-ocean_currents_p12_201510.nc', # (149.5, 814.0)
#          'BathyDuck-ocean_currents_p22_201510.nc'] # (150.5, 740.2)

vbarDir = '/mnt/argus/argus_vbar2nc_transfer_test/vbar200/2015/10-Oct/'
advFn = ['BathyDuck-ocean_currents_p13_201510.nc', # (200.1, 813.4)
         'BathyDuck-ocean_currents_p23_201510.nc', # (199.9, 740.7)
         'BathyDuck-ocean_currents_p83_201510.nc'] # (199.2, 521.4)

advN = len(advFn)
tAdv = [None] * advN
yAdv = [None] * advN
vAdv = [None] * advN

for i in range(advN):
    with Dataset(advDir + advFn[i]) as rootgrp:
        tAdv[i] = rootgrp['time'][:]
        yAdv[i] = rootgrp['alongVelYloc'][0]
        vAdv[i] = rootgrp['alongVel'][:, 0] / 100

with PdfPages('timeStack.pdf') as pdf:
    for fn in glob(vbarDir + '*.nc'):
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