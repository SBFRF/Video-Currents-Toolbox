import sys, glob
sys.path.append('/home/chuan/Documents/Code/cmtb')
import numpy as np
from netCDF4 import Dataset
from testbedutils import sblib as sb
from datetime import datetime

vbarDir = '/mnt/argus/argus_vbar2nc_transfer_test/vbar125/2015/10-Oct/'
windFn = 'http://134.164.129.55/thredds/dodsC/FRF/meteorology/wind/derived/derived.ncml'

with Dataset(windFn) as rootgrp:
    time = rootgrp['time'][:]
    wdspd = rootgrp['windSpeed'][:]

vbarFn = sorted(glob.glob(vbarDir + '*.nc'))[0]
Y, M, D, h, m = vbarFn[-18:-14], vbarFn[-14:-12], vbarFn[-12:-10], vbarFn[-9:-7], vbarFn[-7:-5]
vbarDateTime = np.datetime64('{}-{}-{}T{}:{}'.format(Y, M, D, h, m), 's')

with Dataset(vbarFn) as rootgrp:
    vbarV = rootgrp['vMean'][:][0]
    vbarTime = rootgrp['time'][:]

ind = time.astype('datetime64[s]') == np.datetime64('2015-10-13T11:29:59')

print(vbarDateTime.astype('float64'))
print(vbarTime)
print(time)
time, wdspd, vbarV = sb.timeMatch(time, wdspd, [vbarDateTime.astype('float64')], vbarV)
print(time.astype('datetime64[s]'))
print(wdspd)
