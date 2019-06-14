import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import sys, pickle
sys.path.append('/home/spike/repos')
from ocmLib import loadRSEX, combineStatsOCM
import numpy as np
import h5py, glob
import netCDF4 as nc
from testbedutils import sblib as sb
""" This file is step 3 in a work flow from openAndProcesss.py to ocm_compare.py to make plots
This work is to begin to understand the 'tuning' knobs of OCM methods, work being done in preparation for CIRN Meeting
June '19 
"""

### load appropriate data
fname = '/home/spike/repos/myOCM/figures/results/statistics_AllGauges.pickle'
with open(fname, 'rb') as fid:
    stats = pickle.load(fid)
# load obs
obsFileName= "/home/spike/repos/myOCM/waveCurStats_RSEX17.mat"
epochTime, xFRF, yFRF, measuredU, measuredV, Hs, depthBelowSurface = loadRSEX(obsFileName, idxG=4)

# fileList = glob.glob(
#     '/home/spike/repos/myOCM/data/processed/OCM_150m_{}_{}_Tstep{}_Twin{}_dy{}_dyWin{}_StackRes{}.mat'.format('*', '*',
#                                                                                                               step, win,
#                                                                                                               dy, dyWin,
#                                                                                                               res))
#
# ocmT, meanV, color1, color2, color3, QCspan_vals = combineOCMmats(fileList, yFRF, var1='prob', var2='QCspan', var3='SNR',
#                                                              var4='cispan')
# statOut[step][win][dy][dyWin][res][probFitThreshold] = stats
##################
# am i processing a wave staff
#
# gaugeY= 860
# for gauge in [150, 200, 250]:
#     flist = glob.glob('/mnt/gaia/peeler/argus/argus02b/2017/cx/*Oct*/*vbar{}.mat'.format(gauge))
#     print('doing gauge {}'.format(gauge))
#     # fname = '/mnt/gaia/peeler/argus/argus02b/2017/cx/287_Oct.14/1508020140.Sat.Oct.14_22_29_00.GMT.2017.argus02b.cx.vbar150.mat'
#     for fname in flist:
#         arrays = {}
#         f = h5py.File(fname, 'r')
#         for k, v in f.items():
#             arrays[k] = np.array(v)
#         f.close()
#
#         idxYs = np.argsort(arrays['XYZ'][1])
#         ys = arrays['XYZ'][1][idxYs]
#         raw = arrays['RAW'][idxYs, :]
#         time = nc.num2date(arrays['T'], 'seconds since 1970-01-01')
#
#         plt.figure()
#         plt.title('time stack @ x={} y={}'.format(gauge, gaugeY))
#         plt.pcolormesh(ys, time, raw.T, cmap='Greys')
#         plt.xlim([gaugeY-20, gaugeY+20])
#         plt.ylim([time[0], time[200]])
#         plt.xlabel('xFRF')
#         plt.ylabel('Time')
#         saveFname = '/home/spike/repos/myOCM/figures/rawImagery/TimeStackZoomed_gauge{}_{}'.format(gauge, time[0][0].strftime("%Y%m%dT%H%M%SZ"))
#         plt.savefig(saveFname); plt.close()






######## plot 1:  All comparisons RMSE
plt.figure(figsize=(12,8))
ax1 = plt.subplot(321)
ax1.set_title('Time Step [s]')
labels = []
for key in stats.keys():
    bias, rmse = combineStatsOCM(stats, stepWindow=[key])
    ax1.plot(np.tile(key, len(rmse)), rmse, 'r.', label='RMSE')
    labels.append(key)
ax1.set_ylabel('RMSE [m]')
ax1.set_xticks(labels)
ax1.set_xticklabels(labels)
#####################
ax2 = plt.subplot(323, sharey=ax1)
ax2.set_title('Time Window [s]')
labels = []
for key in stats[256].keys():
    bias, rmse = combineStatsOCM(stats, winWindow=[key])
    ax2.plot(np.tile(key, len(rmse)), rmse, 'r.', label='RMSE')
    labels.append(key)
ax2.set_ylabel('RMSE [m]')
ax2.set_xticks(labels)
ax2.set_xticklabels(labels)
## make box
box = 256
slushX = 3
slushY = 0.05
ax2.plot([box + slushX, box + slushX], [min(rmse) - slushY, max(rmse) + slushY], 'k-')
ax2.plot([box - slushX, box - slushX], [min(rmse) - slushY, max(rmse) + slushY], 'k-')
ax2.plot([box - slushX, box + slushX], [min(rmse) - slushY, min(rmse) - slushY], 'k-')
ax2.plot([box - slushX, box + slushX], [max(rmse) + slushY, max(rmse) + slushY], 'k-')
####################
ax3 = plt.subplot(322, sharey=ax1)
ax3.set_title('Alongshore Step [m]')
labels=[]
for key in stats[256][256].keys():
    bias, rmse = combineStatsOCM(stats, dyWindow=[key])
    ax3.plot(np.tile(key, len(rmse)), rmse, 'r.', label='RMSE')
    labels.append(key)
ax3.set_ylabel('RMSE [m]')
ax3.set_xticks(labels)
ax3.set_xticklabels(labels)
#####################
ax4 = plt.subplot(324, sharey=ax1)
ax4.set_title('Alongshore Window [m]')
labels = []
for key in stats[256][256][5].keys():
    bias, rmse = combineStatsOCM(stats, dyWinWindow=[key])
    ax4.plot(np.tile(key, len(rmse)), rmse, 'r.', label='RMSE')
    labels.append(key)
ax4.set_ylabel('RMSE [m]')
ax4.set_xticks(labels)
ax4.set_xticklabels(labels)
box = 10
slushX = .3
slushY = 0.05
ax4.plot([box + slushX, box + slushX], [min(rmse) - slushY, max(rmse) + slushY], 'k-')
ax4.plot([box - slushX, box - slushX], [min(rmse) - slushY, max(rmse) + slushY], 'k-')
ax4.plot([box - slushX, box + slushX], [min(rmse) - slushY, min(rmse) - slushY], 'k-')
ax4.plot([box - slushX, box + slushX], [max(rmse) + slushY, max(rmse) + slushY], 'k-')
###################################
ax5 = plt.subplot(326, sharey=ax1)
ax5.set_title('Stack interpolation resolution')
labels = []
for key in stats[256][256][5][10].keys():
    bias, rmse = combineStatsOCM(stats, resWindow=[key])
    ax5.plot(np.tile(key, len(rmse)), rmse, 'r.', label='RMSE')
    labels.append(key)
ax5.set_ylabel('RMSE [m]')
ax5.set_xticks(labels)
ax5.set_xticklabels(labels)
####################################################
ax6 = plt.subplot(325, sharey=ax1)
ax6.set_title('probFitWindow')
labels = []
for key in stats[256][256][5][10][0.5].keys():
    bias, rmse = combineStatsOCM(stats, probFitWindow=[key])
    ax6.plot(np.tile(key, len(rmse)), rmse, 'r.', label='RMSE')
    labels.append(key)
ax6.set_ylabel('RMSE [m]')
ax6.set_xticks(labels)
ax6.set_xticklabels(labels)
plt.tight_layout(rect=[0.02, 0.02, 0.98, 0.95])

########################################################################################################################
################################################################################################################################
################## refined by time  window
################################################################################################################################
plt.figure(figsize=(12,8))
ax1 = plt.subplot(321)
ax1.set_title('Time Step [s]')
labels = []
for key in stats.keys():
    bias, rmse = combineStatsOCM(stats, stepWindow=[key], winWindow=[256])
    ax1.plot(np.tile(key, len(rmse)), rmse, 'r.', label='RMSE')
    labels.append(key)
ax1.set_ylabel('RMSE [m]')
ax1.set_xticks(labels)
ax1.set_xticklabels(labels)
#####################
ax2 = plt.subplot(323, sharey=ax1)
ax2.set_title('Time Window [s]')
labels = []
for key in stats[256].keys():
    bias, rmse = combineStatsOCM(stats, winWindow=[256])
    ax2.plot(np.tile(key, len(rmse)), rmse, 'r.', label='RMSE')
    labels.append(key)
ax2.set_ylabel('RMSE [m]')
ax2.set_xticks(labels)
ax2.set_xticklabels(labels)
####################
ax3 = plt.subplot(322, sharey=ax1)
ax3.set_title('Alongshore Step [m]')
labels=[]
for key in stats[256][256].keys():
    bias, rmse = combineStatsOCM(stats, dyWindow=[key], winWindow=[256])
    ax3.plot(np.tile(key, len(rmse)), rmse, 'r.', label='RMSE')
    labels.append(key)
ax3.set_ylabel('RMSE [m]')
ax3.set_xticks(labels)
ax3.set_xticklabels(labels)
#####################
ax4 = plt.subplot(324, sharey=ax1)
ax4.set_title('Alongshore Window [m]')
labels = []
for key in stats[256][256][2.5].keys():
    bias, rmse = combineStatsOCM(stats, dyWinWindow=[key], winWindow=[256])
    ax4.plot(np.tile(key, len(rmse)), rmse, 'r.', label='RMSE')
    labels.append(key)
ax4.set_ylabel('RMSE [m]')
ax4.set_xticks(labels)
ax4.set_xticklabels(labels)
###################################
ax5 = plt.subplot(326, sharey=ax1)
ax5.set_title('Stack interpolation resolution')
labels = []
for key in stats[256][256][2.5][10].keys():
    bias, rmse = combineStatsOCM(stats, resWindow=[key], winWindow=[256])
    ax5.plot(np.tile(key, len(rmse)), rmse, 'r.', label='RMSE')
    labels.append(key)
ax5.set_ylabel('RMSE [m]')
ax5.set_xticks(labels)
ax5.set_xticklabels(labels)
####################################################
ax6 = plt.subplot(325, sharey=ax1)
ax6.set_title('probFitWindow')
labels = []
for key in stats[256][256][2.5][10][2].keys():
    bias, rmse = combineStatsOCM(stats, probFitWindow=[key], winWindow=[256])
    ax6.plot(np.tile(key, len(rmse)), rmse, 'r.', label='RMSE')
    labels.append(key)
ax6.set_ylabel('RMSE [m]')
ax6.set_xticks(labels)
ax6.set_xticklabels(labels)
plt.tight_layout(rect=[0.02, 0.02, 0.98, 0.95])
########################################################################################################################
################################################################################################################################
################## refined by time  window and Alongshore Window
################################################################################################################################
plt.figure(figsize=(12,8))
ax1 = plt.subplot(321)
ax1.set_title('Time Step [s]')
labels = []
for key in stats.keys():
    bias, rmse = combineStatsOCM(stats, stepWindow=[key], winWindow=[256], dyWinWindow=[30])
    ax1.plot(np.tile(key, len(rmse)), rmse, 'r.', label='RMSE')
    labels.append(key)
ax1.set_ylabel('RMSE [m]')
ax1.set_xticks(labels)
ax1.set_xticklabels(labels)
#####################
ax2 = plt.subplot(323, sharey=ax1)
ax2.axis('off')
ax2.set_title('Time Window [s]')
labels = []
for key in stats[256].keys():
    bias, rmse = combineStatsOCM(stats, winWindow=[key])
    ax2.plot(np.tile(key, len(rmse)), rmse, 'r.', label='RMSE')
    labels.append(key)
ax2.set_ylabel('RMSE [m]')
ax2.set_xticks(labels)
ax2.set_xticklabels(labels)
####################
ax3 = plt.subplot(322, sharey=ax1)
ax3.set_title('Alongshore Step [m]')
labels=[]
for key in stats[256][256].keys():
    bias, rmse = combineStatsOCM(stats, dyWindow=[key], winWindow=[256], dyWinWindow=[30])
    ax3.plot(np.tile(key, len(rmse)), rmse, 'r.', label='RMSE')
    labels.append(key)
ax3.set_ylabel('RMSE [m]')
ax3.set_xticks(labels)
ax3.set_xticklabels(labels)
#####################
ax4 = plt.subplot(324, sharey=ax1)
ax4.axis('off')
ax4.set_title('Alongshore Window [m]')
labels = []
for key in stats[256][256][2.5].keys():
    bias, rmse = combineStatsOCM(stats, dyWinWindow=[key], winWindow=[256])
    ax4.plot(np.tile(key, len(rmse)), rmse, 'r.', label='RMSE')
    labels.append(key)
ax4.set_ylabel('RMSE [m]')
ax4.set_xticks(labels)
ax4.set_xticklabels(labels)
###################################
ax5 = plt.subplot(326, sharey=ax1)
ax5.set_title('Stack interpolation resolution')
labels = []
for key in stats[256][256][2.5][10].keys():
    bias, rmse = combineStatsOCM(stats, resWindow=[key], winWindow=[256], dyWinWindow=[30])
    ax5.plot(np.tile(key, len(rmse)), rmse, 'r.', label='RMSE')
    labels.append(key)
ax5.set_ylabel('RMSE [m]')
ax5.set_xticks(labels)
ax5.set_xticklabels(labels)
####################################################
ax6 = plt.subplot(325, sharey=ax1)
ax6.set_title('probFitWindow')
labels = []
for key in stats[256][256][2.5][10][2].keys():
    bias, rmse = combineStatsOCM(stats, probFitWindow=[key], winWindow=[256], dyWinWindow=[30])
    ax6.plot(np.tile(key, len(rmse)), rmse, 'r.', label='RMSE')
    labels.append(key)
ax6.set_ylabel('RMSE [m]')
ax6.set_xticks(labels)
ax6.set_xticklabels(labels)
plt.tight_layout(rect=[0.02, 0.02, 0.98, 0.95])

########################################################################################################################
################################################################################################################################
################## refined by time  window and Alongshore Window, probibility of Fit
################################################################################################################################
plt.figure(figsize=(12,8))
ax1 = plt.subplot(321)
ax1.set_title('Time Step [s]')
labels = []
for key in stats.keys():
    bias, rmse = combineStatsOCM(stats, stepWindow=[key], winWindow=[256], dyWinWindow=[30], probFitWindow=[.75])
    ax1.plot(np.tile(key, len(rmse)), rmse, 'r.', label='RMSE')
    labels.append(key)
ax1.set_ylabel('RMSE [m]')
ax1.set_xticks(labels)
ax1.set_xticklabels(labels)
#####################
ax2 = plt.subplot(323)
ax2.axis('off')
ax2.set_title('Time Window [s]')
labels = []
for key in stats[256].keys():
    bias, rmse = combineStatsOCM(stats, winWindow=[key])
    ax2.plot(np.tile(key, len(rmse)), rmse, 'r.', label='RMSE')
    labels.append(key)
ax2.set_ylabel('RMSE [m]')
ax2.set_xticks(labels)
ax2.set_xticklabels(labels)
####################
ax3 = plt.subplot(322, sharey=ax1)
ax3.set_title('Alongshore Step [m]')
labels=[]
for key in stats[256][256].keys():
    bias, rmse = combineStatsOCM(stats, dyWindow=[key], winWindow=[256], dyWinWindow=[30], probFitWindow=[.75])
    ax3.plot(np.tile(key, len(rmse)), rmse, 'r.', label='RMSE')
    labels.append(key)
ax3.set_ylabel('RMSE [m]')
ax3.set_xticks(labels)
ax3.set_xticklabels(labels)
#####################
ax4 = plt.subplot(324)
ax4.axis('off')
ax4.set_title('Alongshore Window [m]')
labels = []
for key in stats[256][256][2.5].keys():
    bias, rmse = combineStatsOCM(stats, dyWinWindow=[key], winWindow=[256])
    ax4.plot(np.tile(key, len(rmse)), rmse, 'r.', label='RMSE')
    labels.append(key)
ax4.set_ylabel('RMSE [m]')
ax4.set_xticks(labels)
ax4.set_xticklabels(labels)
###################################
ax5 = plt.subplot(326, sharey=ax1)
ax5.set_title('Stack interpolation resolution')
labels = []
for key in stats[256][256][2.5][10].keys():
    bias, rmse = combineStatsOCM(stats, resWindow=[key], winWindow=[256], dyWinWindow=[30], probFitWindow=[.75])
    ax5.plot(np.tile(key, len(rmse)), rmse, 'r.', label='RMSE')
    labels.append(key)
ax5.set_ylabel('RMSE [m]')
ax5.set_xticks(labels)
ax5.set_xticklabels(labels)
####################################################
ax6 = plt.subplot(325)
ax6.set_title('probFitWindow')
ax6.axis('off')
labels = []
for key in stats[256][256][2.5][10][2].keys():
    bias, rmse = combineStatsOCM(stats, probFitWindow=[key], winWindow=[256], dyWinWindow=[30])
    ax6.plot(np.tile(key, len(rmse)), rmse, 'r.', label='RMSE')
    labels.append(key)
ax6.set_ylabel('RMSE [m]')
ax6.set_xticks(labels)
ax6.set_xticklabels(labels)
plt.tight_layout(rect=[0.02, 0.02, 0.98, 0.95])

# ################################################################################################################################
# ################## individual plots as above with bias
# ################################################################################################################################
# plt.figure()
# ax1 = plt.subplot(211)
# ax2 = plt.subplot(212)
# plt.suptitle('Time Step [s]')
# for key in stats.keys():
#     bias, rmse = combineStatsOCM(stats, stepWindow=[key])
#     ax1.plot(np.tile(key, len(rmse)), rmse, 'r.', label='RMSE')
#     ax2.plot(np.tile(key, len(bias)), bias, 'b.', label='Bias')
# ax2.plot(list(ax2.get_xlim()), [0,0], 'k')
# ax1.set_title('RMSE [m]')
# ax2.set_title('Bias [m]')
# plt.tight_layout(rect=[0.05, 0.05, 0.95, 0.95])
# #####################
# plt.figure()
# ax1 = plt.subplot(211)
# ax2 = plt.subplot(212)
# plt.suptitle('Time Window [s]')
# for key in stats[256].keys():
#     bias, rmse = combineStatsOCM(stats, winWindow=[key])
#     ax1.plot(np.tile(key, len(rmse)), rmse, 'r.', label='RMSE')
#     ax2.plot(np.tile(key, len(bias)), bias, 'b.', label='Bias')
# ax2.plot(list(ax2.get_xlim()), [0,0], 'k')
# ax1.set_title('RMSE [m]')
# ax2.set_title('Bias [m]')
# plt.tight_layout(rect=[0.05, 0.05, 0.95, 0.95])
# ####################
# plt.figure()
# ax1 = plt.subplot(211)
# ax2 = plt.subplot(212)
# plt.suptitle('Alongshore Step [m]')
# for key in stats[256][256].keys():
#     bias, rmse = combineStatsOCM(stats, dyWindow=[key])
#     ax1.plot(np.tile(key, len(rmse)), rmse, 'r.', label='RMSE')
#     # ax2.plot(np.tile(key, len(bias)), bias, 'b.', label='Bias')
# ax2.plot(list(ax2.get_xlim()), [0,0], 'k')
# ax1.set_title('RMSE [m]')
# ax2.set_title('Bias [m]')
# plt.tight_layout(rect=[0.05, 0.05, 0.95, 0.95])
# #####################
#
# plt.figure()
# ax1 = plt.subplot(211)
# ax2 = plt.subplot(212)
# plt.suptitle('Alongshore Window [m]')
# for key in stats[256][256][2.5].keys():
#     bias, rmse = combineStatsOCM(stats, dyWinWindow=[key])
#     ax1.plot(np.tile(key, len(rmse)), rmse, 'r.', label='RMSE')
#     ax2.plot(np.tile(key, len(bias)), bias, 'b.', label='Bias')
# ax2.plot(list(ax2.get_xlim()), [0,0], 'k')
# ax1.set_title('RMSE [m]')
# ax2.set_title('Bias [m]')
# plt.tight_layout(rect=[0.05, 0.05, 0.95, 0.95])
# ###################################
# plt.figure()
# ax1 = plt.subplot(211)
# ax2 = plt.subplot(212)
# plt.suptitle('Stack interpolation resolution')
# for key in stats[256][256][2.5][10].keys():
#     bias, rmse = combineStatsOCM(stats, resWindow=[key])
#     ax1.plot(np.tile(key, len(rmse)), rmse, 'r.', label='RMSE')
#     ax2.plot(np.tile(key, len(bias)), bias, 'b.', label='Bias')
# ax2.plot(list(ax2.get_xlim()), [0,0], 'k')
# ax1.set_title('RMSE [m]')
# ax2.set_title('Bias [m]')
# plt.tight_layout(rect=[0.05, 0.05, 0.95, 0.95])
# ####################################################
# plt.figure()
# ax1 = plt.subplot(211)
# ax2 = plt.subplot(212)
# plt.suptitle('probFitWindow')
# for key in stats[256][256][2.5][10][2].keys():
#     bias, rmse = combineStatsOCM(stats, probFitWindow=[key])
#     ax1.plot(np.tile(key, len(rmse)), rmse, 'r.', label='RMSE')
#     ax2.plot(np.tile(key, len(bias)), bias, 'b.', label='Bias')
# ax2.plot(list(ax2.get_xlim()), [0,0], 'k')
# ax1.set_title('RMSE [m]')
# ax2.set_title('Bias [m]')
# plt.tight_layout(rect=[0.05, 0.05, 0.95, 0.95])
#

