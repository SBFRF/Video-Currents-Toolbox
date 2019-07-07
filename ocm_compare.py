"""this script is to take processed data as output from openAndProcess.py and compares to observations
using file waveCurStats_RSEX17.mat"""
import matplotlib
matplotlib.use('Agg')
import glob, pickle, sys
import numpy as np
from ocmLib import timeMatchOCM, loadRSEX, combineOCMmats
import netCDF4 as nc
sys.path.append('/home/spike/repos')
from testbedutils import sblib as sb
from matplotlib import pyplot as plt

#######################################################################

#######################################################################
obsFileName= "/home/spike/repos/myOCM/waveCurStats_RSEX17.mat"
epochTime, xFRF, yFRF, measuredU, measuredV, Hs, depthBelowSurface = loadRSEX(obsFileName, idxG=4)
############################## plot winds and waves for the day

# plt.figure()  ## figure if i have correct data parsed
# plt.plot(obsData['uv']['p'].squeeze()[idxG].squeeze(), '.', label='pressure')
# plt.plot(, '.', label='mab')
# plt.plot()


############## all the data
probFitThreshold = 0.75   ###### threshold used for removing data
mon = '*'
day = '*'
gauge = '*'
saveFname = 'statistics_AllGauges.pickle'
myCmap = 'Reds'
var1 = 'prob'
var2 = 'ci'
var3 = 'SNR'
var4 = 'QCspan'
statOut = {}
gaugeList = []

for step in [256, 128, 64,32][::-1]:
    statOut[step] = {}
    for win in [32, 64, 128, 256, 512]:
        statOut[step][win] = {}
        for dy in [5, 10]:
            statOut[step][win][dy] = {}
            for dyWin in [10, 20, 30]:
                statOut[step][win][dy][dyWin] = {}
                for res in [0.5]:
                    statOut[step][win][dy][dyWin][res] = {}
                    for probFitThreshold in [ 0.75, 0.9, 0.95]:
                        statOut[step][win][dy][dyWin][res][probFitThreshold] = {}
                        for QCspan in [20, 40, 60, 80, 100, 120]:
                            statOut[step][win][dy][dyWin][res][probFitThreshold][QCspan] = {}
                            for gauge in [150, 200]:
                                statOut[step][win][dy][dyWin][res][probFitThreshold][QCspan][gauge] = {}
                                fileList = sorted(glob.glob('/home/spike/repos/myOCM/data/processed/OCM_{}*_{}_{}_Tstep{}_Twin{}_dy{}*_dyWin{}*_StackRes{}*.mat'.format(gauge, mon, day, step, win, dy, dyWin, res)))
                                # how many gauges do i have
                                if gauge == 150:
                                    idxG = 4
                                elif gauge == 200:
                                    idxG = 5
                                epochTime, xFRF, yFRF, measuredU, measuredV, Hs, depthBelowSurface = loadRSEX(
                                    obsFileName, idxG=idxG)

                                ############################################# preprocess for plotting ##########################
                                # pre process data from multiple files for one plot
                                ocmT, meanV, color1, color2, color3, QCspan_vals = combineOCMmats(fileList, yFRF, var1=var1, var2 =var2, var3 =var3, var4 =var4)
                                ####################### Refine data after loaded ###############################################
                                # time match and Threshold, convert to masked arrays
                                color1 = np.ma.masked_array(color1, mask=np.isnan(color1))
                                color2 = np.ma.masked_array(color2, mask=np.isnan(color2))
                                color3 = np.ma.masked_array(color3, mask=np.isnan(color3))
                                QCspan_vals = np.ma.masked_array(QCspan_vals, mask=np.isnan(QCspan_vals))

                                ### refine time by probability of fit above some value
                                filterLogic = np.ma.masked_greater(color1, probFitThreshold).mask | np.ma.masked_less(QCspan_vals, QCspan).mask
                                meanV_plot = meanV[filterLogic]
                                ocmT_plot = ocmT[filterLogic]
                                color2_plot = color2[filterLogic]
                                color3_plot = color3[filterLogic]
                                color4_plot = QCspan_vals[filterLogic]
                                color1_plot = color1[filterLogic]  # this has to be done last becasue it removes comparison data

                                # filterLogic = np.ma.masked_less(color3_plot, SNRThreshold).mask
                                # meanV_plot = meanV_plot[filterLogic]
                                # ocmT_plot = ocmT_plot[filterLogic]
                                # color2_plot = color2_plot[filterLogic]
                                # color3_plot = color3_plot[filterLogic]
                                # color4_plot = color4_plot[filterLogic]
                                # color1_plot = color1_plot[filterLogic]  # this has to be done last becasue it removes comparison data


                                idxMatchedObs, idxMatchedOptical = timeMatchOCM(obsTime=epochTime, opticalTime=ocmT_plot)
                                if np.size(idxMatchedObs) < 2 or np.size(idxMatchedOptical) < 2:
                                    continue
                                meanV_plot = meanV_plot[idxMatchedOptical]
                                ocmT_plot = ocmT_plot[idxMatchedOptical]
                                color1_plot = color1_plot[idxMatchedOptical]
                                color2_plot = color2_plot[idxMatchedOptical]
                                color3_plot = color3_plot[idxMatchedOptical]
                                color4_plot = color4_plot[idxMatchedOptical]
                                measuredV_plot = measuredV[idxMatchedObs]
                                measuredV_plot_time = epochTime[idxMatchedObs]

                                stats = sb.statsBryant(observations=measuredV_plot, models=meanV_plot)
                                ############################################# begin figure #########################################
                                myAlpha = 1
                                scatterMarkerSize=25
                                capsize = 2 # caps for the error bar
                                ebarMarkerSize = 100
                                ebarLineWidth = 1
                                ebarAlpha = 0.1
                                ebarLineStyle = ':'
                                myCmap = 'Reds'
                                figsize=(12, 8)
                                ##########################################################################################
                                fig = plt.figure(figsize=figsize)
                                fig.suptitle('Comparisons at xFRF={} yFRF={}\n tStep={} tWin={}, yStep={} yWin={}\np'
                                             'lotting only data with probility of Fit > {} QCspan < {}'.format(xFRF, yFRF,
                                                                           step,win, dy, dyWin, probFitThreshold, QCspan))
                                ax0 = plt.subplot2grid((7, 8), (0,0), colspan=8, rowspan=2)
                                ax0.plot(nc.num2date(epochTime, 'seconds since 1970-01-01'), measuredV, '.', label='ADV', ms=2)
                                ax0.errorbar(nc.num2date(ocmT_plot, 'seconds since 1970-01-01'), meanV_plot, color2_plot.T, label='OCM', linestyle='', capsize=capsize, alpha=ebarAlpha)
                                ax0.set_xlim([min(nc.num2date(ocmT_plot, 'seconds since 1970-01-01')), max(nc.num2date(ocmT_plot, 'seconds since 1970-01-01'))])
                                ax0.set_ylabel('velocity $[m/s]$')
                                ax0.legend()
                                ax0.set_ylim([-4, 4])


                                ax1 = plt.subplot2grid((7, 8), (2,0), colspan=4, rowspan=4)
                                ax1.errorbar(measuredV_plot, meanV_plot, color2_plot.T, ms=ebarMarkerSize, mfc='Reds',
                                             label='OCM', linestyle='', capsize=capsize, elinewidth=ebarLineWidth, alpha=ebarAlpha)
                                cbar1 = ax1.scatter(measuredV_plot, meanV_plot, c=color1_plot, alpha=myAlpha, cmap=myCmap, s=scatterMarkerSize)
                                cbar11 = plt.colorbar(cbar1)
                                cbar11.set_label(var1)
                                ax1.plot([-2, 2], [-2, 2], '-k')
                                ax1.set_ylabel('Ocm Velocity [m/s]')
                                ax1.set_ylim([-4, 4])
                                ax1.set_xlim([-4, 4])

                                ax3 = plt.subplot2grid((7, 8), (2,4), colspan=4, rowspan=4)
                                ax3.errorbar(measuredV_plot, meanV_plot, color2_plot.T, ms=ebarMarkerSize, mfc='Reds',
                                             label='OCM', linestyle='', capsize=capsize, elinewidth=ebarLineWidth, alpha=ebarAlpha)
                                cbar3 = ax3.scatter(measuredV_plot, meanV_plot, c=color4_plot, alpha=myAlpha, cmap=myCmap, s=scatterMarkerSize) #, norm=matplotlib.colors.LogNorm()
                                cbar22 = plt.colorbar(cbar3)
                                cbar22.set_label(var4)
                                ax3.plot([-2, 2], [-2, 2], '-k')
                                ax3.set_ylim([-4,4])
                                ax3.set_xlim([-4,4])

                                # ax3.set_ylabel('Ocm Velocity [m/s]')

                                # do stats here below
                                ax2 = plt.subplot2grid((7, 8), (6, 0), colspan=8, rowspan=1)
                                ax2.axis('off')
                                ax2.text(0.5, 0.5, 'bias = {:.2f} $[m/s]$\nRMSE = {:.2f} $[m/s]$\n$R^2$ = {:.2f}\nsymetric Slope = {:.2f}\nsample Count = {}'.format(
                                    stats['bias'],stats['RMSE'],stats['r2'],stats['symSlope'], len(stats['residuals'])), fontsize=12,
                                         horizontalalignment='center', verticalalignment='top')


                                plt.tight_layout(rect=[0.02, 0.05, .98, 0.95])
                                ofname = (
                                    "/home/spike/repos/myOCM/figures/results/whichQC_Tstep{}_Twin{}_dyStep{}_dyWin{}_QCspan{}_prob{}_xFRF{}_yFRF{}.png".format(
                                        step, win, dy, dyWin, QCspan, probFitThreshold, xFRF, yFRF))
                                plt.savefig(ofname); plt.close()
                                statOut[step][win][dy][dyWin][res][probFitThreshold][QCspan][gauge] = stats



print('yay!')
# save somestuff
fname = '/home/spike/repos/myOCM/figures/results/'+ saveFname
with open(fname, 'wb') as fid:
    pickle.dump(statOut, fid, protocol=pickle.HIGHEST_PROTOCOL)


# make plot with QC span (or velocity span in the OCM measurments) against obs to see how well the bounding boxes
# are made
# ax2 = plt.subplot(423, sharex=ax1, sharey=ax1)
# # ax2.set_title("Only QC data")
# cbar2 = ax2.scatter(measuredV_plot, meanV_plot, c=color2_plot, alpha=myAlpha, cmap = myCmap, s=scatterMarkerSize)
# cbar22 = plt.colorbar(cbar2)
# cbar22.set_label(var2)
# ax2.plot([-2, 2], [-2, 2], '-k')
# ax2.set_ylabel('Ocm Velocity [m/s]')
# ax22 = plt.subplot(424)
# ax22.hist(color2_plot, bins=40)
#
# ax3 = plt.subplot(425, sharex=ax1, sharey=ax1)
# # ax3.set_title("Only QC data")
# cbar3 = ax3.scatter(measuredV_plot, meanV_plot, c=color3_plot, alpha=myAlpha, cmap=myCmap, s=scatterMarkerSize)
# cbar33 = plt.colorbar(cbar3)
# cbar33.set_label(var3)
# ax3.plot([-2, 2], [-2, 2], '-k')
# ax3.set_ylabel('Ocm Velocity [m/s]')
# ax3.set_ylabel('OCM measurements')
# ax33 = plt.subplot(426)
# ax33.hist(color3_plot, bins=40)
#
# ax4 = plt.subplot(427, sharex=ax1, sharey=ax1)
# # ax4.set_title("Only QC data")
# cbar4 = ax4.scatter(measuredV_plot, meanV_plot, c=color4_plot, alpha=myAlpha, cmap=myCmap, s=scatterMarkerSize)
# cbar44 = plt.colorbar(cbar4)
# cbar44.set_label(var4)
# ax4.plot([-2, 2], [-2, 2], '-k')
# ax4.set_xlabel('observations')
# ax4.set_ylabel('Ocm Velocity [m/s]')
# ax44 = plt.subplot(428)
# ax44.hist(color4_plot, bins=40)
# ax44.set_xlabel('value')
# fig.suptitle("Tstep{}_Twin{}_dy{}_yWindow{}_StackRes{}.mat".format(step, win, dy, dyWin, res))

########################################################################################################################