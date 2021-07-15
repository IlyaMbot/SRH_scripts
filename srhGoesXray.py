#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 02:08:41 2021

@author: svlesovoi
"""

import json
import numpy as NP
import pylab as PL
from astropy.io import fits
from scipy import signal

def hhmm_format(t, pos):
  hh = (int)(t / 3600.);
  t -= hh*3600.;
  mm = (int)(t / 60.);
  return '%02d:%02d' % (hh,mm);

win = signal.windows.gaussian(2199,2155)
rainbowColors = PL.get_cmap('rainbow')

#goesFile = open('xrays-6-hour.json')
goesFile = open('xrays-1-day.json')
goesData = json.load(goesFile)
goesFile.close()

N = len(goesData)
xrays_time = NP.zeros(N//2)
xrays_4_8 = NP.zeros(N//2)
xrays_005_04 = NP.zeros(N//2)

for i in range(N):
    if i%2:
        xrays_4_8[i//2] = goesData[i]['flux']
        hhmm = goesData[i]['time_tag'].split('T')[1].split(':')[0:2]
        xrays_time[i//2] = 3600*int(hhmm[0]) + 60*int(hhmm[1])
        if xrays_time[i//2] > 10*3600:
            xrays_time[i//2] -= 24*3600
    else:
        xrays_005_04[i//2] = goesData[i]['flux']

cF = fits.open('srh_cp_20210412.fits')
srhFreqList = cF[1].data['frequencies']
srhTime = cF[2].data['time'] # necc
srhCorrI = cF[2].data['I']
srhCorrV = cF[2].data['V']
srhFluxI = cF[2].data['flux_I'] # necc
srhMeanFluxI = srhFluxI.mean(axis=0)
srhMeanFluxISmoothed = signal.convolve(srhMeanFluxI,win,mode='same')/win.sum()

t0 = 200

fig = PL.figure()
sub = fig.add_subplot(1,1,1);
sub.set_ylabel('flux');
sub.set_xlabel('UT');
sub.xaxis.set_major_locator(PL.MultipleLocator(1800));
sub.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
sub.xaxis.set_minor_locator(PL.MultipleLocator(600));
sub.set_xlim(3600,6.0*3600)
sub.set_ylim(0,5e-7)

sub.plot(xrays_time[t0:],xrays_4_8[t0:],label='GOES X-Ray 0.1-0.8 nm',color='red',markersize=0.2)
sub.plot(xrays_time[t0:],xrays_005_04[t0:],label='GOES X-Ray 0.05-0.4 nm',color='blue',markersize=0.2)
for freq in range(srhFreqList.shape[0]):
    sub.plot(srhTime[freq],srhCorrI[freq]*1e-4,'.',markersize=0.2,color=rainbowColors(100+(srhFreqList.shape[0] - freq)*20),label='SRH %d MHz'%(srhFreqList[freq]*1e-3))
#    sub.plot(srhTime[freq],srhCorrV[freq]*1e-4)
#sub.plot(srhTime[0],(srhMeanFluxI - srhMeanFluxISmoothed)*5e-9 + 1e-8)
sub.plot([3600,10*3600],[1e-7,1e-7], label='X-ray flare class A')
sub.grid()
sub.legend(markerscale=50)
sub.set_title('SRH and GOES incredible coincidence , %s'%(cF[0].header['DATE-OBS']))
