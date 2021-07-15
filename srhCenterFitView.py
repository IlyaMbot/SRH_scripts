#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 31 09:41:28 2021

@author: sergeyvlesovoi
"""


#from BadaryRAO import BadaryRAO
#from astropy import coordinates
#from sunpy import coordinates as sunpy_coordinates
#from skimage.transform import warp, AffineTransform
#from astropy import constants
#from zeep import Client
#from skimage import measure
import numpy as NP
import pylab as PL;
from astropy.io import fits
import os, fnmatch, sys

def arcmin_format(xy, pos):
    global gArcsecPerPixel
    global gXSize
    return '%2d' % ((xy - gXSize/2) * gArcsecPerPixel / 60);

def findFits(path, pattern):
    result = [];
    for root, dirs, files in os.walk(path):
        for basename in files:
            if fnmatch.fnmatch(basename,pattern):
                if os.path.getsize(os.path.join(root,basename)) > 2880:
                    result.append(os.path.join(root,basename));
    return result;

#fitPath = 'SRH36/cleanMaps/202106/2800/'
#fitPath = 'SRH36/cleanMaps/20210708/'
fitPath = '/home/mariagloba/Pictures/SRH36_pngs/20210714/'
##

iFitNames =  findFits(fitPath + 'I/','*.fit')
vFitNames =  findFits(fitPath + 'V/','*.fit')
iFitNames.sort()
vFitNames.sort()
iImages = []
vImages = []
freqList = []
gXSize = 0
gArcsecPerPixel = 0
gDateObs = []
for fileName in iFitNames:
    fd = fits.open(fileName)
    iImages.append(NP.flipud(fd[0].data))
    freqList.append(fd[0].header['FREQUENC'])
    gXSize = int(fd[0].header['NAXIS1'])
    gArcsecPerPixel = fd[0].header['CDELT1']
    gDateObs.append(fd[0].header['DATE-OBS'].split('.')[0])
    fd.close()
      
for fileName in vFitNames:
    fd = fits.open(fileName)
    vImages.append(NP.flipud(fd[0].data))
    freqList.append(fd[0].header['FREQUENC'])
    gXSize = int(fd[0].header['NAXIS1'])
    gArcsecPerPixel = fd[0].header['CDELT1']
    gDateObs.append(fd[0].header['DATE-OBS'].split('.')[0])
    fd.close()

T0 = 2e5
fig, pl = PL.subplots(figsize=(6,6))
fig.tight_layout()
for i in range(len(iFitNames)):
    levels = NP.linspace(T0,NP.max(iImages[i]),10)
    pl.clear()
    pl.xaxis.set_major_locator(PL.MultipleLocator(128));
    pl.xaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
    pl.xaxis.set_minor_locator(PL.MultipleLocator(32));
    pl.yaxis.set_major_locator(PL.MultipleLocator(128));
    pl.yaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
    pl.yaxis.set_minor_locator(PL.MultipleLocator(32));
    pl.imshow(iImages[i],cmap='hot',vmax=T0,vmin=0)
    pl.contour(iImages[i],cmap='hot',levels=levels,linewidths=0.5)
    pl.set_title('SRH %s,  %s MHz'%(gDateObs[i], freqList[i]))
    pl.set_xlabel('arcmin')
    pl.set_ylabel('arcmin')
    fig.savefig(iFitNames[i].split('.')[0] + '.png')

T0 = 1e5
fig, pl = PL.subplots(figsize=(6,6))
fig.tight_layout()
for i in range(len(vFitNames)):
    pl.clear()
    pl.xaxis.set_major_locator(PL.MultipleLocator(128));
    pl.xaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
    pl.xaxis.set_minor_locator(PL.MultipleLocator(32));
    pl.yaxis.set_major_locator(PL.MultipleLocator(128));
    pl.yaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
    pl.yaxis.set_minor_locator(PL.MultipleLocator(32));
    pl.imshow(vImages[i],cmap='gray',vmax=T0,vmin=-T0)
    levels = NP.linspace(T0,NP.max(vImages[i]),3)
    pl.contour(vImages[i],cmap='gray',levels=levels,linewidths=0.5)
    levels = NP.linspace(NP.min(vImages[i]),-T0,3)
    pl.contour(vImages[i],cmap='gray',levels=levels,linewidths=0.5)
    pl.set_title('SRH %s,  %s MHz'%(gDateObs[i], freqList[i]))
    pl.set_xlabel('arcmin')
    pl.set_ylabel('arcmin')
    fig.savefig(vFitNames[i].split('.')[0] + '.png')
