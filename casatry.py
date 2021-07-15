import os
import matplotlib.pyplot as plt
from astropy.io import fits
from casatasks import exportfits

from casatasks import tclean

filepath = os.path.relpath('3c391_ctm_mosaic_10s_spw0.ms')

print("running tclean, may take a bit...")


tclean(vis=filepath, imagename='try1', imsize=100, cell='10.0arcsec')

print("complete")

'''
filenames = ['first_image.image', 'first_image.mask', 'first_image.model', 
             'first_image.pb', 'first_image.psf', 'first_image.residual']

ff, aa = plt.subplots(2,3, figsize=(18,12))
for ii, name in enumerate(filenames):
  exportfits(imagename=name, fitsimage=name+'.fits', overwrite=True)
  xx,yy = int(ii/3),ii%3
  im = aa[xx,yy].imshow(fits.getdata(name+'.fits')[0,0,:,:])
  plt.colorbar(im, ax=aa[xx,yy])
  aa[xx,yy].set_title(name)
plt.show()
'''