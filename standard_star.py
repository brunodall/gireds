#!/usr/bin/env python

#########################################################################
#                                                                       #
#   ATTENTION!!!                                                        #
#       DO NOT follow on a reduction process unless you are sure about  #
#       the fiber masks in the MDF file. Disregarding this warning will #
#       most certainly lead to major headaches at the final stages of   #
#       the reduction.                                                  #
#                                                                       #
#########################################################################

# Table of images

from pyraf import iraf
import numpy as np
import pyfits as pf
import glob
import time
from aplacos import applylacos

iraf.set(stdimage='imtgmos')

iraf.gemini()
iraf.gemtools()
iraf.gmos()

#iraf.unlearn('gemini')
#iraf.unlearn('gmos')

iraf.task(lacos_spec='/storage/work/gemini_pairs/lacos_spec.cl')

tstart = time.time()

#set directories
iraf.set(caldir='/datassd/gmos/GS-2013A-Q-61/')  # cal files, mainly the biases
iraf.set(rawdir='/datassd/gmos/GS-2013A-Q-61/')  # raw files
iraf.set(procdir='/datassd/gmos/star2013a/')  # processed files

iraf.cd('rawdir')

#iraf.setjd('*fits[0]')

l = glob.glob('*.fits')
l.sort()
idx = np.arange(len(l))

headers = [pf.getheader(i, ext=0) for i in l]
jd = [pf.getheader(i, ext=1)['mjd-obs'] for i in l]

# tolerance in hours
ttol = 2

starobj = 'Feige34'
stdstar = 'feige34'
caldir = 'onedstds$ctionewcal/'

star_idx = [i for i in idx if ((headers[i]['obstype'] == 'OBJECT')\
    &(headers[i]['obsclass'] == 'partnerCal')
    &(headers[i]['object'] != 'Twilight'))]

star = [l[i][:-5] for i in star_idx]

detector = [headers[i]['detector'] for i in star_idx]
observatory = [headers[i]['observat'] for i in star_idx]

arc = [l[i][:-5] for i in idx if\
    ((headers[i]['obstype'] == 'ARC') &\
    (headers[i]['date-obs'] == headers[star_idx[0]]['date-obs']))]

flat = [l[i][:-5] for i in idx if\
    ((headers[i]['obstype'] == 'FLAT') &\
    (abs(jd[i] - jd[star_idx[0]]) <= ttol*1./24.))]

twilight_idx = [i for i in idx if\
    ((headers[i]['obstype'] == 'OBJECT') & \
    (headers[i]['object'] == 'Twilight'))]

twilight = [l[i][:-5] for i in twilight_idx]
twilight_arc = [l[i][:-5] for i in idx if\
    ((headers[i]['obstype'] == 'ARC') & \
    (headers[i]['date-obs'] == headers[twilight_idx[0]]['date-obs']) & \
    (abs(jd[i] - jd[twilight_idx[0]]) <= ttol*1./24.))]

bias = [l[i] for i in idx if headers[i]['obstype'] == 'BIAS']




iraf.gmos.logfile='feige34.log'

iraf.cd('procdir')

# remove previous lists
# no need

# building lists

def range_string(l):
    return (len(l)*'{:4s},').format(*[i[-4:] for i in l])

iraf.gemlist(range=range_string(flat), root=flat[0][:-4], Stdout='flat.list')
iraf.gemlist(range=range_string(arc), root=arc[0][:-4], Stdout='arc.list')
iraf.gemlist(range=range_string(star), root=star[0][:-4], Stdout='star.list')
iraf.gemlist(range=range_string(twilight),
    root=twilight[0][:-4], Stdout='twilight.list')


iraf.gfreduce.bias = 'caldir$'+bias[0]

#######################################################################
#######################################################################
###   Star reduction                                                  #
#######################################################################
#######################################################################

#
#   Flat reduction
#

for i in ['g', 'rg', 'erg']:
    iraf.imdelete(i+'@flat.list')

iraf.gfreduce('@flat.list', slits='header', rawpath='rawdir$', fl_inter='no',
    fl_addmdf='yes', key_mdf='MDF', mdffile='default', weights='no',
    fl_over='yes', fl_trim='yes', fl_bias='yes', trace='yes', t_order=4,
    fl_flux='no', fl_gscrrej='no', fl_extract='yes', fl_gsappwave='no',
    fl_wavtran='no', fl_novl='no', fl_skysub='no', reference='',
    recenter='yes', fl_vardq='yes')

for i in ['g', 'rg', 'erg']:
    iraf.imdelete(i+'@twilight.list')

iraf.gfreduce('@twilight.list', slits='header', rawpath='rawdir$',
    fl_inter='no', fl_addmdf='yes', key_mdf='MDF',
    mdffile='default', weights='no',
    fl_over='yes', fl_trim='yes', fl_bias='yes', trace='yes', recenter='no',
    fl_flux='no', fl_gscrrej='no', fl_extract='yes', fl_gsappwave='no',
    fl_wavtran='no', fl_novl='no', fl_skysub='no',
    reference='erg'+flat[0], fl_vardq='yes')
#
#   Response function
#


for i, j in enumerate(flat):

    iraf.imdelete(j+'_response')
    iraf.gfresponse('erg'+j+'.fits', out='erg'+j+'_response',
        skyimage='erg'+twilight[i], order=95, fl_inter='no', func='spline3',
        sample='*', verbose='yes')

#   Arc reduction
#

for i in arc:
    iraf.imdelete('*'+i+'.fits')

iraf.gfreduce('@arc.list', slits='header', rawpath='rawdir$', fl_inter='no',
    fl_addmdf='yes', key_mdf='MDF', mdffile='default', weights='no',
    fl_over='yes', fl_trim='yes', fl_bias='yes', trace='no', recenter='no',
    fl_flux='no', fl_gscrrej='no', fl_extract='yes', fl_gsappwave='no',
    fl_wavtran='no', fl_novl='no', fl_skysub='no', reference='erg'+flat[0],
    fl_vardq='yes')


#   Finding wavelength solution
#   Note: the automatic identification is very good
#

for i in arc:
    iraf.gswavelength('erg'+i, function='chebyshev', nsum=15, order=4,
        fl_inter='no', nlost=5, ntarget=20, aiddebug='s', threshold=5,
        section='middle line')

#
#   Apply wavelength solution to the lamp 2D spectra
#

    iraf.imdelete('terg'+i+'.fits')
    iraf.gftransform('erg'+i, wavtran='erg'+i, outpref='t', fl_vardq='yes')

##
##   Actually reduce star
##

iraf.delete('std')
iraf.delete('sens.fits')

for i in star:
    iraf.delete('*g'+i+'.fits')
    
    iraf.gfreduce(
        i, slits='header', rawpath='rawdir$', fl_inter='no',
        fl_addmdf='yes', key_mdf='MDF', mdffile='default', weights='no',
        fl_over='yes', fl_trim='yes', fl_bias='yes', trace='no', recenter='no',
        fl_flux='no', fl_gscrrej='no', fl_extract='no', fl_gsappwave='no',
        fl_wavtran='no', fl_novl='yes', fl_skysub='no', fl_vardq='yes')

    #    applylacos('rg{:s}.fits'.format(i), clobber=True)
    iraf.gemcrspec('rg{:s}.fits'.format(i), out='lrg'+i, sigfrac=0.32, 
         niter=4, fl_vardq='yes')
         
    iraf.gfreduce(
        'lrg'+i+'.fits', slits='header', rawpath='./', fl_inter='no',
        fl_addmdf='no', key_mdf='MDF', mdffile='default',
        fl_over='no', fl_trim='no', fl_bias='no', trace='no', recenter='no',
        fl_flux='no', fl_gscrrej='yes', fl_extract='yes', fl_gsappwave='yes',
        fl_wavtran='yes', fl_novl='no', fl_skysub='yes',
        reference='erg'+flat[0], weights='no',
        wavtraname='erg'+arc[0], response='erg'+flat[0]+'_response.fits',
        fl_vardq='yes')
#
#   Apsumming the stellar spectra
#
    iraf.imdelete('astexlrg'+i)
    iraf.gfapsum(
        'stexlrg'+i, fl_inter='no', lthreshold=400.,
        reject='avsigclip')
#
#   Building sensibility function
#

iraf.delete('std')
iraf.delete('sens.fits')

iraf.gsstandard(
    (len(star)*'astexlrg{:s}.fits,').format(*star), starname=stdstar,
    observatory='Gemini-South', sfile='std', sfunction='sens', caldir=caldir)
#
#   Apply flux calibration to galaxy
#
#
##iraf.imdelete('cstexlrg@objr4.list')
#
##iraf.gscalibrate('stexlrg@objr4.list',sfunction='sens.fits',fl_ext='yes',extinct='onedstds$ctioextinct.dat',observatory='Gemini-South',fluxsca=1)
#
##
##   Create data cubes
##
#
#
##for i in objs:
##  iraf.imdelete('d0.1cstexlrg'+i+'.fits')
##  iraf.gfcube('cstexlrg'+i+'.fits',outpref='d0.1',ssample=0.1,fl_atmd='yes',fl_flux='yes')
#
##
## Combine cubes
##
#
#
##iraf.imdelete('am2306-721r4_wcsoffsets.fits')
##iraf.imcombine('d0.1cstexlrgS20141113S00??.fits[1]',output='am2306-721r4_wcsoffsets.fits',combine='average',reject='sigclip',masktype='badvalue',lsigma=2,hsigma=2,offset='wcs',outlimits='2 67 2 48 100 1795')
#

tend = time.time()

print('Elapsed time in reduction: {:.2f}'.format(tend - tstart))

