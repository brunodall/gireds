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
#from aplacos import applylacos
import matplotlib.pyplot as plt

# Read starinfo.dat and define structured array
infofile = '/home/eslamba/Dropbox/ic/reducao/starinfo.dat'
nstar = sum(1 for line in open(infofile))
infoname = ['obj', 'stdcal', 'caldir', 'flux', 'esoflux']
infofmt = ['|S25','|S25', '|S25', '|S25', '|S25']
ninfo = len(infoname)
starinfo = np.zeros(nstar,  dtype={'names':infoname, 'formats':infofmt})

with open(infofile, 'r') as arq:
    for i in range(nstar):
        linelist = arq.readline().split()
        for j in range(ninfo):
            starinfo[i][j] = linelist[j]


#set directories
iraf.set(caldir='/home/eslamba/ic/reducao/GN-2014B-Q-87_GS-2013A-Q-56/')  # cal files, mainly the biases
iraf.set(rawdir='/home/eslamba/ic/reducao/GN-2014B-Q-87_GS-2013A-Q-56/')  # raw files
iraf.set(procdir='/home/eslamba/ic/reducao/NGC5728/procSTAR_15-01-16/')  # processed files
locrep = '/home/eslamba/Dropbox/ic/reducao/15-01-16/ctionewcal_flux/' # star flux files

iraf.cd('rawdir')

iraf.set(stdimage='imtgmos')

iraf.gemini()
iraf.gemtools()
iraf.gmos()

iraf.unlearn('gemini')
iraf.unlearn('gmos')

iraf.gmos.logfile='ltt7379.log'
iraf.task(lacos_spec='/home/eslamba/LAcosmic/lacos_spec.cl')

tstart = time.time()


#iraf.setjd('*fits[0]')

l = glob.glob('*.fits')
l.sort()
idx = np.arange(len(l))

headers = [pf.getheader(i, ext=0) for i in l]
jd = [pf.getheader(i, ext=1)['mjd-obs'] for i in l]


#   Define lists
# flat tolerance in hours
ttol = 2
# boas toletance in dias
ttol_bias = 100

# Estrela nao esta como uma lista de listas
star_idx = [i for i in idx if ((headers[i]['obstype'] == 'OBJECT')\
    &(headers[i]['obsclass'] == 'partnerCal')
    &(headers[i]['object'] != 'Twilight'))]

star = [l[i][:-5] for i in star_idx]

arc = [[l[i][:-5] for i in idx if\
       ((headers[i]['obstype'] == 'ARC') &\
        (headers[i]['date-obs'] == headers[j]['date-obs'])&\
        (headers[i]['observat'] == headers[j]['observat'])&\
        (headers[i]['centwave'] == headers[j]['centwave']))]
       for j in star_idx]

flat = [[l[i][:-5] for i in idx if\
        ((headers[i]['obstype'] == 'FLAT') &\
         (abs(jd[i] - jd[j]) <= ttol*1./24.)&\
         (headers[i]['observat'] == headers[j]['observat'])&\
         (headers[i]['centwave'] == headers[j]['centwave']))]
        for j in star_idx]

twilight = [[l[i][:-5] for i in idx if\
            ((headers[i]['obstype'] == 'OBJECT') & \
             (headers[i]['object'] == 'Twilight')&\
             (headers[i]['observat'] == headers[j]['observat'])&\
             (headers[i]['centwave'] == headers[j]['centwave']))]
            for j in star_idx]
 
bias = [[l[i][:-5] for i in idx if\
        ((headers[i]['obstype'] == 'BIAS')&\
         (abs(jd[i] - jd[j]) <= ttol_bias)&\
         (headers[i]['observat'] == headers[j]['observat']))]
        for j in star_idx]

detector = [headers[i]['detector'] for i in star_idx]
observatory = [headers[i]['observat'] for i in star_idx]


# Get star info
starinfo_idx = [ [i for i, m in enumerate(starinfo['obj']) \
                 if m==headers[j]['object']][0] for j in star_idx]

starobj = [starinfo[i]['obj'] for i in starinfo_idx]
stdstar = [starinfo[i]['stdcal'] for i in starinfo_idx]
caldir = [starinfo[i]['caldir'] for i in starinfo_idx]
starflux = [[starinfo[i]['flux'] for i in starinfo_idx]]


iraf.cd('procdir')

# remove previous lists
# no need

# building lists

def range_string(l):
    return (len(l)*'{:4s},').format(*[i[-4:] for i in l])

for j in range(len(star)):
    iraf.gemlist(range=range_string(flat[j]), root=flat[j][0][:-4], 
        Stdout='flat'+str(j)+'.list')
    iraf.gemlist(range=range_string(arc[j]), root=arc[j][0][:-4],
        Stdout='arc'+str(j)+'.list')
    iraf.gemlist(range=range_string([star[j]]), root=[star[j]][0][:-4], 
        Stdout='star'+str(j)+'.list')
    iraf.gemlist(range=range_string(twilight[j]),
        root=twilight[j][0][:-4], Stdout='twilight'+str(j)+'.list')


#######################################################################
#######################################################################
###   Star reduction                                                  #
#######################################################################
#######################################################################

#
#   Flat reduction
#


for m in range(len(star)):
    for i in ['g', 'rg', 'erg']:
        iraf.imdelete(i+'@flat'+str(m)+'.list')

    # (B) - 'fl_over=no';
    iraf.gfreduce('@flat'+str(m)+'.list', slits='header', 
        rawpath='rawdir$', fl_inter='no',
        fl_addmdf='yes', key_mdf='MDF', mdffile='default', weights='no',
        fl_over='no', fl_trim='yes', fl_bias='yes', trace='yes', t_order=4,
        fl_flux='no', fl_gscrrej='no', fl_extract='yes', fl_gsappwave='no',
        fl_wavtran='no', fl_novl='no', fl_skysub='no', reference='',
        recenter='yes', fl_vardq='yes', bias='caldir$'+bias[m][0])


#
#   Twilight reduction
#
for m in range(len(star)):
    for i in ['g', 'rg', 'erg']:
        iraf.imdelete(i+'@twilight'+str(m)+'.list')

    # (B) - 'fl_over=no';
    iraf.gfreduce('@twilight'+str(m)+'.list', slits='header', rawpath='rawdir$',
        fl_inter='no', fl_addmdf='yes', key_mdf='MDF',
        mdffile='default', weights='no',
        fl_over='no', fl_trim='yes', fl_bias='yes', trace='yes', recenter='no',
        fl_flux='no', fl_gscrrej='no', fl_extract='yes', fl_gsappwave='no',
        fl_wavtran='no', fl_novl='no', fl_skysub='no',
        reference='erg'+flat[m][0], fl_vardq='yes', bias='caldir$'+bias[m][0])

#
#   Response function
#
for m in range(len(star)):
    for i, j in enumerate(flat[m]):
        iraf.imdelete('erg'+j+'_response')

        iraf.gfresponse('erg'+j+'.fits', out='erg'+j+'_response',
            skyimage='erg'+twilight[m][i], order=95, fl_inter='no', func='spline3',
            sample='*', verbose='yes')

#
#   Arc reduction
#
for m in range(len(star)):
    for i in arc[m]:
        iraf.imdelete('*'+i+'.fits')

    # (B) - 'fl_over=yes'; 'fl_bias=no'; #bias+,over- -> ERR
    # (B) - Eh preciso subtrair o bias para criar VAR/DQ
    # (B) - Da para usar bias-,over+  -> OK
    # (B) - Usarei 'fl_over=yes'; 'fl_bias=no', fl_vardq- -> OK
    iraf.gfreduce('@arc'+str(m)+'.list', slits='header', rawpath='rawdir$', fl_inter='no',
        fl_addmdf='yes', key_mdf='MDF', mdffile='default', weights='no',
        fl_over='yes', fl_trim='yes', fl_bias='no', trace='no', recenter='no',
        fl_flux='no', fl_gscrrej='no', fl_extract='yes', fl_gsappwave='no',
        fl_wavtran='no', fl_novl='no', fl_skysub='no', reference='erg'+flat[m][0],
        fl_vardq='no', bias='caldir$'+bias[m][0])

#
#   Finding wavelength solution
#   Note: the automatic identification is very good
#

for m in range(len(star)):
    for i in arc[m]:
        iraf.gswavelength('erg'+i, function='chebyshev', nsum=15, order=4,
            fl_inter='no', nlost=5, ntarget=20, aiddebug='s', threshold=5,
            section='middle line')

#
#   Apply wavelength solution to the lamp 2D spectra
#

for m in range(len(star)):
    for i in arc[m]:
        # (B) Fiz fl_vardq-, pois mudei acima tambem --------
        iraf.imdelete('terg'+i+'.fits')
        iraf.gftransform('erg'+i, wavtran='erg'+i, outpref='t', fl_vardq='no')

##
##   Actually reduce star
##

#-------------------------------------------------------------------------------------
for m, i in enumerate(star):
    iraf.imdelete('*g'+i+'.fits')

    # (B) - 'fl_over=no';
    iraf.gfreduce(
        i, slits='header', rawpath='rawdir$', fl_inter='no',
        fl_addmdf='yes', key_mdf='MDF', mdffile='default', weights='no',
        fl_over='no', fl_trim='yes', fl_bias='yes', trace='no', recenter='no',
        fl_flux='no', fl_gscrrej='no', fl_extract='no', fl_gsappwave='no',
        fl_wavtran='no', fl_novl='yes', fl_skysub='no', fl_vardq='yes',
        bias='caldir$'+bias[m][0])

    iraf.gemcrspec('rg{:s}.fits'.format(i), out='lrg'+i, sigfrac=0.32, 
         niter=4, fl_vardq='yes')

for m, i in enumerate(star):
    iraf.imdelete('*xlrg'+i+'.fits')

    iraf.gfreduce(
        'lrg'+i+'.fits', slits='header', rawpath='./', fl_inter='no',
        fl_addmdf='no', key_mdf='MDF', mdffile='default',
        fl_over='no', fl_trim='no', fl_bias='no', trace='no', recenter='no',
        fl_flux='no', fl_gscrrej='yes', fl_extract='yes', fl_gsappwave='yes',
        fl_wavtran='yes', fl_novl='no', fl_skysub='yes',
        reference='erg'+flat[m][0], weights='no',
        wavtraname='erg'+arc[m][0], response='erg'+flat[m][0]+'_response.fits',
        fl_vardq='yes')

# (( TEMP ))
#for m, i in enumerate(star):
#    # (B) Fiz fl_vardq-, pois mudei acima tambem --------
#    iraf.imdelete('texlrg'+i+'.fits')
#    iraf.gftransform('exlrg'+i, wavtran='erg'+arc[m][0], outpref='t', fl_vardq='yes')

    #
    #   Apsumming the stellar spectra
    #
for m, i in enumerate(star):
    iraf.imdelete('astexlrg'+i)
    iraf.gfapsum(
        'stexlrg'+i, fl_inter='no', lthreshold=400.,
        reject='avsigclip')
#-------------------------------------------------------------------------------------

#
#   Building sensibility function
#

for m in range(len(star)):
    iraf.delete('std'+str(m))
    iraf.delete('sens'+str(m)+'.fits')


for m in range(len(star)):
    iraf.gsstandard(
        'astexlrg'+star[m], starname=stdstar[m],
        observatory=observatory[m], sfile='std'+str(m),
        sfunction='sens'+str(m), caldir=caldir[m])


'''
# Create list of list of star (to be used in gsstandard)
# Nao funcionou

tmpstar = [star_idx, star]
star_idx0 = [[tmpstar[0][i] for i in range(len(star)) if ((star[i]==j))] 
            for j in list(set(star))]
star0 = [[tmpstar[1][i] for i in range(len(star)) if ((star[i]==j))] 
            for j in list(set(star))]

for m in range(len(star)):
    iraf.delete('std'+str(m)+'_'+str(i))
    iraf.delete('sens'+str(m)+'.fits')

for m in range(len(star0)):
    for i, n in enumerate(star0[m]):
        iraf.gsstandard(
            'astexlrg'+n, starname=stdstar[m],
            observatory=observatory[m], sfile='std'+str(m)+'_'+str(i),
            sfunction='sens'+str(m), caldir=caldir[m])
'''


# (B)Test calibration
import matplotlib.pyplot as plt

starflux = [locrep+starinfo[i]['flux'] for i in starinfo_idx]
for m, i in enumerate(star):
    iraf.delete('castexlrg'+i+'.fits')

    iraf.gscalibrate(
        'astexlrg'+i, sfuncti='sens'+str(m),
        observa=observatory[m], fluxsca=1)

    starcal_fl = pf.getdata('castexlrg'+i+'.fits', ext=2)
    starcal_he = pf.getheader('castexlrg'+i+'.fits', ext=2)
    starcal_wl = starcal_he['crval1'] + np.arange(starcal_he['naxis1'])*starcal_he['cdelt1']
    starrep = np.loadtxt(starflux[m], unpack=True)

    # Need: baixar arquivos de todas estrelas com contagem em fluxo
    #       ou converter o arquivo do gmos de mag>fluxo
    # Need: limitar range(x) do arquivo do repositorio
    # Need: limitar range(y) do plot. 
    plt.close('all')
    plt.plot(starcal_wl, starcal_fl, 'b')
    plt.plot(starrep[0], starrep[1], 'r', lw=1.56)
    plt.xlim(starcal_wl[0], starcal_wl[-1])
    plt.ylim(ymax=max(starrep[1])*1.3)
    plt.savefig('calib_'+starobj[m]+'.eps')


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

