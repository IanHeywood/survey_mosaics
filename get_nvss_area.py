#----------------------------------------------------------------------
#
# Define a sky area using min/max RA/Dec and this script will produce 
# an NVSS image that covers the area.
#
# It invokes both wget and some of the Montage tools from the command 
# line, so make sure they are installed and in your path.
# 
# It also requires the nvss_fields.p pickle.
#
# ian.heywood@csiro.au
# 11 April 2017
#
#----------------------------------------------------------------------

import pickle
import numpy
import os 
import glob
import string
import random
import pyfits
from multiprocessing import Pool
from astLib import astCoords as ac
from optparse import OptionParser

parser = OptionParser()
parser.add_option('--x0',dest='ra_min',help='Minimum Right Ascension in degrees',default=160.0)
parser.add_option('--x1',dest='ra_max',help='Maximum Right Ascension in degrees',default=165.0)
parser.add_option('--y0',dest='dec_min',help='Minimum Declination in degrees',default=-9.0)
parser.add_option('--y1',dest='dec_max',help='Maximum Declination in degrees',default=1.0)
parser.add_option('-j','--j',dest='ncpu',help='Number of parallel reprojections to invoke',default=8)
parser.add_option('-o','--out',dest='out',help='Name of final FITS image',default='my_nvss_mosaic.fits')
parser.add_option('-c','--cleanup',dest='cleanup',action='store_true',help='Remove temporary files (default = False)',default=False)

(options,args) = parser.parse_args()
ra_min = float(options.ra_min)
ra_max = float(options.ra_max)
dec_min = float(options.dec_min)
dec_max = float(options.dec_max)
ncpu = int(options.ncpu)
out = options.out
cleanup = options.cleanup

fields = pickle.load(open('nvss_fields.p','rb'))


#----------------------------------------------------------------------
# FUNCTION DEFINITIONS
#----------------------------------------------------------------------


def get_nvss(id):
	syscall = 'wget ftp://nvss.cv.nrao.edu/pub/nvss/MAPS/'+id+'.gz'
	os.system(syscall)


def in_area(ra,dec,ra_min,ra_max,dec_min,dec_max):
	inside = False
	if ra > ra_min*0.9 and ra < ra_max*1.1:
		if dec > dec_min*0.9 and dec < dec_max*1.1:
			inside = True
	return inside


def tempname(size=12,chars=string.ascii_uppercase+string.digits+string.ascii_lowercase):
        return ''.join(random.choice(chars) for _ in range(size))


def reprofits(inpimg):
        outputimg = 'repro/'+inpimg.rstrip('.fits')+'.repro.fits'
        cmd = 'mProject '+inpimg+' '+outputimg+' template.hdr'
        print inpimg,'--->',outputimg
        os.system(cmd)


def fixMontageHeaders(infile,outfile,axes):
	inphdu = pyfits.open(infile)
	inphdr = inphdu[0].header
	outhdu = pyfits.open(outfile,mode='update')
	outhdr = outhdu[0].header
	keywords = ['CTYPE','CRVAL','CDELT','CRPIX']
	for axis in axes: 
		for key in keywords:
			inkey = key+str(axis)
			outkey = key+str(axis)
			afterkey = key+str(axis-1)
			xx = inphdr[inkey]
			outhdr.set(outkey,xx,after=afterkey)
	outhdr.set('BUNIT',inphdr['BUNIT'],after=outkey)
	outhdr.set('BMAJ',inphdr['BMAJ'],after='BUNIT')
	outhdr.set('BMIN',inphdr['BMIN'],after='BMAJ')
	outhdr.set('BPA',inphdr['BPA'],after='BMIN')
	outhdu.flush()


#----------------------------------------------------------------------


tempdir = tempname()
os.system('mkdir -p '+tempdir+'/repro')

# ra_min_h = ac.decimal2hms(ra_min,delimiter=':')
# ra_max_h = ac.decimal2hms(ra_max,delimiter=':')

downloads = []

for field in fields:
	f_ra_h = field[1:3]
	f_ra_m = field[3:5]
	f_dec_sign_str = field[5]
	if f_dec_sign_str == 'M':
		f_dec_sign = '-'
	else:
		f_dec_sign = '' 
	f_dec_d = field[6:8]
	f_ra_hms = f_ra_h+':'+f_ra_m+':00'
	f_dec_dms = f_dec_sign+f_dec_d+':00:00'
	f_ra = ac.hms2decimal(f_ra_hms,delimiter=':')
	f_dec = ac.dms2decimal(f_dec_dms,delimiter=':')
	if in_area(f_ra,f_dec,ra_min,ra_max,dec_min,dec_max):
		downloads.append(field)

for dl in downloads:
	get_nvss(dl)
	os.system('gunzip '+dl+'.gz')
	os.system('mv '+dl+' '+dl+'.fits')
	os.system('mSubimage -p '+dl+'.fits '+dl+'_crop.fits 0 0 1024 1024')
	os.system('rm '+dl+'.fits')
	os.system('mv '+dl+'_crop.fits '+tempdir)

os.chdir(tempdir)
os.system('mImgtbl . images.tbl')
os.system('mMakeHdr images.tbl template.hdr')

fitslist = glob.glob('*_crop.fits')
proclist = []
for infits in fitslist:
	opfits = 'repro/'+infits.replace('.fits','_repro.fits')
	proclist.append(infits)
pool = Pool(processes=ncpu)
pool.map(reprofits,proclist)

os.chdir('repro')
os.system('mImgtbl . images.tbl')
os.system('mAdd images.tbl ../template.hdr '+out)
os.system('mv '+out+' ../../')
os.chdir('../../')
fixMontageHeaders(tempdir+'/'+downloads[0]+'_crop.fits',out)
if cleanup:
	os.system('rm -rf '+tempdir)
