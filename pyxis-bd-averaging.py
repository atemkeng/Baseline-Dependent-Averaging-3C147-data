"""
Module import here
"""
import os
import sys
import MSResampler
import Pyxis
from Pyxis.ModSupport import *
import mqt
import pyrap.tables
from pyrap.tables import table
import numpy as np
import ms
import imager
import scipy.special
import pylab
import pyfits
import Tigger

"""
     Conversion modules
"""
Radian2Arcmin = lambda x: x * 180 * 60 /np.pi;
Radian2deg = lambda x: x * 180./np.pi
Radian2min = lambda x:x*12*60/np.pi;
#v.DESTDIR_Template = '${OUTDIR>/}plots-${MS:BASE}${-stage<STAGE}'
def givephasecenter (msname = None):

	"""Gets  MS Dec and RA.
                Return couple: (Dec, RA) in (radian,radian)
        """;
	msname = msname or v.MS
        tab = ms.ms(msname, "FIELD").getcol("PHASE_DIR");
        dec = Radian2Arcmin (tab[0,0,1]);
        ra = Radian2min(tab[0,0,0]);
        return tab[0,0,1], tab[0,0,0];

def changePhaseCenter (msname=None, LSM=None, SUBLSM=None, DESTDIR=None):

	""" Change the phase centre of the sky model
	to feed the one of your MS
	""";
	DESTDIR = DESTDIR or "."
	msname = msname or v.MS
	dec0, ra0 = givephasecenter(msname);
	dec_d = math.floor(Radian2deg (dec0)); 
	dec_m = Radian2deg (dec0) - dec_d 
	ra_d = math.floor(Radian2deg (ra0))
	ra_m = Radian2deg (ra0) - ra_d 
	SUBLSM = SUBLSM or "%s-ra%dh%dmin0s-dec%ddeg%dmin0s.lsm.html"%(LSM.split(".")[0],ra_d,ra_m,dec_d,dec_m)
	info("******* %s phase center at ra=%dh%dmin0s, dec=%ddeg%darmin0arcsec *******"%(SUBLSM,ra_d,ra_m,dec_d,dec_m))
	command = "tigger-convert -f %s %s --recenter=%s,%dh%dm0s,%dd%dm0s"%(LSM,SUBLSM,"j2000",ra_d,ra_m,dec_d,dec_m)
	os.system(command)
	v.SUBLSM = II("${DESTDIR}/$SUBLSM")
	if os.path.exists(DESTDIR) is not None:
		makedir("$DESTDIR")
	if DESTDIR != "." :
		if os.path.exists(v.SUBLSM) is not None:
			x.mv("$SUBLSM $DESTDIR");
		else:
			x.sh("rm -fr ${DESTDIR}/$SUBLSM");	
			x.mv("$SUBLSM $DESTDIR");
	return v.SUBLSM;
    

def antennasDiameter(msname = None):
	"""
	this function gives the diameter in meter of the interferometer antenna
	""";
	msname = msname or v.MS;
        tab = table(msname+"/ANTENNA")
	tab = tab.getcol("DISH_DIAMETER");
	#info("***** Antenna size %f meter*******"%tab[0])
	print "",tab[0]
	stop
        return tab[0];

def makems (name=None, destdir = None, Ntime=8,
	                                        integration=60,dec=-45,ra=0,startfq=1400,nchan=1,chanwidth=10,starttime=None):
	"""Makes an MS using the specified parameters.
	Ntime: total synthesis duration
	integration: timeslot in seconds
	dec: declination, in degrees
	freq0: starting freq in Hz
	nchan: number of channels
	chanwidth: channel width, in Hz
	"""
	destdir = destdir or "TEMPLE"
	name,destdir = interpolate_locals("name destdir");
	# for anttab in conf,"${conf}_ANTENNA","Layouts/${conf}_ANTENNA":
	# 	if exists(anttab):
	# 		anttab = II(anttab);
	# 		break;
	# else:
	# 	abort("configuration $conf not found");
	# if not name:
	# 	name = os.path.basename(anttab);
	# 	if name.endswith("_ANTENNAS"):
	# 		name = name.rsplit("_",1)[0];
	msname = "%s-%.1fs-%.1fMHz.MS"%(name,integration,chanwidth*1e-6);
	info("ms $msname");
	if exists(msname):
		x.sh("rm -fr $msname");
	conffile = II("makems.${msname:BASE}.cfg");
	file(conffile,"w").write(II("""
WriteAutoCorr=Tlo
StartFreq=%g
StepFreq=%g
NFrequencies=$nchan
WriteImagerColumns=T
StepTime=$integration
#TileSizeRest=1 0
NParts=1
MSDesPath=.
AntennaTableName=ANTENNA
Declination=$dec
NBands=1
RightAscension=$ra
StartTime=$starttime
MSName=$msname
NTimes=$Ntime 
 #TileSizeFreq=16
"""%(startfq,chanwidth)));
	info("""creating $msname: ${Ntime} timslots, ${integration}s integration, Dec=$dec, Ra=$ra
                                                          $nchan channels of $chanwidth Hz starting at $startfreq Hz, Nband 2""");
        # run makems
	x.makems(conffile);
	if exists(msname+"_p0") and not exists(msname):
		x.mv("${msname}_p0 $msname");
	if os.path.exists(destdir) is not None:
		makedir("$destdir")
	try:
		x.mv ("$msname $destdir")
	except:
                x.sh("rm -fr $destdir/$msname");
		x.mv ("$msname $destdir")
		info("overight and saved $msname")
	v.MS = II("${destdir}/$msname");
	return v.MS

def mscriticalsampling (msname=None, DISHDIA = None):
	"""
	Work out MS critical sampling time of the Longest baseline;
	This should be the low resolution MS maximun integration
	""";
	msname = msname or v.MS;
	DISHDIA = DISHDIA or antennasDiameter(msname)
	tab = ms.ms(msname);
	times = sorted(set(tab.getcol("TIME")));
	dt = times[1]-times[0];
	tottime = times[-1]-times[0]+dt;
	uvw = tab.getcol("UVW");
	maxbl = math.sqrt((uvw[:,:2]**2).sum(1).max());
	# max baseline sweeps out a circle of length pi*d over 24h
	arclen = math.pi*maxbl*(tottime/(24.*3600));
	# critical sampling is at half the dish size
	nsamp = arclen/(DISHDIA/2)
	# corresponding sampling interval
	critint = tottime/nsamp;
	info("(NB: critical sampling interval for this baseline is %.3fs, diam=%f, *****PLEASE USED THIS INTEGRATION FOR THE LOW RESOLUTION MS)"%(critint,DISHDIA))
	return critint;

def msinfo (msname=None):
	"""Gets various MS info.
	Return values: NUM_TIMESLOTS,NUM_CHANNELS,
	MAX_BASELINE,TOTAL_TIME,TIMESLOT_SIZE,FIRST_TIMESTAMP 
	""";
	msname = msname or v.MS;
	tab = ms.ms(msname);
	chanwidths = ms.ms(msname,"SPECTRAL_WINDOW").getcol("CHAN_WIDTH",0,1);
	nchan = chanwidths.size;
	times = sorted(set(tab.getcol("TIME")));
	dt = times[1]-times[0];
	tottime = times[-1]-times[0]+dt;
	ntime = len(times);
	uvw = tab.getcol("UVW");
	maxbl = math.sqrt((uvw[:,:2]**2).sum(1).max());
	return ntime,nchan,maxbl,tottime,dt,chanwidths[0,0];


def simulate_imaging_bd_3c147 (hiresms=None,inputcolumn="DATA",outputcolumn="CORRECTED_DATA",time0=0,\
				 dtime=None,dfreq=None,freq0=0,nfreq=None,ntime=None):
	"""
		The function make use of baseline dependent averaging on the 3C147 real data observed on the
		2013/01/27/01:00:28.5, at declination 49.51.07.23356, ra 05:42:36.137916
		NB: We have removed the fictive baseline of index 7
	"""	

	# make an instance of the class containing the method for bd-averaging
	mshi = MSResampler.MSResampler(hiresms+"/",column=inputcolumn,time0=time0,ntime=ntime,freq0=freq0,nfreq=nfreq)
	# BD-averaging, giving the integration time of the shortest baseline, dtime and 
	# the number of uv-frequency bins dfreq to average
	arrays = mshi.bd_averaging (dtime,dfreq)
	# arrays is of size (p,q,dtimepq,data)
	# take the number of time bins of the longest baseline and make low res timeslots
	p,q,dtlong,arr = arrays[200] ;
	Ntime = arr.shape[1]#3238
	# take the integration time of the longest baseline and make low res integration time
	integration = dtlong;
	# parameters for low res
	dec = "49.51.07.23356"
	ra =  "05:42:36.137916"
	nchan = arr.shape[2]*3;
	chanwidth = dfreq*1000000;
	starttime = "2013/01/27/01:00:28.5"
	antenna = "ANTENNA"
	startfq = 1329e6;
	namel = "3C147-lores"
	
	# make the measurment set corresponding to the low res 
        lowresms = makems (name=namel,Ntime=Ntime,integration=integration,\
      		dec=dec,ra=ra,startfq=startfq,nchan=nchan,chanwidth=chanwidth,starttime=starttime);
	# # save visibilities
    	MSResampler.save_visibility_arrays (lowresms,arrays,column=outputcolumn)
  	imager.npix= 2048
	imager.cellsize = "2arcsec"
	imager.stokes   = "I"
	imager.weight   = "natural"
	imager.wprojplanes = 128
	#cleaning options
	imager.niter = 10000000
	#imager.threshold = "5mJy"
	imager.CLEAN_ALGORITHM = "csclean"
	
	imager.make_image(msname = lowresms, column = 'CORRECTED_DATA',restore = True, dirty = False, weight = "natural");

