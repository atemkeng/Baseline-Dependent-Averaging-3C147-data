import numpy as np
from pyrap.tables import table
import os
import sys
import pylab

def save_visibility_arrays (msname,arrays,column="CORRECTED_DATA"):
  """Saves a set of visibility arrays to the MS.
  arrays is a list of (p,q,data) tuples, where p & q are antenna indices,
  and data is the corresponding array of shape ntime,nfreq,ncorr.
  The shape must match the shape of the MS, or errors will result. 
  """
  # add the flagging column to the MS structure
  os.system("addbitflagcol %s"%msname)
  # Load the MS
  tab = table(msname,readonly=False,ack=False);
  # read in data column and antenna indices, fill data with 0
  data = tab.getcol(column)
  a1 = tab.getcol("ANTENNA1")
  a2 = tab.getcol("ANTENNA2")
  data.fill(0.);
  # read spectral windows ID
  data_desc_id = tab.getcol("DATA_DESC_ID");
  bitflagrow = tab.getcol("BITFLAG_ROW");
  bitflagrow.fill(0);
  bitflag = tab.getcol("BITFLAG");
  bitflag.fill(0);
  ##tab.close();
  # get the integration for the longest baseline and the shortest we know the position in the arrays
  p,q,dtlong,arr1 = arrays[100] ;
  p,q,dtshort,arr2 = arrays[51];
  # number of timeslots, you can still take this from arr1 by doing this arr1.shape[0]
  ntimeslots = data[(a1==0)&(a2==1)].shape[0];
  #ntimeslots = arr1.shape[1]
  # number of data point of the longest baseline if dtshort is given for the shortest baseline
  # We need this to know how we save the data into MS, we have to avoid phase shift when flagging
  # notice that for one integration dtshort, we have to flag nbdatalongbase*constant data points
  nbdatalongbase = int(np.ceil(dtshort/dtlong));
  # keep the indexes to flag	
  flagarray = [];
  for p,q,dtime,arr in arrays:
     # concatenate the data of all the bands to one band
     # Notice that MS is of a unique band
     dataav =  np.concatenate((arr[0],arr[1],arr[2]),axis=1)
     # work out the data to be flag for this baseline
     nbdatapq = int(np.ceil(dtshort/dtime));
     # give us the step between flag data and observing data
     step_index = (nbdatalongbase-nbdatapq)/2;
     datapq = data[(a1==p)&(a2==q)].copy()
     bitflagrowpq = bitflagrow[(a1==p)&(a2==q)].copy()
     bitflagpq = bitflag[(a1==p)&(a2==q)].copy()
     # Initialize variable that define the starting area to flag
     in_time0 = 0;
     # Initialize variable that define the ending area to flag
     in_time1 = in_time0+step_index;
     # for idw in range(arr.shape[0]):
     # while in_timesw0 <= dataav.shape[0]-nbdatapq:
     # iteration for the average data
     for in_timesw0 in range(0,dataav.shape[0],nbdatapq):
	 print "==== Baseline Dependent Averaging: Baseline (%d,%d), save the data ===="%(p,q)
         # make timeslice corresponding to the flagging area
	 timesliceflag = slice(in_time0,in_time1);
	 # flag the data in the range timesliceflag
	 bitflagrowpq[timesliceflag] = 1;
	 bitflagpq[timesliceflag] = 1;
	 print "",bitflagpq[timesliceflag]
	 ##flagarray.append((p,q,in_time0,in_time1));
	 #os.system("flag-ms.py -S %i,%i -T %i:%i -c %s"%(
	 #				p,q,in_time0,(in_time0+step_index),msname));
	
	 # make timeslice corresponding to the range to save data
         timeslicedata = slice(in_time1+1,in_time1+nbdatapq);
	 # this contained the range to extract the average data
	 timeslicesw = slice(in_timesw0,in_timesw0+nbdatapq);
	 # dataav.shape[0] can be shorter that ntimeslots if not longest baselines
	 try :
                datapq[timeslicedata] = dataav[timeslicesw];
         except:
		# in the case where data[(a1==p)&(a2==q),...][timeslicedata].shape[0] is longer that dataav[timeslicesw].shape[0]
		# make data[(a1==p)&(a2==q),...][timeslicedata].shape[0] equal to dataav[timeslicesw].shape[0]
		new_index = len(dataav[timeslicesw]);
                timeslicedata = slice(in_time1,in_time1+new_index)
		# save the troncatenate data
		print "in_time1,in_time1+new_index",in_time1,in_time1+new_index
                datapq[timeslicedata] = dataav[timeslicesw];    
	 # change the indexes
	 in_time0 = in_time1 +1 + nbdatapq;
	 n_time1 = in_time0+step_index;
	 # flag the second area
	 timesliceflag = slice(in_time0,in_time1);
	 #flagarray.append((p,q,in_time0,in_time1))
	 bitflagrowpq[timesliceflag] = 1;
	 bitflagpq[timesliceflag] = 1;
	 #flagrow[(a1==p)&(a2==q),...][timesliceflag] = 1;
	 # change indexes
	 in_time0 = in_time1 +1;#+ step_index;
	 in_time1 = in_time0 + step_index;
     timesliceflag = slice(in_time0,ntimeslots);
     bitflagrowpq[timesliceflag] = 1;
     bitflagpq[timesliceflag] = 1;
     data[(a1==p)&(a2==q)] = datapq.copy()
     bitflagrow[(a1==p)&(a2==q)] = bitflagrowpq.copy()
     bitflag[(a1==p)&(a2==q)] = bitflagpq.copy()
  # save changes
  tab.putcol(column,data)
  tab.putcol("BITFLAG_ROW",bitflagrow)
  tab.putcol("BITFLAG",bitflag)
  print "save succeded"
  tab.close();


class MSResampler (object):
    """Class for reading and resampling data from an MS"""

    def __init__ (self,msname,column="DATA",time0=0,ntime=None,freq0=0,nfreq=None):
      """Inits with given MS and column name.
      If time0/ntime and freq0/nfreq is given, handles only a subset of the data,
      otherwise reads all channels and timeslots.
      """;
      self.msname = msname;
      tab = table(msname,ack=False,readonly=False)
      self.A0 = A0 = tab.getcol('ANTENNA1')
      self.A1 = A1 = tab.getcol('ANTENNA2')
      data = tab.getcol(column);
      #spectral window id 
      data_desc_id = tab.getcol("DATA_DESC_ID")
      if nfreq is None:
        nfreq = data.shape[1]-freq0;
      #self.data = data = data[:,freq0:freq0+nfreq,:];
      self.data = data.copy();
      self.data_desc_id = data_desc_id.copy()
      self.nfreq = nfreq;
      self.freq0=freq0
      print "Visibility column shape:",data.shape
      self.na = na = np.max(A1)+1
      # do not consider auto-correlation
      self.nbl = (na*(na-1))/2 
      self.ncorr = data.shape[2]
      self.nbins = data.shape[0]
      self.UVW = tab.getcol("UVW")
      # get number of timeslots. This assumes a regular MS staructure (i.e. first baseline is
      # present)
      #ntimes = (data[(A0==0)&(A1==1)]).shape[0]
      # actual number of timeslot to use given by time0 and ntime
      self.ntime = ntime #or ntimes - time0;
      self.time0 = time0;

      # get frequency and wavelength (per channel) from SPECTRAL_WINDOW subtable
      t2 = table(tab.getkeyword("SPECTRAL_WINDOW"),readonly=False)
      self.channels = t2.getcol("CHAN_FREQ",0);
      self.freqs = t2.getcol("CHAN_FREQ",0)[0,freq0:freq0+nfreq];
      # print "Frequencies are",self.freqs;
      self.wavel = 3e8/self.freqs;
      # print "Wavelengths are",self.wavel;
      t2.close()
      tab.close()
   
    def bd_averaging (self,dtime,dfreq):
      """Downsamples data using baseline dependent averaging. The compression is done over a boxcar window of size dfreqxdtimepq
      Returns list of (antenna1,antenna2,dtimepq,averaged_data) tuples, one per each baseline, where
      averaged_data has the shape (ntime/dtimepq,nfreq/dfreq,ncorr)."""
      # this is the frequency  output shape
      # Notice that the time output depends on the baseline length
      nfreq1 = self.nfreq/dfreq;
      # make a list of per-baseline output arrays 
      result = [];
      # shortest baseline
      pqs = [2,8]; 
      # longest baseline
      pql = [4,17];
      # shortest baseline length
      lengthsh = np.sqrt((self.UVW[(self.A0==2)&(self.A1==8)]**2).sum(axis=1))[0];
      # loop over each baseline
      # Evaluate the integration time for all baseline, given the one of the shortest baseline
      for p in range(self.na):
        for q in range(p+1,self.na):
          # extract full array for this baseline, apply subset in time
          input_index = (self.A0==p)&(self.A1==q)
	  # extract the data and the  spectral windows indexes for this baseline
          data = self.data[input_index].copy();
          data_desc_id = self.data_desc_id[input_index].copy()
	  # uvw bin for this baseline
	  uvw=self.UVW[(self.A0==p)&(self.A1==q)]
          # the baselines with 0 data are fictives in VLA this observations. Do not consisder 
	  # One option is to remove these fictives Antennas from the Antenna table
	  if len(data)!=0:
                # evaluate the integration time of this baseline, NB : dtime is the integration time for the shortest baseline
		# length for this baseline
		lengthpq = np.sqrt((uvw**2).sum(axis=1))[0];
                # compresion time or integration time
		dtimepq = np.ceil(dtime/(lengthpq/lengthsh));
		print "==== Baseline Dependent Averaging: Baseline (%d,%d), integration =%ds this may take some time ===="%(p,q,dtimepq)
		# number of spectral windows
		nbsw = data_desc_id.max()+1;
		# this is the ouput shape in time
		ntime1 = int((data.shape[0]/nbsw)/dtimepq);
		# work out the total number of bins to avg for this baseline
		ntime = int(ntime1*dtimepq);
		# prepare list for result of each spectral window
		resultsw = []
		# extract the data of each spectract window, NB: number of spectral window = data_desc_id.max()+1;
		for idw in range(nbsw):
			# extract data for this spectral window
			dataid = data[data_desc_id==idw]; 
	  		dataid = dataid[self.time0:self.time0+ntime,self.freq0:self.freq0+self.nfreq,:];
          		# reshape so that e.g. ntime becomes ntime1,dtime, so that we can then reduce over the second axis
          		dataid = dataid.reshape((ntime1,dtimepq,nfreq1,dfreq,self.ncorr))
			# keep this spectral windows
                        # take mean over the new axes, this becomes our result
			resultsw.append(dataid.mean(3).mean(1))
          	# save the result
          	result.append((p,q,dtimepq,np.array(resultsw)));
      return result;
