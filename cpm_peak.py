# Code for detecting the CPM peaks and recording them according to the path traced by ASTROSAT. 
# Version 1.1
# Written by Dhruv Bal and Nikita Khatiya
# Date: 21 August 2018

#----------------------------------------------------------------------------------------------------------#

from __future__ import division
import os, warnings
import numpy as np
from astropy.io import fits,ascii
from astropy.table import Table
from astropy.modeling import models,fitting
from matplotlib import pyplot as plt
from scipy.signal import find_peaks,peak_widths
from itertools import groupby
from operator import itemgetter
import sys
import collections
import sqlite3
from lmfit.models import GaussianModel

#----------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------------#

# Connect to the required databases 
conn = sqlite3.connect('/home/czt/Desktop/Dhruv/New/All_data_l1.db')
c = conn.cursor()
conn1 = sqlite3.connect('/home/czt/Desktop/Dhruv/New/Asc.db')
c1 = conn1.cursor()
conn2 = sqlite3.connect('/home/czt/Desktop/Dhruv/New/Dec.db')
c2 = conn2.cursor()

# Generate the tables with the required columns
c1.execute('''CREATE TABLE asc (Ind real, File_ID real ,MJD real, Latitude real, Longitude real, CPM real)''')
conn1.commit()
c2.execute('''CREATE TABLE dec (Ind real, File_ID real ,MJD real, Latitude real, Longitude real, CPM real)''')
conn2.commit()

tab=Table.read('/home/czt/Desktop/Dhruv/New/Fid.txt',format='ascii')
ext=np.array(tab[0][:])

#----------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------------#

for i in range(1,len(ext)):
	
	# Extracting the data from the database with the help of a query
	c.execute('''select Ind, File_ID,MJD,Latitude,Longitude,CPM from data where File_ID=(?)''',(ext[i],))
	ext1= np.array(c.fetchall())
	ind=ext1[:,0]
	fid=ext1[:,1]
	mjd=ext1[:,2]
	lat=ext1[:,3]
	lon=ext1[:,4]
	cpm=ext1[:,5]
	
	# Detect the CPM peaks 
	peaks, properties = find_peaks(cpm, prominence=30, width=10) # Record the CPM peak value indices
	w=peak_widths(cpm,peaks,rel_height=0.8)
	pd=np.diff(peaks)
	pd1=np.where(pd<1000)[0]
	if len(pd1)>0:
		rp=[]
		res1=[]
		for k,l in groupby (enumerate(pd1),lambda(i,x):i-x):
			i1=map(itemgetter(1),l)
			i2=np.append(i1,i1[-1]+1)
			ix=peaks[i2]
			fix=np.where(cpm[ix]==max(cpm[ix]))[0]
			res1=np.append(res1,ix[fix])
			rp=np.append(rp,ix)
		res1=res1.astype(int)
		rp=rp.astype(int)
		res2=np.setxor1d(peaks, rp)
		fp=np.sort(np.append(res1,res2))
	else:
		fp=peaks
	fp=fp.astype(int)
	
	# Extract the part of orbit's data where it traces the SAA region   
	for z in range (len(fp)):
		if fp[z]>=650:
			ind1=ind[fp[z]-650:fp[z]+650]
			mjd1=mjd[fp[z]-650:fp[z]+650]
			fid1=fid[fp[z]-650:fp[z]+650]
			data1=cpm[fp[z]-650:fp[z]+650]
			lat1=lat[fp[z]-650:fp[z]+650]
			lon1=lon[fp[z]-650:fp[z]+650]
		elif fp[z]+650<=len(cpm)-650:
			ind1=ind[fp[z]-650:len(cpm)-1]
			mjd1=mjd[fp[z]-650:len(cpm)-1]
			fid1=fid[fp[z]-650:len(cpm)-1]
			data1=cpm[fp[z]-650:len(cpm)-1]
			lat1=lat[fp[z]-650:len(cpm)-1]
			lon1=lon[fp[z]-650:len(cpm)-1]
		else:
			ind1=ind[0:fp[z]+650]
			mjd1=mjd[0:fp[z]+650]
			fid1=fid[0:fp[z]+650]
			data1=cpm[0:fp[z]+650]
			lat1=lat[0:fp[z]+650]
			lon1=lon[0:fp[z]+650]

		# Fit the CPM data and find the peak value and it's corresponding coordinates in latitude and longitude
		if len(data1)>0:
			if max(data1)>50:
				mod = GaussianModel() 
				pars = mod.guess(data1, x=lon1)
				out = mod.fit(data1, pars, x=lon1) # Fit the CPM data with a gaussian curve
				if max(out.best_fit)>50:
					med=np.where(out.best_fit==max(out.best_fit))[0] # Find the index of the peak count
					sub= np.diff(lat1)
					maj1=np.where(sub>0)[0]
					maj2=np.where(sub<0)[0]					
					ind2=ind1[med]
					fid2=fid1[med]
					mjd2=mjd1[med]
					lat2=lat1[med]
					lon2=lon1[med]
					fit2=out.best_fit[med]
					
					# Check if the orbit traces an Ascending or descending path based on its latitude values					
					if len(maj1)>(len(lat1)/2+1):
						dat=np.transpose(np.array([ind2,fid2,mjd2,lat2,lon2,fit2]))
						c1.executemany('''insert into asc values (?,?,?,?,?,?)''', map(tuple, dat.tolist()))
						conn1.commit()
					if len(maj2)>(len(lat1)/2+1):
						dat=np.transpose(np.array([ind2,fid2,mjd2,lat2,lon2,fit2]))
						c2.executemany('''insert into dec values (?,?,?,?,?,?)''', map(tuple, dat.tolist()))
						conn2.commit()
c.close()
conn.close() 
c1.close()
conn1.close() 
c2.close()
conn2.close() 
	
#----------------------------------------------------------------------------------------------------------#
