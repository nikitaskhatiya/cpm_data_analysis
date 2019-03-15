# Code for detecting and removing noise peaks from the CPM data. 
# Version 1.0
# Written by Dhruv Bal and Nikita Khatiya
# Date: 1 August 2018

#----------------------------------------------------------------------------------------------------------#

from __future__ import division
import os, warnings
import numpy as np
from astropy.io import fits,ascii
from astropy.table import Table
from matplotlib import pyplot as plt
from scipy.signal import find_peaks,peak_widths
import sys
import collections
import sqlite3

#----------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------------#

# Connect to the required databases
conn = sqlite3.connect('/home/czt/Desktop/Dhruv/New/All_data.db')
c = conn.cursor()
conn1 = sqlite3.connect('/home/czt/Desktop/Dhruv/New/All_data_l1.db')
c1 = conn1.cursor()

# Generate the table and declare the column names  
c1.execute('''CREATE TABLE data (Ind real, File_ID real, Time real,MJD real, Latitude real, Longitude real, CPM real,Veto_Q1 real,Veto_Q2 real,Veto_Q3 real,Veto_Q4 real)''')
conn1.commit()

# Read the File IDs of all the MKF files
tab=Table.read('/home/czt/Desktop/Dhruv/New/Fid.txt',format='ascii')
ext=np.array(tab[0][:])

#----------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------------#

for i in range(1,len(ext)):
	print ext[i]
	
	# Import the required data from the database 
	c.execute('''select Ind, File_ID,Time,MJD,Latitude,Longitude,CPM, Veto_Q1, Veto_Q2, Veto_Q3, Veto_Q4 from data where File_ID=(?)''',(ext[i],))
	ext1= np.array(c.fetchall())
	cpm=ext1[:,6]
	
	# Detect noise peaks that are present in the CPM data
	peaks, properties = find_peaks(cpm, prominence=250, width=(1,75))
	w=peak_widths(cpm,peaks,rel_height=0.5)
	for j in range (len(w[2])):
		
		# Replace the noise peaks with interpolated counts
		z=np.arange(int(math.floor(w[2][j]))-2,int(math.ceil(w[3][j]))+2,1)
		z=z.astype(int)
		cpm[z]=0
		cpm[cpm==0]='nan'
		ok=-np.isnan(cpm)
		xp=ok.ravel().nonzero()[0]
		fp=cpm[-np.isnan(cpm)]
		x=np.isnan(cpm).ravel().nonzero()[0]
		cpm[np.isnan(cpm)]=np.interp(x,xp,fp)
	ext1[:,6]=cpm
	
	# Write the data to a new database
	c1.executemany('''insert into data values (?,?,?,?,?,?,?,?,?,?,?)''', map(tuple, ext1.tolist()))
	conn1.commit()
c.close()
conn.close() 
c1.close()
conn1.close() 

#----------------------------------------------------------------------------------------------------------#	
