# Code for clubbing the common orbits i.e. orbits which nearly trace the same path. 
# Version 1.0
# Written by Dhruv Bal and Nikita Khatiya
# Date: 18 November 2018

#----------------------------------------------------------------------------------------------------------#

import sqlite3
from astropy.table import Table,Column,vstack
import numpy as np
from matplotlib import pyplot as plt
import itertools
from itertools import groupby
from operator import itemgetter
import operator

#----------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------------#

# Connect to the required database and import data from text files 
tab1=Table.read('/home/czt/Desktop/Dhruv/New/Median_A.txt',format='ascii')
lata=np.array(tab1[0][:])
lona=np.array(tab1[1][:])
conn = sqlite3.connect('/home/czt/Desktop/Dhruv/New/All_data_l1.db')
c = conn.cursor()
conn1=sqlite3.connect('/home/czt/Desktop/Dhruv/New/asc_test1.db')
c1=conn1.cursor()

#----------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------------#

# 
for i in range(len(lata)):
	
	# Import the data, based on the median value it crosses, from the database 
	c.execute('''select Ind,File_ID from data where Latitude >=(?) AND Latitude <(?) and Longitude >=(?) AND Longitude <(?) ''',(lata[i],lata[i+1],lona[i]-0.05,lona[i+1]+0.05))
	ext= np.array(c.fetchall())
	
	if len(ext)>0:
		nam=str('asc%d'%i)
		print nam
		
		# Create a table in the new database and declare the column headers
		c1.execute('''CREATE TABLE %s ( Ind real,File_ID real,MJD real, Time real, Latitude real, Longitude real, CPM real)'''%(nam),)
		conn1.commit() 
		ind=ext[:,0].astype(int)
		fid=ext[:,1].astype(int)
		
		# Find the indices of unique orbits for further analysis
		pd=np.diff(ind)
		pd1=np.where(np.abs(pd)>3000)[0]
		if len(pd1)>0:
			fp=np.append(ind[pd1],ind[pd1[-1]+1])
			fd=np.append(fid[pd1],fid[pd1[-1]+1])
			for j in range(len(fp)):
				
				# Extract the entire orbit containing the full Gaussian distribution of CPM data
				l1=fp[j]-2000
				l2=fp[j]+1950
				c.execute('''select Ind,File_ID,MJD,Time,Latitude,Longitude,CPM from data where Ind >=(?) AND Ind<=(?) AND File_ID=(?)''',(l1,l2,fd[j]))
				dat=np.array(c.fetchall())
				mlat=dat[:,4]
				mlon=dat[:,5]
				mlat=mlat.astype(float)
				mlon=mlon.astype(float)
				indl=np.where((mlon>=-95)&(mlon<=10))[0]
				indl=indl.astype(int)
				dlat=mlat[indl]
				dif=np.diff(dlat)
				maj=np.where(dif>0)[0]
				
				# Check if the extracted orbit is either an Ascending or Descending orbit								
				if len(maj)>(len(dlat)/2+1) and len(mlat)>2000:
					
					# Write the data to the new database 
					c1.executemany('''insert into %s values (?,?,?,?,?,?,?)'''%(nam),map(tuple, dat.tolist()))
					conn1.commit()
c1.close()
conn1.close()

#----------------------------------------------------------------------------------------------------------#
