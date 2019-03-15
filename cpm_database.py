# Code for generating a database of all the CPM data. 
# Version 1.0
# Written by Dhruv Bal and Nikita Khatiya
# Date: 25th July 2018

#----------------------------------------------------------------------------------------------------------#

from __future__ import division
import os, warnings
import numpy as np
from astropy.io import fits,ascii
from matplotlib import pyplot as plt
import sys
import collections
import sqlite3

#----------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------------#

# Reading the names of all the MKF files available
f = open('filename.txt', 'r')
x =np.array( f.readlines())
f.close()
x=np.array(map(lambda z: z[:-1],x))
temp=0

#----------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------------#

# Creating a databse and connecting to it
conn=sqlite3.connect('/home/czt/Desktop/Dhruv/New/All_data.db')
c = conn.cursor()

# Declaring a table with the columns required in it
c.execute('''CREATE TABLE data (Ind real,File_ID real, Time real,MJD real, Latitude real, Longitude real, CPM real,Veto_Q1 real,Veto_Q2 real,Veto_Q3 real,Veto_Q4 real)''')
conn.commit()

#----------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------------#

namespace=globals()
for i in range(len(x)):
	
	# Importing the required parameters from the MKF Files
	namespace['hdu%d'%i]=fits.open(x[i])
	mjdrefi=namespace['hdu%d' %i][0].header['MJDREFI']
	mjdreff=namespace['hdu%d' %i][0].header['MJDREFF']
	namespace['t%d'%i]=np.transpose(np.array(namespace['hdu%d'%i][1].data['TIME']))
	namespace['lat%d'%i]=np.transpose(np.array(namespace['hdu%d'%i][1].data['EARTHLAT']))
	namespace['lon%d'%i]=np.transpose(np.array(namespace['hdu%d'%i][1].data['EARTHLON']))
	namespace['cpm%d'%i]=np.transpose(np.array(namespace['hdu%d'%i][1].data['CPM_Rate']))
	namespace['vet1q%d'%i]=np.transpose(np.array(namespace['hdu%d'%i][1].data['Q1_VetoCounter']))
	namespace['vet2q%d'%i]=np.transpose(np.array(namespace['hdu%d'%i][1].data['Q2_VetoCounter']))
	namespace['vet3q%d'%i]=np.transpose(np.array(namespace['hdu%d'%i][1].data['Q3_VetoCounter']))
	namespace['vet4q%d'%i]=np.transpose(np.array(namespace['hdu%d'%i][1].data['Q4_VetoCounter']))
	namespace['hdu%d'%i].close()
	
	# Processing the data according to the format required by the database file 
	t=namespace['t%d'%i]
	offset=mjdrefi+mjdreff
	mjd=(namespace['t%d'%i]/86400)+offset # Generating the MJD values for corresponding time indices
	fil_id=x[i][47:52]
	fid=np.transpose(np.repeat(fil_id,len(namespace['t%d'%i])))
	temp1=temp+len(namespace['t%d'%i])
	ind=np.transpose(np.arange(temp,temp1,1))
	temp=temp1
	ind=np.reshape(ind,(len(ind),1))
	t=np.reshape(namespace['t%d'%i],(len(t),1))
	mjd=np.reshape(mjd,(len(t),1))
	fid=np.reshape(fid,(len(t),1))
	lat=np.reshape(namespace['lat%d'%i],(len(t),1))
	lon=np.reshape(namespace['lon%d'%i],(len(t),1))
	cpm=np.reshape(namespace['cpm%d'%i],(len(t),1))
	vet1=np.reshape(namespace['vet1q%d'%i],(len(t),1))
	vet2=np.reshape(namespace['vet2q%d'%i],(len(t),1))
	vet3=np.reshape(namespace['vet3q%d'%i],(len(t),1))
	vet4=np.reshape(namespace['vet4q%d'%i],(len(t),1))
	
	# Creating a matrix of all the data required to be written to the database
	mat=np.concatenate((ind,fid,t,mjd,lat,lon,cpm,vet1,vet2,vet3,vet4),axis=1)
	
	# Writing to the database
	c.executemany('''insert into data values (?,?,?,?,?,?,?,?,?,?,?)''', map(tuple, mat.tolist()))
	conn.commit()
	
c.close()
conn.close() 

#----------------------------------------------------------------------------------------------------------#
