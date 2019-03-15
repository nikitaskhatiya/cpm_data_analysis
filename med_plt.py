# Code for detecting the CPM peaks and recording them according to the path traced by ASTROSAT. 
# Version 1.1
# Written by Dhruv Bal and Nikita Khatiya
# Date: 21 August 2018

#----------------------------------------------------------------------------------------------------------#

from __future__ import division
import os, warnings
import numpy as np
from astropy.modeling import models,fitting
from astropy.io import fits,ascii
from astropy.table import Table,Column
from itertools import groupby
from operator import itemgetter
from scipy import interpolate as itp
from scipy.interpolate import UnivariateSpline
from scipy.signal import find_peaks,peak_widths
from scipy.optimize import curve_fit
import scipy.stats
import sys
import collections
import sqlite3
from lmfit.models import GaussianModel
from sklearn.cluster import KMeans
from matplotlib import pyplot as plt

#----------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------------#

# Connect to the database
conn1 = sqlite3.connect('/home/czt/Desktop/Dhruv/New/A.db')
c1 = conn1.cursor()
conn2 = sqlite3.connect('/home/czt/Desktop/Dhruv/New/D.db')
c2 = conn2.cursor()

output_file_path1='/home/czt/Desktop/Dhruv/New/Median_A.txt'
output_file_path2='/home/czt/Desktop/Dhruv/New/Median_D.txt'

#----------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------------#

# Importing the CPM peak data from the database
c1.execute('''select Ind, File_ID,MJD,Latitude,Longitude,CPM from asc where Longitude>=-60 and Longitude<=-40''')
ext1= np.array(c1.fetchall())
c2.execute('''select Ind, File_ID,MJD,Latitude,Longitude,CPM from dec where Longitude>=-70 and Longitude<=-30''')
ext2= np.array(c2.fetchall())

# Storing the latitude and longitude values of peaks in separate variables
lon1A=ext1[:,4]
lat1A=ext1[:,3]
lon1D=ext2[:,4]
lat1D=ext2[:,3]
lon1A=np.reshape(lon1A,(np.shape(lon1A)[0],1))
lat1A=np.reshape(lat1A,(np.shape(lat1A)[0],1))
A21=np.where(lat1A>-4)[0]
nlatA1=lat1A[A21]
nlonA1=lon1A[A21]
X1=np.concatenate((nlonA1, nlatA1), axis=1)
A22=np.where(lat1A<=-4)[0]
nlatA2=lat1A[A22]
nlonA2=lon1A[A22]
X2=np.concatenate((nlonA2, nlatA2), axis=1)

lon1D=np.reshape(lon1D,(np.shape(lon1D)[0],1))
lat1D=np.reshape(lat1D,(np.shape(lat1D)[0],1))
X3=np.concatenate((lon1D, lat1D), axis=1)

#----------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------------#

# Defining the number of clusters required by the K_Mean algorithm
Kmean1 = KMeans(n_clusters=4)
Kmean2 = KMeans(n_clusters=1)
Kmean3 = KMeans(n_clusters=1)

# Using the cluster centers to fit the data
Kmean1.fit(X1)
Kmean2.fit(X2)
Kmean3.fit(X3)

# Storing the points of the fit for better accuracy
new_row1= np.array([[-53,-3.97656],[-40.7419,4.25]]) 
new_array1=np.vstack([Kmean1.cluster_centers_, new_row1])
kc1=new_array1[new_array1[:, 0].argsort()]
kx1 = kc1[:,0]
ky1 = kc1[:,1]

new_row2= np.array([[-53.04,-3.83656],[-50.4395,-6.03906]])
new_array2=np.vstack([Kmean2.cluster_centers_, new_row2])
kc2=new_array2[new_array2[:, 0].argsort()]
kx2 = kc2[:,0]
ky2 = kc2[:,1]


new_row3= np.array([[-48.3065,-6.14062],[-32.8226,4.28646]])
new_array3=np.vstack([Kmean3.cluster_centers_, new_row3])
kc3=new_array3[new_array3[:, 0].argsort()]
kx3 = kc3[:,0]
ky3 = kc3[:,1]

# Fitting the data and storing the points of the fit in variables
mytck1,myu1=itp.splprep([kx1,ky1])
xnew1,ynew1= itp.splev(np.linspace(0,1,400),mytck1)
i1=np.where(ynew1>-3.83691)[0]
mytck2,myu2=itp.splprep([kx2,ky2],k=2)
xnew2,ynew2=itp.splev(np.linspace(0,1,200),mytck2)
fint=np.append(ynew2,ynew1[i1])
fint=np.sort(fint)[::-1]
finn=np.append(xnew2,xnew1[i1])
finn=np.sort(finn)[::-1]
fint=np.reshape(fint,(np.shape(fint)[0],1))
finn=np.reshape(finn,(np.shape(finn)[0],1))
fX3=np.concatenate((fint, finn), axis=1)
fX3=fX3[fX3[:, 0].argsort()]
mytck3,myu3=itp.splprep([kx3,ky3],k=1)
xnew3,ynew3=itp.splev(np.linspace(0,1,550),mytck3)

#----------------------------------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------------------------------#

# Exporting the points of the fit to a text file
tab1=Table()
tab2=Table()
tab1.add_column(Column(data=fX3[:, 0],name='Latitude'))
tab1.add_column(Column(data=fX3[:, 1],name='Longitude'))	
tab2.add_column(Column(data=ynew3,name='Latitude'))
tab2.add_column(Column(data=xnew3,name='Longitude'))	

sys.stdout = open(output_file_path1, 'w')
tab1.pprint(max_lines=-1, max_width=-1)
sys.stdout = open(output_file_path2, 'w')
tab2.pprint(max_lines=-1, max_width=-1)

#----------------------------------------------------------------------------------------------------------#
