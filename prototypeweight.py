# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 18:01:48 2018

@author: Taylor
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from netCDF4 import Dataset
import pickle
import sys
import os
from datetime import datetime,timedelta
from pygeodesy.ellipsoidalVincenty import LatLon
import cartopy.crs as ccrs
import cartopy.feature as ft
import wrf
from scipy import interpolate
import pdb
#mainDir="C:\\Python\\verif"
#dirs=[x[0] for x in os.walk(mainDir)]
#pkldirs=[loc for loc in dirs if 'd03' in dirs]
#pklfiles=[n for n in os.listdir(loc) if '.pkl' in n for loc in pkldirs]
#n=0
#m=0
#with open(pkldirs[n]+'\\'+list(set(pklfiles))[m], 'rb') as f:
#    print(list(set(pklfiles))[m])
#    data = pickle.load(f)
#    print(list(data))
#    print(len(data))

lbot=964
rbot=1014
ltop=593
rtop=642
startdate=datetime(2001,12,31,23)
targetdate=datetime(2005,1,13,12)
extent=[-76, -72, 39, 42]
diff=targetdate-startdate
days, seconds = diff.days, diff.seconds

hours = (days * 24 + seconds // 3600)-27
hoursend=hours+48
mainDir="D:\\Verification\\nycDEP_wrf_5min"
dirs=[x[0] for x in os.walk(mainDir)]
ncdirs=[walk for walk in dirs if 'precip-5min' in walk]
ncfiles='prcp_d03.nc'
n=0
m=0
nct=Dataset(ncdirs[n]+'\\'+ncfiles)

# class verifWeight(self,dist=10):
#    def __init__(self,dist=10):
#        self.dist = dist

#    def latslons(self, lat, lon):
#        
#    def gridSize(self, x, y):
#        self.x = x
#        self.y = y
        
        
         
test1=Dataset('https://cida.usgs.gov/thredds/dodsC/stageiv_combined?Total_precipitation_surface_1_Hour_Accumulation[{}:1:{}][{}:1:{}][{}:1:{}],lon[{}:1:{}][{}:1:{}],lat[{}:1:{}][{}:1:{}],time[{}:1:{}]'.format(hours,hoursend,ltop,rtop,lbot,rbot,lbot,rbot,ltop,rtop,lbot,rbot,ltop,rtop,hours,hoursend))
lonsverif=test1['lon'][:]
latsverif=test1['lat'][:]
verifTimes = [startdate+timedelta(hours=n) for n in np.array(test1['time'][:])]

latmid=(int(latsverif.shape[0]/2),int(latsverif.shape[1]/2))
lonmid=(int(lonsverif.shape[0]/2),int(lonsverif.shape[1]/2))
ll = LatLon(latsverif[latmid], lonsverif[lonmid])
newlatlon=ll.destination(20000,45)
londiff=newlatlon.lon-lonsverif[lonmid]
latdiff=newlatlon.lat-latsverif[latmid]
precipdata=test1['Total_precipitation_surface_1_Hour_Accumulation'][:]
dir1="D:\\Verification\\"
latslons=Dataset(dir1+'geo_em.d03.nc')
lats=latslons['XLAT_M'][:]
lons=latslons['XLONG_M'][:]

times=[''.join([nct['Times'][m][n].decode("utf-8") for n in range(len(nct['Times'][m]))]) for m in range(len(nct['Times'][:]))]
strptimes = [datetime.strptime(n, '%Y-%m-%d_%H:%M:%S') for n in times]
timeArray=[np.where([strptimes[n].minute==0 for n in range(len(strptimes))])]

nctPrecip=nct['RAINNC'][timeArray[0][0]]
nctTimes= [strptimes[q] for x,q in enumerate(timeArray[0][0])]
lons=np.squeeze(lons)
lats=np.squeeze(lats)
states = ft.NaturalEarthFeature(category='cultural',scale='50m',facecolor='none',
                                            name='admin_1_states_provinces_lines')

for n in range(len(test1['time'])):
    grid_z0 = interpolate.griddata((lonsverif.flatten(),latsverif.flatten()),precipdata[n].T.flatten(), (lons.squeeze(),lats.squeeze()), method='nearest')
    for x in range(0,lons.shape[0]):
        for y in range(0,lons.shape[1]):
            
            indtestlon=np.where(np.logical_and(lons<=lons[y,x]+londiff,lons>=lons[y,x]-londiff))

            indcomblon=np.vstack((indtestlon[0][:],indtestlon[1][:])).T
    
            indtestlat=np.where(np.logical_and(lats<=lats[y,x]+latdiff,lats>=lats[y,x]-latdiff))
    
            indcomblat=np.vstack((indtestlat[0][:],indtestlat[1][:])).T
        
            pointsinrange=indcomblat[np.squeeze(np.array([np.where(np.prod(indcomblat==n, axis = -1)) for n in indcomblon if np.where(np.prod(indcomblat==n, axis = -1))[0].size>0]))]
    
            weightArray = np.ones_like(pointsinrange[:,0],dtype=np.float)
            setArray = pointsinrange-[y,x]

            testones=np.random.randint(0,high=20,size=lons.shape)
            testones[pointsinrange[:,0],pointsinrange[:,1]]=100

            for m in range(len(pointsinrange)):
                wt=float(np.max(abs(setArray[m])))
                if wt == 0:
                    pass
                else:
                    weightArray[m]=1./wt
    
            

            ax1 = plt.axes(projection=ccrs.PlateCarree())
            ax1.coastlines(resolution=('10m'), zorder=4)
            ax1.add_feature(ft.BORDERS,alpha=0.7,zorder=3)
            ax1.set_extent(extent,crs=ccrs.PlateCarree())
            ax1.add_feature(ft.RIVERS,alpha=0.4,edgecolor='gray', zorder=3)
            ax1.add_feature(states, edgecolor='gray', zorder=3)
            ax1.add_feature(ft.LAKES,alpha=0.5,facecolor='gray', zorder=3)
            c=ax1.pcolormesh(lons,lats,testones,transform=ccrs.PlateCarree(),cmap=plt.cm.Blues)
            #c1=ax1.pcolormesh(lons[0],lats[0],grid_z0,vmin=0,vmax=10,transform=ccrs.PlateCarree(),cmap=plt.cm.Blues)
            #plt.savefig('C:/Python/img/{}{}.png'.format(y,x),bbox_inches='tight')
            plt.show()
            c.remove()

                
