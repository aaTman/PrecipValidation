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

lbot=970
rbot=1010
ltop=595
rtop=640
startdate=datetime(2001,12,31,23)
targetdate=datetime(2005,1,13,13)
extent=[-76, -72, 39, 42]
diff=targetdate-startdate
days, seconds = diff.days, diff.seconds

hours = days * 24 + seconds // 3600
hoursend=hours+48
mainDir="C:\\Python\\wrf_precip\\nycDEP_wrf_5min"
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
lons=test1['lon'][:]
lats=test1['lat'][:]

latmid=(int(lats.shape[0]/2),int(lats.shape[1]/2))
lonmid=(int(lons.shape[0]/2),int(lons.shape[1]/2))
ll = LatLon(lats[latmid], lons[lonmid])
newlatlon=ll.destination(20000,45)
londiff=newlatlon.lon-lons[lonmid]
latdiff=newlatlon.lat-lats[latmid]
precipdata=test1['Total_precipitation_surface_1_Hour_Accumulation'][:]

for x in range(0,lons.shape[0]):
    for y in range(0,lons.shape[1]):
#        a=datetime.now()
        indtestlon=np.where(np.logical_and(lons<=lons[x,y]+londiff,lons>=lons[x,y]-londiff))
#        print('indtestlon {}.{:1.0f}'.format((datetime.now()-a).seconds, (datetime.now()-a).microseconds/1000))
#        a=datetime.now()
        indcomblon=np.vstack((indtestlon[0][:],indtestlon[1][:])).T

        indtestlat=np.where(np.logical_and(lats<=lats[x,y]+latdiff,lats>=lats[x,y]-latdiff))

        indcomblat=np.vstack((indtestlat[0][:],indtestlat[1][:])).T
    
        np.squeeze(np.array([np.where(np.prod(indcomblat==n, axis = -1)) for n in indcomblon if np.where(np.prod(indcomblat==n, axis = -1))[0].size>0]))
       
        pointsinrange=indcomblat[np.squeeze(np.array([np.where(np.prod(indcomblat==n, axis = -1)) for n in indcomblon if np.where(np.prod(indcomblat==n, axis = -1))[0].size>0]))]


        #test1['Total_precipitation_surface_1_Hour_Accumulation'][0,pointsinrange[:,1],pointsinrange[:,0]]
        
#        
#        states = ft.NaturalEarthFeature(category='cultural',scale='50m',facecolor='none',
#                                            name='admin_1_states_provinces_lines')
#        fig=plt.figure(figsize=(10,10))
#        testones=np.random.randint(0,high=20,size=lons.shape)
#        testones[pointsinrange[:,0],pointsinrange[:,1]]=100
#        print(x,y)
#        print(pointsinrange[:,0]-x,pointsinrange[:,1]-y)
        weightArray = np.ones_like(pointsinrange[:,0],dtype=np.float)
        setArray = pointsinrange-[x,y]
        for n in range(len(pointsinrange)):
            wt=float(np.max(abs(setArray[n])))
            if wt == 0:
                pass
            else:
                weightArray[n]=1./wt

#        print(pointsinrange)
#        print(setArray)
#        print(weightArray)
        
#        testones[pointsinrange[:,0],pointsinrange[:,1]]*=weightArray
#        ax1 = plt.subplot(2, 1, 1, projection=ccrs.PlateCarree())
#        ax1.coastlines(resolution=('10m'), zorder=4)
#        ax1.add_feature(ft.BORDERS,alpha=0.7,zorder=3)
#        ax1.set_extent(extent,crs=ccrs.PlateCarree())
#        ax1.add_feature(ft.RIVERS,alpha=0.4,edgecolor='gray', zorder=3)
#        ax1.add_feature(states, edgecolor='gray', zorder=3)
#        ax1.add_feature(ft.LAKES,alpha=0.5,facecolor='gray', zorder=3)
#        c=ax1.pcolormesh(lons[:],lats[:],testones,transform=ccrs.PlateCarree(),cmap=plt.cm.Blues)
#        plt.colorbar(c)
        #plt.savefig('C:/Python/img/{}{}.png'.format(x,y),bbox_inches='tight')
#        plt.show()
#        plt.close(fig)        
#        states = ft.NaturalEarthFeature(category='cultural',scale='50m',facecolor='none',
#                                            name='admin_1_states_provinces_lines')
#        fig=plt.figure(figsize=(10,10))
#        testones=np.random.randint(0,high=20,size=lons.shape)
#        ax1 = plt.subplot(2, 1, 1, projection=ccrs.PlateCarree())
#        ax1.coastlines(resolution=('10m'), zorder=4)
#        ax1.add_feature(ft.BORDERS,alpha=0.7,zorder=3)
#        ax1.set_extent([-76, -72, 39, 42],crs=ccrs.PlateCarree())
#        ax1.add_feature(ft.RIVERS,alpha=0.4,edgecolor='gray', zorder=3)
#        ax1.add_feature(states, edgecolor='gray', zorder=3)
#        ax1.add_feature(ft.LAKES,alpha=0.5,facecolor='gray', zorder=3)
#        c=ax1.pcolormesh(lons[:],lats[:],test1['Total_precipitation_surface_1_Hour_Accumulation'][0].T,transform=ccrs.PlateCarree(),cmap=plt.cm.Blues)
#        plt.colorbar(c)
#        plt.show()
#        plt.close(fig)