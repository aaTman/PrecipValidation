# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 18:01:48 2018

@author: Taylor
"""

import numpy as np
from netCDF4 import Dataset
import os
from datetime import datetime,timedelta
from pygeodesy.ellipsoidalVincenty import LatLon
from scipy import interpolate
import h5py

def h5gen(data,verif,wps,vps,lats,lons,dates,accum=None):
    if accum == None:
        pass
    
    h5 = h5py.File('D:\\Verification\\{}accum.h5'.format(str(accum)),'w')
    
    dset = h5.require_dataset('date',shape=dates.shape,dtype='f')
    dset.attrs['units'] = 'Hours Since December 31, 2001 23z'
    dset.attrs['name'] = 'Dates'
    dset[...]=dates
    
    dset = h5.require_dataset('lats',shape = lats.shape,dtype='f')
    dset.attrs['name']='latitude'
    dset.attrs['units']='degrees north'
    dset[...]=lats
    dset = h5.require_dataset('lons',shape=lons.shape,dtype='f')
    dset[...]=lons
    dset.attrs['name']='longitude'
    dset.attrs['units']='degrees east'
    
    dset=h5.require_dataset("verif",shape = verif.shape,dtype='f')
    dset.attrs['units'] = 'mm'
    dset.attrs['name'] = 'Stage IV Linearly Interpolated 1 km Accumulated Precipitation'
    dset[...]=verif
    
    dset=h5.require_dataset("wrf",shape = data.shape,dtype='f')
    dset.attrs['units'] = 'mm'
    dset.attrs['name'] = 'WRF 1 km Accumulated Precipitation'
    dset[...]=data
    
    dset=h5.require_dataset("wrfpercentile",shape = data.shape,dtype='f')
    dset.attrs['units'] = 'percent'
    dset.attrs['name'] = 'WRF 1 km Accumulated Precipitation'
    dset[...]=wps
    
    dset=h5.require_dataset("verifpercentile",shape = data.shape,dtype='f')
    dset.attrs['units'] = 'percent'
    dset.attrs['name'] = 'WRF 1 km Accumulated Precipitation'
    dset[...]=vps
    
    h5.close()

## Accessing local dirs
plotdirs=[]
dir1="D:\\Verification\\"
latslons=Dataset(dir1+'geo_em.d03.nc')
lats=latslons['XLAT_M'][:]
lons=latslons['XLONG_M'][:]
mainDir="D:\\Verification\\nycDEP_wrf_5min"
for root, dirs, files in os.walk("D:\Verification\Plots", topdown=False):
   for name in dirs:
      plotdirs.append(os.path.join(root, name))
    
dirs=[x[0] for x in os.walk(mainDir)]
wrfdirs=[walk for walk in dirs if 'precip-5min' in walk]
wrffiles='prcp_d03.nc'

## Grid for Stage IV closest to WRF
lbot=964
rbot=1014
ltop=593
rtop=642


        
## Hour analysis
accumu=24
## To place into arrays for entire dataset analysis

wrfDist = np.zeros((10,(49-accumu),148,148))
verifDist = np.zeros((10,(49-accumu),148,148))
hourTimes = np.zeros((10,(49-accumu)))
c=0
d=0
## Walks through the directories to gather the data

for dirs in wrfdirs:
    
    ## Assigns nc wrf file to variable
    nct=Dataset(dirs+'\\'+wrffiles)
    ## Extract times
    times=[''.join([nct['Times'][m][n].decode("utf-8") for n in range(len(nct['Times'][m]))]) for m in range(len(nct['Times'][:]))]
    
    ## Convert times to hourly and datetime format
    strptimes = [datetime.strptime(n, '%Y-%m-%d_%H:%M:%S') for n in times]
    timeArray=[np.where([strptimes[n].minute==0 for n in range(len(strptimes))])]
    nctTimes= [strptimes[q] for x,q in enumerate(timeArray[0][0])]
    
    ## Take precip arrays with the time data
    modelPrecip=nct['RAINNC'][timeArray[0][0]]
    
    ## Sets the arrays from total accum to hourly precip
    modelPrecip=np.array([modelPrecip[n]-modelPrecip[n-1] for n in range(len(modelPrecip)) if n>0])
    
    ## Generates hourly, 6 hourly, or 24 hourly data depending on accumu set earlier in script
    modelPrecip=np.array([np.sum(modelPrecip[n-accumu:n],axis=0) for n in range(accumu,len(nctTimes))])
    
    
    ## Starting date for Stage IV record
    startdate=datetime(2001,12,31,23)
    
    ## Datetime which Stage IV data should start at
    targetdate=nctTimes[0]
    
    ## Find hour to start the Stage IV thredds pull from and the end time
    diff=targetdate-startdate
    days, seconds = diff.days, diff.seconds
    hours = (days * 24 + seconds // 3600)
    hoursend=hours+48
    
    ## Correction for Stage IV thredds issue
    timedata=Dataset('https://cida.usgs.gov/thredds/dodsC/stageiv_combined?time[0:1:146114]')
    hours=np.where(timedata['time'][:]==hours)[0][0]
    hoursend=np.where(timedata['time'][:]==hoursend)[0][0]
    
    ## Pull Stage IV data
    test1=Dataset('https://cida.usgs.gov/thredds/dodsC/stageiv_combined?Total_precipitation_surface_1_Hour_Accumulation[{}:1:{}][{}:1:{}][{}:1:{}],lon[{}:1:{}][{}:1:{}],lat[{}:1:{}][{}:1:{}],time[{}:1:{}]'.format(hours,hoursend,ltop,rtop,lbot,rbot,lbot,rbot,ltop,rtop,lbot,rbot,ltop,rtop,hours,hoursend))
    lonsverif=test1['lon'][:]
    latsverif=test1['lat'][:]
    
    ## Sanity check to make sure times are correct
    verifTimes = [startdate+timedelta(hours=n) for n in np.array(test1['time'][:])]
    print(verifTimes[0],nctTimes[0])
    

    
    ## Stage IV precip data
    verifdata=test1['Total_precipitation_surface_1_Hour_Accumulation'][:]
    n1=48
    hourTimes[d] = np.linspace(hours,hoursend,49)[accumu-1:n1]
    verifnhr=np.array([np.sum(np.array(verifdata[n-accumu:n]),axis=0) for n in range(accumu,len(nctTimes))])
    
    ## Reshapes the Stage IV grid to WRF grid shape linearly (least noise in regridded histogram compared to cubic and nearest neighbor)
    sivRegrid = np.array([np.array(interpolate.griddata((lonsverif.flatten(),latsverif.flatten()),n.T.flatten(), (lons.squeeze(),lats.squeeze()), method='linear')) for n in verifdata])
    
    ## Same as previous for WRF, calculates accumulated precip array for Stage IV data
    sivRegrid=np.array([np.sum(np.array(sivRegrid[n-accumu:n]),axis=0) for n in range(accumu,len(nctTimes))])

    ## Builds array for all events
    wrfDist[d] = modelPrecip
    verifDist[d] = sivRegrid

    ## Print max of each array for each event
    print(np.max(modelPrecip))
    print(np.max(sivRegrid))
    
    ## Increments for the big arrays to store all cases
    c+=len(modelPrecip)
    d+=1
    
percwrf = wrfDist.flatten()
percverif=verifDist.flatten()

wrfsort=np.argsort(percwrf)
verifsort=np.argsort(percverif)

wrfsortsort=np.argsort(wrfsort)
verifsortsort=np.argsort(verifsort)

wrfpercent=(percwrf[wrfsort]/np.max(percwrf))*100
verifpercent=(percverif[verifsort]/np.max(percverif))*100

wps=wrfpercent[wrfsortsort]
vps=verifpercent[verifsortsort]
wps=wps.reshape(wrfDist.shape)
vps=vps.reshape(verifDist.shape)
h5gen(wrfDist,verifDist,wps,vps,lats,lons,hourTimes,accum=accumu)

## Everything below this is my sandbox, nothing important yet!
    
## This is for smoothing/weighting function which might be removed
#    latmid=(int(latsverif.shape[0]/2),int(latsverif.shape[1]/2))
#    lonmid=(int(lonsverif.shape[0]/2),int(lonsverif.shape[1]/2))
#    ll = LatLon(latsverif[latmid], lonsverif[lonmid])
#    newlatlon=ll.destination(25000,45)
#    londiff=newlatlon.lon-lonsverif[lonmid]
#    latdiff=newlatlon.lat-latsverif[latmid]
    
# class verifWeight(self,dist=10):
#    def __init__(self,dist=10):
#        self.dist = dist

#    def latslons(self, lat, lon):
#        
#    def gridSize(self, x, y):
#        self.x = x
#        self.y = y
        
    
#    extent=[np.min(lons),np.max(lons),np.min(lats),np.max(lats)]
#    states = ft.NaturalEarthFeature(category='cultural',scale='50m',facecolor='none',
#                                                name='admin_1_states_provinces_lines')
#    for q in range(len(sivRegrid)):
#        
#        fig=plt.figure(figsize=(5,5))
#        ax1 = plt.axes([0.1,0.7,0.8,0.7],projection=ccrs.PlateCarree())
#        ax2 = plt.axes([0.1,0.69,1.6,0.05])
#        ax3 = plt.axes([.9,0.7,0.8,0.7],projection=ccrs.PlateCarree())
#
#        ax1.coastlines(resolution=('10m'), zorder=4)
#        ax1.add_feature(ft.BORDERS,alpha=0.7,zorder=3)
#        ax1.set_extent(extent,crs=ccrs.PlateCarree())
#        ax1.add_feature(ft.RIVERS,alpha=0.4,edgecolor='gray', zorder=3)
#        ax1.add_feature(states, edgecolor='gray', zorder=3)
#        ax1.add_feature(ft.LAKES,alpha=0.5,facecolor='gray', zorder=3)
#        c1=ax1.pcolormesh(lons[0],lats[0],verif6Dist[c+q],vmin=0,vmax=np.max((np.max(sivRegrid),np.max(modelPrecip))),transform=ccrs.PlateCarree(),cmap=plt.cm.Blues)
#        
#        ax3.coastlines(resolution=('10m'), zorder=4)
#        ax3.add_feature(ft.BORDERS,alpha=0.7,zorder=3)
#        ax3.set_extent(extent,crs=ccrs.PlateCarree())
#        ax3.add_feature(ft.RIVERS,alpha=0.4,edgecolor='gray', zorder=3)
#        ax3.add_feature(states, edgecolor='gray', zorder=3)
#        ax3.add_feature(ft.LAKES,alpha=0.5,facecolor='gray', zorder=3)
#        ax3.pcolormesh(lons[0],lats[0],wrf6Dist[c+q],vmin=0,vmax=np.max((np.max(sivRegrid),np.max(modelPrecip))),transform=ccrs.PlateCarree(),cmap=plt.cm.Blues)
#        
#        cbar=plt.colorbar(c1,orientation='horizontal',pad=0.02,aspect=30,cax=ax2)
#        cbar.set_label('6-Hourly Precip Accumulation')
#        #c1=ax1.pcolormesh(lons[0],lats[0],sivRegrid,vmin=0,vmax=10,transform=ccrs.PlateCarree(),cmap=plt.cm.Blues)
#        #plt.savefig('C:/Python/img/{}{}.png'.format(y,x),bbox_inches='tight')
#        plt.show()
#        plt.close(fig) 

#    lons=np.squeeze(lons)
#    lats=np.squeeze(lats)

#    exceedGrid=np.array(np.zeros_like(modelPrecip))
#    validArray=np.zeros((49))
#for n in range(len(test1['time'])):
#    distArray=np.zeros((148*148))
#    wrfPrecip=np.array(modelPrecip[n])
#    sivRegrid = np.array(interpolate.griddata((lonsverif.flatten(),latsverif.flatten()),verifdata[n].T.flatten(), (lons.squeeze(),lats.squeeze()), method='linear'))
#    diffGrid=wrfPrecip-sivRegrid
#
#    cc=0
#    for x in range(0,lons.shape[0]):
#        for y in range(0,lons.shape[1]):
#            
#            indtestlon=np.where(np.logical_and(lons<=lons[y,x]+londiff,lons>=lons[y,x]-londiff))
#
#            indcomblon=np.vstack((indtestlon[0][:],indtestlon[1][:])).T
#    
#            indtestlat=np.where(np.logical_and(lats<=lats[y,x]+latdiff,lats>=lats[y,x]-latdiff))
#    
#            indcomblat=np.vstack((indtestlat[0][:],indtestlat[1][:])).T
#            
#            pointsinrange=indcomblat[np.squeeze(np.array([np.where(np.prod(indcomblat==n, axis = -1)) for n in indcomblon if np.where(np.prod(indcomblat==n, axis = -1))[0].size>0]))]
#
#            verifGrid=diffGrid[pointsinrange[:,0],pointsinrange[:,1]]
#            weightArray = np.ones_like(verifGrid,dtype=np.float)
#            zeroThreshArray = np.zeros_like(verifGrid,dtype=np.float)
#            setArray = pointsinrange-[y,x]
#            
#            #testones[pointsinrange[:,0],pointsinrange[:,1]]=100
#            
#            for m in range(len(verifGrid)):
#                wt=float(np.max(abs(setArray[m])))
#                
#                #zeroThreshArray[m]=1
#                if wt == 0:
#                    pass
#                else:
#                    weightArray[m]=1./(wt+0.1)
#            stackArray=np.vstack((np.expand_dims(weightArray,axis=0),setArray.T)).T
#            resx=int(np.max(stackArray[:,1])-np.min(stackArray[:,1])+1)
#            intArray = np.array((stackArray[:,1]-np.min(stackArray[:,1]),stackArray[:,2]-np.min(stackArray[:,2]))).T.astype(int)
#            resy=int(np.max(stackArray[:,2])-np.min(stackArray[:,2])+1)
##            
##            try:
#            resArray=np.zeros((resx,resy))
#            d=0
#            c=0
#            for xx,yy in intArray:
#                
#                if sivRegrid[xx,yy]>0.01:
#                    
#                    zeroThreshArray[c]=1
#                resArray[xx,yy]=stackArray[c,0]
#                c+=1
#                    
#            exceedGrid[n,x,y] = np.sum(zeroThreshArray)/len(zeroThreshArray)
#
#                #resArray=np.reshape(stackArray[:,0],(resx,resy))
##            except IndexError:
##                pdb.set_trace()
#            cc+=1
#

#    
#            #verifGrid*=weightArray 
#            #distArray[cc] = np.sum(abs(verifGrid))
#    #validArray[n] = np.mean(distArray)
#            #diffGrid[pointsinrange[:,0],pointsinrange[:,1]]=diffGrid[pointsinrange[:,0],pointsinrange[:,1]] * weightArray
#            
#            

                
