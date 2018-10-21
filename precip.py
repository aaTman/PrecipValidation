# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 15:49:17 2018

@author: Taylor
"""

import numpy as np
from netCDF4 import Dataset
import os
from datetime import datetime,timedelta
from pygeodesy.ellipsoidalVincenty import LatLon
from scipy import interpolate
from concurrent.futures import ProcessPoolExecutor
import xarray as xr

class precip_breakdown(object):
    
    def __init__(self, hour=1, directory=None):
        self.hour = hour
        if directory is not None:    
            self.directory = directory
        else:
            self.directory=self.get_path()
        self.latstr = 'XLAT'
        self.lonstr = 'XLONG'
        self.lats, self.lons = self.get_latlons()
        self.verif = self.load_verif()
        self.lbot=964
        self.rbot=1014
        self.ltop=593
        self.rtop=642
        self.veriflats=self.verif['lat'].T
        self.veriflons=self.verif['lon'].T
    def regrid(self):
 
        model_regrid = interpolate.griddata((self.lons.stack(z=('south_north', 'west_east')).reset_index('z'),self.lats.stack(z=('south_north', 'west_east')).reset_index('z')),self.get_hourly_precip().T.stack(z=('south_north', 'west_east')).reset_index('z'), (self.veriflons,self.veriflats), method='linear')
        return model_regrid
    
    def load_filelist(self):


        filenames=[n for n in os.listdir(self.directory) if os.path.isfile(os.path.join(self.directory, n))]
        return filenames   
    
    def load_file(self,ll=False):
                
            filenames=self.load_filelist()
            if ll==True:

                ll_file = Dataset(filenames[0])
                return ll_file
    def get_hourly_precip(self):
        
        bignc=self.load_file()
        model_precip_hourly = bignc['RAINNC'][:]
        return model_precip_hourly
    
    def load_verif(self):
        timedata=xr.open_dataset('https://cida.usgs.gov/thredds/dodsC/stageiv_combined?Total_precipitation_surface_1_Hour_Accumulation[0:1:146114][{}:1:{}][{}:1:{}],time[0:1:146114],lon[0:1:146114][{}:1:{}][{}:1:{}],lat[0:1:146114][{}:1:{}][{}:1:{}]'.format(self.ltop,self.rtop,self.lbot,self.rbot,self.lbot,self.rbot,self.ltop,self.rtop,self.lbot,self.rbot,self.ltop,self.rtop))
        return timedata

    def get_latlons(self):
        
        lats=self.load_file()[self.latstr][0]
        lons=self.load_file()[self.lonstr][0]
        
        return lats, lons
    
    def get_path(self):
        plotdirs=[]
        for root, dirs, files in os.walk(os.getcwd(), topdown=False):
           for name in dirs:
              plotdirs.append(os.path.join(root, name))
        return plotdirs
    


    def get_times(self,ncfile):
        
        ## Extract times
        times=[''.join([ncfile['Times'][m][n].decode("utf-8") for n in range(len(ncfile['Times'][m]))]) for m in range(len(nct['Times'][:]))]
    
        ## Convert times to hourly and datetime format
        strptimes = [datetime.strptime(n, '%Y-%m-%d_%H:%M:%S') for n in times]
        timeArray=[np.where([strptimes[n].minute==0 for n in range(len(strptimes))])]
        nctTimes= [strptimes[q] for x,q in enumerate(timeArray[0][0])]
    
    ## Take precip arrays with the time data
    modelPrecip=ncFILE['RAINNC'][timeArray[0][0]]
    
    ## Sets the arrays from total accum to hourly precip
    modelPrecip=np.array([modelPrecip[n]-modelPrecip[n-1] for n in range(len(modelPrecip)) if n>0])
    
    ## Generates hourly, 6 hourly, or 24 hourly data depending on accumu set earlier in script
    modelPrecip=np.array([np.sum(modelPrecip[n-accumu:n],axis=0) for n in range(accumu,len(nctTimes))])
        
def run_loop():
    '''
    Runs through the directory and collects data from each file
    This will probably need to be batchified...

    '''
            