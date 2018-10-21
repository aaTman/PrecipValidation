# -*- coding: utf-8 -*-
"""
Created on Sat Oct 13 23:50:23 2018

@author: Taylor
"""
class h5(hour=1):
    def __init__():
        self.filetype = hour
    def h5_gen(self, data, verif, wps, vps, latstr, lonstr, dates):
        '''
        Builds h5 files for the precentiles and precip data based on the length 
        of accumulation being generated.
        '''
        
        if self.hour == None:
            pass
        
        h5 = h5py.File(self.direct+'\\{}accum.h5'.format(str(self.hour)),'w')
        
        dset = h5.require_dataset('date',shape=dates.shape,dtype='f')
        dset.attrs['units'] = 'Hours Since December 31, 2001 23z'
        dset.attrs['name'] = 'Dates'
        dset[...]=dates
        
        dset = h5.require_dataset('lats',shape = .shape,dtype='f')
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