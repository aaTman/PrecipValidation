# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 21:30:22 2018

@author: Taylor
"""
import numpy as np
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import h5py
from netCDF4 import Dataset
import os
import sys

class modelData:
    
    mainDir = os.path.dirname(os.path.realpath(__file__))
    
    def __init__(self):
        pass
    
    def unPickle(self):
        with open(mainDir+, 'rb') as f:
            data = pickle.load(f)
            
    