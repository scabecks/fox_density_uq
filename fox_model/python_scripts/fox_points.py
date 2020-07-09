#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 21:17:42 2020

@author: scottca
"""
import os
import pandas as pd
import rasterio as rio
from time import strptime

def precip_6_month(date):
    idx = r_list.index('wc2.1_2.5m_prec_' + date + '.tif')
    r0 = rio.open('worldclim/' + r_list[idx])#.read(1)
    profile = r0.profile
    r0 = r0.read(1)
    #r0[np.isnan(r0)] = 0
    r1 = rio.open('worldclim/' + r_list[idx - 1]).read(1)
    #r2[np.isnan(r2)] = 0
    r2 = rio.open('worldclim/' + r_list[idx - 2]).read(1)
    #r2[np.isnan(r2)] = 0
    r3 = rio.open('worldclim/' + r_list[idx - 3]).read(1)
    #r3[np.isnan(r3)] = 0
    r4 = rio.open('worldclim/' + r_list[idx - 4]).read(1)
    #r4[np.isnan(r4)] = 0
    r5 = rio.open('worldclim/' + r_list[idx - 5]).read(1)
    #r5[np.isnan(r5)] = 0
    r6 = rio.open('worldclim/' + r_list[idx - 6]).read(1)
    #r6[np.isnan(r6)] = 0
    rs = r0 + r1 + r2 + r3 + r4 + r5 +r6
    with rio.open('temp.tif','w', **profile) as dst:
        dst.write(rs, 1)
        
         
#Read in foxes data
foxes = pd.read_csv('fox_densities_march2018-enso.csv')
foxes = foxes[foxes['month'].notna()]
foxes['month_int'] = foxes['month'].apply(lambda x: f'{strptime(x, "%b").tm_mon:02d}') 
foxes['r_index'] = foxes.apply(lambda x: x['year_of_data_collection'] + '-' + str(x['month_int']), axis = 1)
foxes = foxes[foxes['lon'].notna()]
foxes['worldclim'] = float(0)
foxes.index = range(len(foxes))
foxes['coords'] = [(x,y) for x, y in zip(foxes.lon, foxes.lat)]

r_list = os.listdir("worldclim")
r_list.sort()


for index, row in foxes.iterrows():
    coords = row['coords']
    precip_6_month(row['r_index'])
    src = rio.open('temp.tif')
    for val in src.sample([(coords)]):
        val[0]
    foxes.at[index, 'worldclim'] = val
    os.remove('temp.tif')

foxes.to_csv('fox_densities_march2018-enso_worldclim.csv')
