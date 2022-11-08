# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 09:12:14 2022

@author: Riccardo
"""

import pandas as pd

data4digits = pd.read_csv('NACA4digits.csv')
data4digits = data4digits.drop(['Cl', 'Cm'], axis=1)

data4digits['ID'] = data4digits['ID'].str.replace('NACA', '')

# make the new columns using string indexing
data4digits['max_cam'] = data4digits['ID'].str[0:1]
data4digits['max_pos'] = data4digits['ID'].str[1:2]
data4digits['thick']   = data4digits['ID'].str[2:4]

# get rid of the extra variable (if you want)
data4digits.drop('ID', axis=1, inplace=True)

# removing leading zeros from thick column
data4digits['thick'] = data4digits['thick'].str.lstrip('0')

data4digits.to_csv('NACA4digits_preprocessed.csv')