# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 09:12:14 2022

@author: Riccardo
"""

import pandas as pd

data5digits = pd.read_csv('CL_Dataset5digits.csv')
data5digits = data5digits.drop(['Cl', 'Cm'], axis=1)

data5digits['ID'] = data5digits['ID'].str.replace('NACA', '')

# make the new columns using string indexing
data5digits['des_cl'] = data5digits['ID'].str[0:1]
data5digits['max_pos'] = data5digits['ID'].str[1:3]
data5digits['thick']   = data5digits['ID'].str[3:5]

# get rid of the extra variable (if you want)
data5digits.drop('ID', axis=1, inplace=True)

# removing leading zeros from thick column
data5digits['thick'] = data5digits['thick'].str.lstrip('0')

data5digits = data5digits[~data5digits['des_cl'].astype(str).str.startswith('1')]
data5digits = data5digits[~data5digits['des_cl'].astype(str).str.startswith('3')]
data5digits = data5digits[~data5digits['des_cl'].astype(str).str.startswith('4')]
data5digits = data5digits[~data5digits['des_cl'].astype(str).str.startswith('5')]
data5digits = data5digits[~data5digits['des_cl'].astype(str).str.startswith('6')]

data5digits.to_csv('NACA5digits_preprocessed_new.csv')