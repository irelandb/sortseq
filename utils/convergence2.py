# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 18:08:10 2015

@author: bill
"""
import numpy as np
import glob

fnnames = glob.glob('/home/bill/Documents/energymatrix/mscS4815/results/analyzeddata/mut2minus10*/*_emat*')
condname = glob.glob('/home/bill/Documents/energymatrix/mscS4815/results/analyzeddata/conditional40*_*/*_emat*')
l = np.genfromtxt(fnnames[0]).shape[1]
emats = np.zeros([4,l,len(fnnames)])
ematvis = np.zeros([4*l,len(fnnames)])
condemats = np.zeros([4,l,len(condname)])
for i,n in enumerate(condname):
    condemats[:,:,i] = np.genfromtxt(n)
    
for i,name in enumerate(fnnames):
    emats[:,:,i] = np.genfromtxt(name)[:,:l]
    ematvis[:,i] = np.genfromtxt(name)[:,:l].ravel(order='F')