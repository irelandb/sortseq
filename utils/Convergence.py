# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 23:29:13 2015

@author: bill
"""
from __future__ import division

import sys, os
sys.path.append('/home/bill/Documents/energymatrix')
sys.path.append(os.path.expanduser('~/sortseq/utils'))
import scipy as sp
import matplotlib.pyplot as plt
#import seaborn as sns
import MCMC_utils
import ConfigParser
import glob
import numpy as np
import csv
import readuniqueseqssingleend
import os
import plottingutils
import glob
import pymc
from matplotlib.mlab import find
from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost



        
numbins = 4
config = ConfigParser.RawConfigParser()
cfgname = os.path.expanduser('~/sortseq/cfgfiles/mscL.cfg')
config.read(cfgname)
datafnbase = config.get('Input','data_fnbase')
#barcodefn = config.get('Input','barcodefn')
#datafnbase = '/home/bill/Documents/energymatrix/35_full-cAMP/'
barcodefn = '/home/bill/Documents/energymatrix/mscS4815/barcodes.csv'
maindir = os.path.expanduser('~/Documents/energymatrix/mscS4815')
resultsfn = glob.glob(maindir + '/results/mscLpromoter*.sql')
infofn = glob.glob(maindir + '/runsdetails/mscLpromoter*.txt')
resultsfn.sort()
infofn.sort()

#Super Crude Error Checking!!
for i, name in enumerate(resultsfn):
    fn = name.split('/results/')[1].split('.sql')[0]
    if fn not in infofn[i]:
        sys.exit('details do not match files!')
        



#masteroutputdir = os.path.expanduser('~/Documents/energymatrix/mscS4815/results/analyzeddata/')

start_dict = {}
csvfile = open(barcodefn,'r')
reader = csv.DictReader(csvfile)
for row in reader:
    start_dict[row['experiment_name']] = row['startseq']

db = pymc.database.sqlite.load(resultsfn[0])
t = db.trace('emat')[0]
mut_region_length = len(t[0,:])
emeans = np.zeros([len(resultsfn),4,mut_region_length])
#emeans = np.zeros([len(resultsfn),4,int(info_test['mut_region_length'])])
estdMCMC = np.zeros_like(emeans)
evis = np.zeros([len(resultsfn),4*mut_region_length])
# file with optimized matrix
for i, fn in enumerate(resultsfn):
    #info_dict = plottingutils.readinfo(infofn[i])
    #fnname = fn.split('/results/')[1].split('.sql')[0]
    # load the matrix, set each column to zero
    #outputdir = os.path.join(masteroutputdir,fnname)
    #try:
    #    os.makedirs(outputdir)
    #except:
    #    print 'already done'
    #    continue
    burn_in = 1000
    db = pymc.database.sqlite.load(fn)
    # emat_mean = db.emat.stats()['mean']
    if len(db.trace('emat')[:]) == 3000:
        ematchain = db.trace('emat')[burn_in:]
        emeans[i,:,:] = MCMC_utils.fix_matrix_gauge(np.mean(ematchain,axis=0))
        if emeans[i,0,0] > 0:
            emeans[i,:,:] = emeans[i,:,:]*-1
        #emeans[i,:,:] = emeans[i,:,:] - emeans[i,:,:].min(axis=0)
        evis[i,:] = emeans[i,:,:].ravel(order='F')
        estdMCMC[i,:,:] = np.std(ematchain,axis=0)
'''        
std1 = np.std(emeans[:len(resultsfn)/2,:,:],axis=0)
emean1 = np.mean(emeans[:len(resultsfn)/2,:,:],axis=0)
std2 = np.std(emeans[len(resultsfn)/2:,:,:],axis=0)
emean2 = np.mean(emeans[len(resultsfn)/2:,:,:],axis=0)
convarr = abs(emean1-emean2) / np.sqrt(std1**2 + std2**2)
'''

      