# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 00:17:36 2015

@author: bill
"""
from __future__ import division
import sys
sys.path.append('/home/bill/Documents/energymatrix/mscsstuff/')
sys.path.append('/home/bill/Documents/energymatrix/lacstardata/')
sys.path.append('/home/bill/Documents/DiversityData/')
import MCMC_utils
import ConfigParser
import glob
import numpy as np
import csv

import matplotlib.pyplot as plt
import readuniqueseqssingleend
config = ConfigParser.RawConfigParser()
config.read('/home/bill/Documents/energymatrix/mscS4815/mscL.cfg')
mut_region_start = config.getint('Input','mut_region_start')
mut_region_length = config.getint('Input','mut_region_length')
data_fnbase = config.get('Input','data_fnbase')
expname = config.get('Input','expname')
fnnames = glob.glob(data_fnbase + expname + '*.fasta')
fnnames.sort()
barcodefn = config.get('Input','barcodefn')

start_dict = {}
end_dict = {}
csvfile = open(barcodefn,'r')
reader = csv.DictReader(csvfile)
for row in reader:
    start_dict[row['experiment_name']] = row['startseq']
    end_dict[row['experiment_name']] = row['endseq']
    
import glob
fnnames = glob.glob('/home/bill/Documents/energymatrix/mscS4815/infofootprints/*conditional*')

for name in fnnames:
    info = np.load(name).transpose()
    titlename = name.split('conditional')[1].strip('.npy')
    start = int(start_dict[titlename])
    end = int(end_dict[titlename])    
    if 'renorm' in name:
        titlename = titlename + ' renorm'
        plt.imshow(info,interpolation='nearest',extent=[start,end,end,start],vmin=0,vmax=.7)
    else:
        plt.imshow(info,interpolation='nearest',extent=[start,end,end,start])
    plt.title(titlename + ' conditional')
    plt.xlabel('Position')
    plt.xlabel('Conditional Position')
    plt.savefig('/home/bill/Documents/energymatrix/mscS4815/infofootprints/' + titlename + '.pdf')
    plt.close()
    
    