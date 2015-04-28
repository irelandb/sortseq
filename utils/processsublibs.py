# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 20:08:42 2015

@author: bill
"""
from __future__ import division
import sys,os
import seaborn as sns
sns.set_context('notebook', font_scale=1.5)
sys.path.append('/home/bill/sortseq/utils')
sys.path.append('/home/bill/Documents/energymatrix/lacstardata')
sys.path.append(os.path.expanduser('~/kinneylab_software/'))
import CorrectNSBseqinfo
import glob

import numpy as np
import readcollatedlasso
import matplotlib.pyplot as plt

names = glob.glob('35_sub+cAMP*.fasta.unique')
names.sort()


    
numseq = []   
seqs = []
for fn in names:
    sequences,ns = readcollatedlasso.collatedmat(fn)
    ns = list(ns)
    seqs.append(sequences)
    numseq.append(ns)

stopped = 0
seqs_list = []


for z,seq in enumerate(seqs[0][0:300]):
    n = []
    for u in range(1,len(seqs)):
        try:
            n.append(numseq[u][next(i for i in range(0,len(seqs[u])) if seq == seqs[u][i])])
        except StopIteration:
            n.append(0)
    seqs_list.append([seq] + n)

pmax = [np.max(seqs_list[i][1:])/np.sum(seqs_list[i][1:]) for i in range(len(seqs_list))]
plt.hist(pmax)
plt.title('35_sub-cAMP Sequence Distribution')
plt.xlabel('Max percentage in one bin')
plt.ylabel('Number of sequecnes')


