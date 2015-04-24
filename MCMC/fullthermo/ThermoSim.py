# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 16:43:46 2015

@author: bill
"""
from __future__ import division
import pymc
import sys, os
sys.path.append(os.path.expanduser('~'))
sys.path.append(os.path.expanduser('~/sortseq'))
sys.path.append(os.path.expanduser('~/sortseq/MCMC/fullthermo'))
sys.path.append(os.path.expanduser('~/sortseq/utils'))

import ThermoSimUtils
import numpy as np
import scipy as sp

import MCMC_utils_mscs
N = 10000
seqsQ = ThermoSimUtils.genseqs('CAG',N)
seqsR = ThermoSimUtils.genseqs('CAG',N)
fn = os.path.expanduser('~/sortseq/MCMC/fullthermo/etest.txt')
seq_matR = ThermoSimUtils.genmat(seqsR,'1Point')
seq_matQ = ThermoSimUtils.genmat(seqsQ,'1Point')
energies = ThermoSimUtils.genenergies(fn,fn,seqsR,seqsQ,-5,-5,-6,-.1)
myscale = np.abs(energies)*.2
expression = energies + np.random.normal(scale=myscale,size=len(energies))
inds = sp.argsort(expression)
batch_vec = np.zeros(N)
batchvectemp = [i for i in range(0,4) for z in range(0,int(np.ceil(len(energies)/4)))]
for q,i in enumerate(inds):
        batch_vec[i] = batchvectemp[q]
mut_region_lengthQ = len(seqsQ[0])
mut_region_lengthR = len(seqsR[0])
        
        
ematR_0 = MCMC_utils_mscs.fix_matrix_gauge(sp.randn(4,mut_region_lengthR))       
ematQ_0 = MCMC_utils_mscs.fix_matrix_gauge(sp.randn(4,mut_region_lengthQ))

gamma = pymc.Uniform('gamma',-10,10,observed=False)
sR = pymc.Uniform('sR',-10,0,observed=False)
sQ = pymc.Uniform('sQ',-10,0,observed=False)
R_0 = pymc.Uniform('R_0',-2,0,observed=False)

@pymc.stochastic(observed=True,dtype=int)
def sequencesR(value=seq_matR):
    return 0

@pymc.stochastic(observed=True,dtype=int)
def sequencesQ(value=seq_matQ):
    return 0

@pymc.stochastic(observed=True,dtype=int)
def batches(value=batch_vec):
    return 0

@pymc.stochastic(dtype=float)
def ematQ(value=ematQ_0):
    return 0

@pymc.stochastic(dtype=float)
def ematR(seqQ=sequencesQ,seqR=sequencesR,b=batches,value=ematR_0,eQ=ematQ,g = gamma, scaleR = sR, scaleQ = sQ,numR = R_0):
    n_seqs = len(b)
    MI, f_reg = ThermoSimUtils.compute_MI_orig(seqQ,seqR,b,eQ,value,g,scaleR,scaleQ,numR)
    return n_seqs*MI        
