# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 23:55:52 2015

@author: bill
"""
from __future__ import division
import emcee
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
energies = ThermoSimUtils.genenergies(fn,fn,seqsR,seqsQ,-5,5,6,.1)
myscale = np.abs(energies)*.2
expression = energies + np.random.normal(scale=myscale,size=len(energies))
inds = sp.argsort(expression)
batch_vec = np.zeros(N)
batchvectemp = [i for i in range(0,4) for z in range(0,int(np.ceil(len(energies)/4)))]
for q,i in enumerate(inds):
        batch_vec[i] = batchvectemp[q]
mut_region_lengthQ = len(seqsQ[0])
mut_region_lengthR = len(seqsR[0])

nwalkers = 58
ntemps = 1
TotalVals = np.zeros([ntemps,nwalkers,4*(mut_region_lengthR + mut_region_lengthQ)+1])   
LR = 4*mut_region_lengthR
LQ = 4*mut_region_lengthQ
for w in range(ntemps):  
    for i in range(0,nwalkers):        
        ematR_0 = MCMC_utils_mscs.fix_matrix_gauge(sp.randn(4,mut_region_lengthR)).ravel(order='F')  
        ematQ_0 = MCMC_utils_mscs.fix_matrix_gauge(sp.randn(4,mut_region_lengthQ)).ravel(order='F')
        g_0 = np.random.rand(1)*-10
        #sR_0 = np.random.rand(1)*-10
        #sQ_0 = np.random.rand(1)*-10
        #Rin_0 = np.random.rand(1)*-1
        TotalVals[w,i,:LR] = ematR_0
        TotalVals[w,i,LR:LR+LQ] = ematQ_0
        #TotalVals[i,-4:] = np.array([g_0,sR_0,sQ_0,Rin_0]).transpose()
        TotalVals[w,i,-1] = g_0
    
def lnprob(vals):
    MI, f_reg = ThermoSimUtils.compute_MI_origemcee(seq_matQ,seq_matR,batch_vec,vals[LR:LR+LQ].reshape([4,mut_region_lengthQ],order='F'),vals[:LR].reshape([4,mut_region_lengthR],order='F'),vals[-1])
    return MI
    
def logp(x):
    return 0.0
    
sampler = emcee.PTSampler(ntemps, nwalkers, len(TotalVals[0,0,:]), lnprob,logp)

sampler.run_mcmc(TotalVals,10)



