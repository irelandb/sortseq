# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 22:31:36 2015

@author: bill
"""
import MCMC_utils

import numpy as np
emat40A = np.genfromtxt('/home/bill/Documents/energymatrix/mscS4815/results/analyzeddata/conditional40A_40/conditional40A_40_emat_mean.txt')
emat40T = np.genfromtxt('/home/bill/Documents/energymatrix/mscS4815/results/analyzeddata/conditional40T_40/conditional40T_40_emat_mean.txt')

wtseq = 'TTATTGTTTACCCTTGTCAG'

wtmat = MCMC_utils.seq2mat(wtseq)

emat40Anorm = emat40A - (wtmat*emat40A).max(axis=0)
emat40Tnorm = emat40T - (wtmat*emat40T).max(axis=0)