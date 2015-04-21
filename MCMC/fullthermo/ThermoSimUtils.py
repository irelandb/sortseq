# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 16:16:55 2015

@author: bill
"""
import numpy as np
import sys
sys.path.append('/home/bill/kinneylab_software')
import scipy as sp
import mpmath as mp
import nsbestimator as nsb
def genseqs(wtseq,numseq):
    seq_dict = {'A':0,'C':1,'G':2,'T':3}
    parr = np.zeros([len(wtseq),4])
    for i,let in enumerate(wtseq):
        parr[i,:] = np.roll(np.array([.7,.1,.1,.1]),seq_dict[let])
    bp = ['A','C','G','T']
    seqs = []
    for i in range(numseq):
            s = [np.random.choice(bp,1,p=parr[q,:])[0] for q in range(len(wtseq))]
            seqs.append(''.join(s))
    return seqs
    
def seq2mat(seq): # Returns a 4 by L matrix where L is the length of the seq, and each 1 represents which base is at that position
    seq_dict = {'A':0,'C':1,'G':2,'T':3}
    mat = sp.zeros((4,len(seq)),dtype=int)
    for i,bp in enumerate(seq):
        mat[seq_dict[bp],i] = 1
    return mat
def seq2mat2(seq): # returns which parameters are true for the sequence, where each parameter is a possible pairing
    pair_dict = {'AA':0,'AC':1, 'AG':2,'AT':3,'CA':4,'CC':5,'CG':6,'CT':7,'GA':8,'GC':9,'GG':10,'GT':11,'TA':12,'TC':13,'TG':14,'TT':15}
    pairlist = []
    lseq = len(seq)
    index = 0
    for i,bp in enumerate(seq):
        for z in range(i+1,lseq):
            pairlist.append(index*16 + pair_dict[bp + seq[z]])
            index = index +1
    return pairlist
def seq2mat3(seq):
    pair_dict = {'A':0,'C':1, 'G':2,'T':3}
    trilist = []
    lseq = len(seq)
    index = 0
    for i,bp in enumerate(seq):
        for z in range(i+1,lseq-1):
            for q in range(z+1,lseq):
                trilist.append(index*64 + pair_dict[bp]*16 + pair_dict[seq[z]]*4 + pair_dict[seq[q]])
                index = index +1
    return trilist    
def genmodels(fn,seqs,modeltype): #Parses seqs and model type then calculates and returns energies
    emat = np.genfromtxt(fn,skiprows=1)
    f = open(fn)
    mattype = f.read()[:6]
    energies = np.zeros(len(seqs))
    N = len(seqs)
    mut_region_length = len(seqs[0])
    if modeltype == mattype:
        if mattype == '1Point':
            lasso_mat = sp.sparse.lil_matrix((N,4*len(seqs[0])))
            for i,s in enumerate(seqs):
                seq_mat = seq2mat(s)
                energies[i] = (seq_mat*emat).sum()
                lasso_mat[i,:] = seq_mat.ravel(order='F')
        elif mattype == '2Point':
            lasso_mat = sp.sparse.lil_matrix((N,round(sp.misc.comb(mut_region_length,2))*16))
            for i,s in enumerate(seqs):
                lasso_mat[i,seq2mat2(s)] = 1
                energies[i] = (lasso_mat[i,:]*(emat.ravel())).sum()
        elif mattype == '3Point':
            lasso_mat = sp.sparse.lil_matrix((N,round(sp.misc.comb(mut_region_length,3))*64))
            for i,s in enumerate(seqs):
                lasso_mat[i,seq2mat3(s)] = 1
                energies[i] = (lasso_mat[i,:]*(emat.ravel())).sum()
        else:
            raise TypeError('Model Name Must be 1Point, 2Point, or 3Point')
    else:
        if mattype == '1Point':
            for i,s in enumerate(seqs):
                seq_mat = seq2mat(s)
                energies[i] = (seq_mat*emat).sum()
        elif mattype == '2Point':
            for i,s in enumerate(seqs):
                seq_mat = np.zeros(round(sp.misc.comb(mut_region_length,2))*16)
                seq_mat[seq2mat2(s)] = 1
                energies[i] = (seq_mat*(emat.ravel())).sum()
        elif mattype == '3Point':
            for i,s in enumerate(seqs):
                seq_mat = np.zeros(round(sp.misc.comb(mut_region_length,3))*64)
                seq_mat[seq2mat3(s)] = 1
                energies[i] = (seq_mat*(emat.ravel())).sum()
        else:
            raise TypeError('Model Name Must be 1Point, 2Point, or 3Point')
        if modeltype == '1Point':
            lasso_mat = sp.sparse.lil_matrix((N,4*len(seqs[0])))
            for i,s in enumerate(seqs):
                seq_mat = seq2mat(s)
                lasso_mat[i,:] = seq_mat.ravel(order='F')
        elif modeltype == '2Point':
            lasso_mat = sp.sparse.lil_matrix((N,round(sp.misc.comb(mut_region_length,2))*16))
            for i,s in enumerate(seqs):
                lasso_mat[i,seq2mat2(s)] = 1
        elif modeltype == '3Point':
            lasso_mat = sp.sparse.lil_matrix((N,round(sp.misc.comb(mut_region_length,3))*64))
            for i,s in enumerate(seqs):
                lasso_mat[i,seq2mat3(s)] = 1
            
    return sp.sparse.csr_matrix(lasso_mat),energies

def genmat(seqs,modeltype): #Parses seqs and model type then calculates and returns energies
    N = len(seqs)
    mut_region_length = len(seqs[0])
    if modeltype == '1Point':
            seq_mat = np.zeros([4,len(seqs[0]),N])
            for i,s in enumerate(seqs):
                seq_mat[:,:,i] = seq2mat(s)
    '''
    elif modeltype == '2Point':
            lasso_mat = sp.sparse.lil_matrix((N,round(sp.misc.comb(mut_region_length,2))*16))
            for i,s in enumerate(seqs):
                lasso_mat[i,seq2mat2(s)] = 1
    elif modeltype == '3Point':
            lasso_mat = sp.sparse.lil_matrix((N,round(sp.misc.comb(mut_region_length,3))*64))
            for i,s in enumerate(seqs):
                lasso_mat[i,seq2mat3(s)] = 1
    '''
    return seq_mat

def genenergies(fnR,fnQ,seqsR,seqsQ,gamma,sQ,sR,R0): #Parses seqs and model type then calculates and returns energies R is transcription factor, Q is RNAP
    ematR = np.genfromtxt(fnR,skiprows=1)
    ematQ = np.genfromtxt(fnQ,skiprows=1)
    fR = open(fnR)
    fQ = open(fnQ)
    mattype = fR.read()[:6] #mattype must be the same
    #mattypeQ = fQ.read()[:6]
    energies = np.zeros(len(seqsQ))
    N = len(seqsQ)
    mut_region_lengthQ = len(seqsQ[0])
    mut_region_lengthR = len(seqsR[0])
    
    if mattype == '1Point':
            for i,s in enumerate(seqsR):
                seq_matR = seq2mat(s)
		seq_matQ = seq2mat(seqsQ[i])
		RNAP = (seq_matQ*ematQ).sum()*sQ
		TF = (seq_matR*ematR).sum()*sR + R0
                energies[i] = -RNAP + mp.log(1 + mp.exp(-TF - gamma)) - mp.log(1 + mp.exp(-TF))
    '''
    elif mattype == '2Point':
            for i,s in enumerate(seqs):
                seq_mat = np.zeros(round(sp.misc.comb(mut_region_length,2))*16)
                seq_mat[seq2mat2(s)] = 1
                energies[i] = (seq_mat*(emat.ravel())).sum()
    elif mattype == '3Point':
            for i,s in enumerate(seqs):
                seq_mat = np.zeros(round(sp.misc.comb(mut_region_length,3))*64)
                seq_mat[seq2mat3(s)] = 1
                energies[i] = (seq_mat*(emat.ravel())).sum()
    '''
    return energies

def compute_MI_orig(seq_matQ,seq_matR,batches,ematQ,ematR,gamma,sR,sQ,R_0):
    # preliminaries
    n_seqs = len(batches)
    n_batches = int(batches.max()) + 1 # assumes zero indexed batches
    n_bins = 1000
    
    #energies = sp.zeros(n_seqs)
    f = sp.zeros((n_batches,n_seqs))
    
    # compute energies
    # for i in range(n_seqs):
    #     energies[i] = sp.sum(seqs[:,:,i]*emat)
    # alternate way
    energies = np.zeros(n_seqs)
    for i in range(n_seqs):
    	RNAP = (seq_matQ[:,:,i]*ematQ).sum()*sQ
    	TF = (seq_matR[:,:,i]*ematR).sum()*sR + R_0
    	energies[i] = -RNAP + mp.log(1 + mp.exp(-TF - gamma)) - mp.log(1 + mp.exp(-TF))


    # sort energies
    inds = sp.argsort(energies)
    for i,ind in enumerate(inds):
        f[batches[ind],i] = 1.0/n_seqs # batches aren't zero indexed
        

    # bin and convolve with Gaussian
    f_binned = sp.zeros((n_batches,n_bins))
    
    for i in range(n_batches):
        f_binned[i,:] = sp.histogram(f[i,:].nonzero()[0],bins=n_bins,range=(0,n_seqs))[0]
    #f_binned = f_binned/f_binned.sum()
    f_reg = sp.ndimage.gaussian_filter1d(f_binned,0.04*n_bins,axis=1)
    f_reg = f_reg/f_reg.sum()

    # compute marginal probabilities
    p_b = sp.sum(f_reg,axis=1)
    p_s = sp.sum(f_reg,axis=0)

    # finally sum to compute the MI
    MI = 0
    for i in range(n_batches):
        for j in range(n_bins):
            if f_reg[i,j] != 0:
                MI = MI + f_reg[i,j]*sp.log2(f_reg[i,j]/(p_b[i]*p_s[j]))
    print MI
    return MI,f_reg

def compute_MIlasso(seqs,batches,matrixcoefs,n_bins):
    # preliminaries
    n_seqs = len(batches)
    n_batches = int(batches.max()) + 1 # assumes zero indexed batches
    
    #energies = sp.zeros(n_seqs)
    f = sp.zeros((n_batches,n_seqs))
    
    # compute energies
    # for i in range(n_seqs):
    #     energies[i] = sp.sum(seqs[:,:,i]*emat)
    # alternate way
    energies = np.zeros(n_seqs)
    energiestemp = seqs.multiply(matrixcoefs).sum(axis=1)
    for i in range(0,n_seqs):
	energies[i] = energiestemp[i]
    #for i in range(0,n_seqs):    
    #    energies[i] = (np.array(matrixcoefs.todense())*np.array(seqs[i,:].todense())).sum(1)[0]

        
    


    # sort energies
    inds = sp.argsort(energies)
    for i,ind in enumerate(inds):
        f[batches[ind],i] = 1.0/n_seqs # batches aren't zero indexed
        

    # bin and convolve with Gaussian
    f_binned = sp.zeros((n_batches,n_bins))
    MI = 0
    for i in range(n_batches):
        f_binned[i,:] = sp.histogram(f[i,:].nonzero()[0],bins=n_bins,range=(0,n_seqs))[0]
    #f_binned = f_binned/f_binned.sum()
    #f_reg = scipy.ndimage.gaussian_filter1d(f_binned,0.04*n_bins,axis=1)
    #f_reg = f_reg/f_reg.sum()

    # compute marginal probabilities
    p_b = sp.sum(f_binned,axis=1)
    originalent = nsb.S(p_b,n_seqs,n_batches)
    p_b = p_b/n_seqs
    p_s = sp.sum(f_binned,axis=0)/n_seqs
    H_nsb = [nsb.S(f_binned[:,i],f_binned[:,i].sum(),len(f_binned[:,i])) for i in range(0,n_bins)]
    H_mean = (p_s*H_nsb).sum()
    MI = originalent - H_mean
    # finally sum to compute the MI
    #MI = 0
    #for i in range(n_batches):
    #    for j in range(n_bins):
    #        if f_reg[i,j] != 0:
    #            MI = MI + f_reg[i,j]*sp.log2(f_reg[i,j]/(p_b[i]*p_s[j]))
    return MI

    
    
    
