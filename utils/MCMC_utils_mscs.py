import scipy as sp
import scipy.ndimage
import numpy as np

def compute_energies(seqs,batches,emat):
    """seqs: matrix of sequences, should be 4xLxN
    batches: vector of batches
    emat: energy matrix, 4xL"""
    dot = emat[:,:,sp.newaxis]*seqs
    energies = dot.sum(0).sum(0)
    return energies

def zero_matrix(emat):
    # set the smallest element of each column to zero
    for j in range(emat.shape[1]):
        emat[:,j] = emat[:,j] - emat[:,j].min() + 1
        
    return emat

def fix_matrix_gauge(emat):
    """Fix gauge of an energy matrix such that the average value
    of each column is zero (columns correspond to positions), and
    overall matrix norm is equal to 1."""
    # fix mean
    for j in range(emat.shape[1]):
        emat[:,j] = emat[:,j] -sp.mean(emat[:,j])
    # fix inner product
    emat = emat/sp.sqrt(sp.sum(emat*emat))
    return emat

def seq2mat(seq):
    seq_dict = {'A':0,'C':1,'G':2,'T':3}
    mat = sp.zeros((4,len(seq)),dtype=int)
    for i,bp in enumerate(seq):
        mat[seq_dict[bp],i] = 1
    return mat

def compute_MI(seqs,batches,emat):
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
    dot = emat[:,:,sp.newaxis]*seqs
    energies = dot.sum(0).sum(0)


    # sort energies
    inds = sp.argsort(energies)
    for i,ind in enumerate(inds):
        f[batches[ind],i] = 1.0/n_seqs # batches aren't zero indexed
        

    # bin and convolve with Gaussian
    f_binned = sp.zeros((n_batches,n_bins))
    
    for i in range(n_batches):
        f_binned[i,:] = sp.histogram(f[i,:].nonzero()[0],bins=n_bins,range=(0,n_seqs))[0]
    #f_binned = f_binned/f_binned.sum()
    f_reg = scipy.ndimage.gaussian_filter1d(f_binned,0.04*n_bins,axis=1)
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

def determine_sign_of_emat(emat,wt_seq):
    """determine what the correct sign is for an energy matrix. We will
    use the assumption that the wild type sequence must be better
    binding than a random sequence.

    INPUTS:
    emat: energy matrix
    wt_seq: wild type sequence of energy matrix

    OUTPUT:
    emat: energy matrix with correct sign
    """
    n_rand = 1000 # number of random sequences to check
    e_rand = sp.zeros(n_rand)
    # convert sequence to matrix
    seq_mat = seq2mat(wt_seq)
    e_wt = sp.sum(emat*seq_mat)

    for i in range(n_rand):
        seq_rand = sp.zeros((4,len(wt_seq)))

        for j in range(len(wt_seq)):
            seq_rand[sp.random.randint(4),j] = 1
        e_rand[i] = sp.sum(emat*seq_rand)
    if e_wt < sp.mean(e_rand):
        return emat
    else:
        return -emat

def get_PSSM_from_weight_matrix(emat,factor):
    """get position specific scoring matrix from weight matrix. There
    is an undetermined scale factor, which JBK suggests manually
    adjusting until getting a reasonable information content (say 1
    bit per bp).

    Assumes that negative energies -> better binding.
    """
    
    # need to reverse sign for PSSM
    emat = -emat
    # set lowest element to zero
    emat = emat - emat.min(axis=0)
    # exponentiate
    p = sp.exp(factor*emat)
    p = p/p.sum(axis=0)
    return p

def compute_PSSM_self_information(p):
    """Compute self information of a PSSM. See the wikipedia page on
    PSSMs, for instance."""
    return -sp.sum(p*sp.log(p))
