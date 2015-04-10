import scipy as sp
import scipy.ndimage
import readcollatedmscs
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

def load_seqs_batches(data_fn,mut_region_start,mut_region_length):
    N = 0
    f = open(data_fn)
    for line in f:
        if line.strip():
            N = N + 1
    f.close()
    print N

    # instantiate batch and sequence variables
    batch_vec = sp.zeros(N,dtype=int)
    seq_mat = sp.zeros((4,mut_region_length,N),dtype=int)

    f = open(data_fn)
    for i, line in enumerate(f):
        if line.strip():
            sb = line.split(',')
            batch_vec[i] = int(sb[1])
            seq_mat[:,:,i] = seq2mat(sb[0][mut_region_start:mut_region_start+mut_region_length])
    f.close()
    batch_vec = batch_vec-batch_vec.min()

    return seq_mat,batch_vec

def load_unique_seqs_batches(data_fn1,data_fn2,data_fn3,data_fn4,mut_region_start,mut_region_length):
    """Load in unique sequence-batche pairs from data file.

    INPUTS:
    data_fn: csv file containing sequence, batch
    mut_region_start: sequence index corresponding to start of ROI
    mut_region_length: self-evident

    OUTPUTS:
    seq_mat: 4xmut_region_lengthxN matrix containing sequence information
    batch_vec: N length vector containing batches
    """
    
    seq_mat1,numseq1 = readcollatedmscs.collatedmat(data_fn1)
    seq_mat2,numseq2 = readcollatedmscs.collatedmat(data_fn2)
    seq_mat3,numseq3 = readcollatedmscs.collatedmat(data_fn3)
    seq_mat4,numseq4 = readcollatedmscs.collatedmat(data_fn4)
    
    numseq1 = np.array([0 for i in range(0, len(numseq1))])
    numseq2 = np.array([1 for i in range(0, len(numseq2))])
    numseq3 = np.array([2 for i in range(0, len(numseq3))])
    numseq4 = np.array([3 for i in range(0, len(numseq4))])
    batch_vec = np.hstack((numseq1,numseq2,numseq3,numseq4))
    
    # instantiate batch and sequence variables
    
    seq_mat = np.dstack((seq_mat1,seq_mat2,seq_mat3,seq_mat4))
    print 'seq_mat done'
   
    
    
    return seq_mat,batch_vec

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
