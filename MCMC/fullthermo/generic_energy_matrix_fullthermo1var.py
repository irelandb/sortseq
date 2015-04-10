import scipy as sp
import pymc
import scipy.ndimage
import MCMC_utils_mscs
import MCMC_utils
import ConfigParser
import os
import numpy as np
import laccfg

import readcollatedremoveoligos
import MCMC_utils_1var

config = ConfigParser.RawConfigParser()
# everything that needs to be changed goes here ########

##################
config.read(os.path.expanduser(laccfg.cfg_fn))

mut_region_start = config.getint('Input','mut_region_start')
mut_region_length = config.getint('Input','mut_region_length')
seq_start = config.getint('Input','seqstart')
seq_end = config.getint('Input','seqend')
datab0 = config.get('Input','data_fn1')
datab01 = config.get('Input','data_fn2')
datab02 = config.get('Input','data_fn3')
datab03 = config.get('Input','data_fn4')
datab04 = config.get('Input','data_fn5')

# initialize random energy matrix
emat_0 = MCMC_utils_mscs.fix_matrix_gauge(sp.randn(4,mut_region_length + 21)) #RNAP emat from 0 to mut_region length, CRP emat from mut_region length to mut_region_length + 20, all other parameters in last collumn
emat_0[:,-1] = np.array([-5,-3,2.5,2.5]) #interaction energy, Additive constant to CRP matrix, multiplicative constant for CRP, mutliplicative constant for RNAP
# load in the data


numseq = [[] for i in range(0,4)]
seq_mat = [[] for i in range(0,4)]

sequences0 = readcollatedremoveoligos.collatedmat(datab0)
sequences1 = readcollatedremoveoligos.collatedmat(datab01)
sequences2 = readcollatedremoveoligos.collatedmat(datab02)
sequences3 = readcollatedremoveoligos.collatedmat(datab03)
sequences4 = readcollatedremoveoligos.collatedmat(datab04)

print 'sequences loaded'
print len(sequences0)
seq_start2 = 29
seq_end2 = 49

sequences0_2 = [sequences0[i][seq_start2:seq_end2] for i in range(0,len(sequences0))] 
sequences0 = [sequences0[i][mut_region_start:mut_region_start + mut_region_length] for i in range(0,len(sequences0))]

sequences1_2 = [sequences1[i][seq_start2:seq_end2] for i in range(0,len(sequences1))]
sequences1 = [sequences1[i][mut_region_start:mut_region_start + mut_region_length] for i in range(0,len(sequences1))]

sequences2_2 = [sequences2[i][seq_start2:seq_end2] for i in range(0,len(sequences2))]
sequences2 = [sequences2[i][mut_region_start:mut_region_start + mut_region_length] for i in range(0,len(sequences2))]

sequences3_2 = [sequences3[i][seq_start2:seq_end2] for i in range(0,len(sequences3))]
sequences3 = [sequences3[i][mut_region_start:mut_region_start + mut_region_length] for i in range(0,len(sequences3))]

sequences4_2 = [sequences4[i][seq_start2:seq_end2] for i in range(0,len(sequences4))]
sequences4 = [sequences4[i][mut_region_start:mut_region_start + mut_region_length] for i in range(0,len(sequences4))]




    

   
batchvec0 = [0 for i in range(0,len(sequences0))]
batchvec1 = [1 for i in range(0,len(sequences1))]
batchvec2 = [2 for i in range(0,len(sequences2))]
batchvec3 = [3 for i in range(0,len(sequences3))]
batchvec4 = [4 for i in range(0,len(sequences4))]
batch_vec_temp = batchvec0 + batchvec1 + batchvec2 + batchvec3 + batchvec4

seqs = sequences0 + sequences1 + sequences2 + sequences3 + sequences4

seqs2 = sequences0_2 + sequences1_2 + sequences2_2 + sequences3_2 + sequences4_2

#batch_vec_temp = [batch_vec_temp[i] for i in range(0,len(seqs)) if seqs[i].count('A') > 2 and len(seqs[i]) == mut_region_length]
print len(batch_vec_temp)
#seqs = [seqs[i] for i in range(0,len(seqs)) if seqs[i].count('A') > 2 and len(seqs[i]) == mut_region_length]
print len(seqs)

seq_mat_temp = np.empty([4,len(seqs[1]),len(seqs)])
seq_mat_temp2 = np.empty([4,len(seqs2[1]),len(seqs2)])

for i, line in enumerate(seqs):
    seq_mat_temp[:,:,i] = MCMC_utils.seq2mat(line)
    
for i, line in enumerate(seqs2):
    seq_mat_temp2[:,:,i] = MCMC_utils.seq2mat(line)
# Run matrix on only section of data



# shuffle the elements of seq_mat and batch_vec. This will prevent
# spuriously high mutual information values
print len
index_shuf = range(len(batch_vec_temp))
sp.random.shuffle(index_shuf)
seq_mat = sp.zeros([4,len(seq_mat_temp[0,:,0]),len(seq_mat_temp[0,0,:])],dtype = 'int')
seq_mat2 = sp.zeros([4,len(seq_mat_temp2[0,:,0]),len(seq_mat_temp2[0,0,:])],dtype = 'int')
batch_vec = sp.zeros_like(batch_vec_temp)
for i, i_s in enumerate(index_shuf):
    seq_mat[:,:,i] = seq_mat_temp[:,:,i_s]
    seq_mat2[:,:,i] = seq_mat_temp2[:,:,i_s]
    batch_vec[i] = batch_vec_temp[i_s]
    




###############################################################################
############### Stochastics go here ###########################################
###############################################################################
@pymc.stochastic(observed=True,dtype=int)
def sequences(value=seq_mat):
    return 0
@pymc.stochastic(observed=True,dtype=int)
def sequences2(value=seq_mat2):
    return 0
@pymc.stochastic(observed=True,dtype=int)
def batches(value=batch_vec):
    return 0

'''
def emat(value=emat_0):
    return 0
    
@pymc.stochastic(dtype=float)
def emat2(value=emat_0):
    return 0
        
@pymc.stochastic(dtype=float)
def rnapnum(value=2000):
    return 0
    
@pymc.stochastic(dtype=float)
def repnum(s=sequences,b=batches, em = emat, em2 = emat2, p = rnapnum, value=100):
    n_seqs = len(b)
    MI, f_reg = MCMC_utils_experiment.compute_MI(s,b,em, em2, value,p)
    return n_seqs*MI
'''
@pymc.stochastic(dtype=float)
def emat(s=sequences,s2 = sequences2,b=batches,value=emat_0):
    n_seqs = len(b)
    MI, f_reg = MCMC_utils_1var.compute_MI(s,s2,b,value)
    return n_seqs*MI


