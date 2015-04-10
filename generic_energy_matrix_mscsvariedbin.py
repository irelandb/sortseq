import scipy as sp
import pymc
import scipy.ndimage
import MCMC_utils_mscs
import MCMC_utils
import ConfigParser
import mscscfg
import os
import numpy as np

import readcollatedremoveoligos
import readuniqueseqssingleend
config = ConfigParser.RawConfigParser()
# everything that needs to be changed goes here ########

##################
config.read(mscscfg.cfg_fn)

mut_region_start = config.getint('Input','mut_region_start')
mut_region_length = config.getint('Input','mut_region_length')
seq_start = config.getint('Input','seqstart')
seq_end = config.getint('Input','seqend')
datab0 = config.get('Input','data_fn1')
datab01 = config.get('Input','data_fn2')
datab02 = config.get('Input','data_fn3')
datab03 = config.get('Input','data_fn4')

# initialize random energy matrix
emat_0 = MCMC_utils_mscs.fix_matrix_gauge(sp.randn(4,mut_region_length))
# load in the data

numseq = [[] for i in range(0,4)]
seq_mat = [[] for i in range(0,4)]

sequences0,u = readuniqueseqssingleend.collatedmat(datab0)
sequences1,u = readuniqueseqssingleend.collatedmat(datab01)
sequences2,u = readuniqueseqssingleend.collatedmat(datab02)
sequences3,u = readuniqueseqssingleend.collatedmat(datab03)
numbins = 4
sequences = [sequences0,sequences1,sequences2,sequences3]

print 'sequences loaded'
print len(sequences0)
for i in range(0,numbins):
    sequences[i] = [sequences[i][z][seq_start:seq_end] for z in range(0,len(sequences[i])) if sequences[i][z][seq_start:seq_end].count('N') == 0 and len(sequences[i][z][seq_start:seq_end]) == (seq_end - seq_start)]

batch_vec_temp = []
seqs = []
for i in range(0,numbins):
    tempseqs = list(set(sequences[i]))
    seqs = seqs + [tempseqs[z][mut_region_start:mut_region_start + mut_region_length] for z in range(0,len(tempseqs))]
    batch_vec_temp = batch_vec_temp + [i for z in range(0,len(tempseqs))]

batch_vec_temp = np.array(batch_vec_temp)
'''  
a = list(set(sequences0))

a = [a[i][mut_region_start:mut_region_start + mut_region_length] for i in range(0,len(a))]

a1 = list(set(sequences1))
a1 = [a1[i][mut_region_start:mut_region_start + mut_region_length] for i in range(0,len(a1))]

a2 = list(set(sequences2))
a2 = [a2[i][mut_region_start:mut_region_start + mut_region_length] for i in range(0,len(a2))]

a3 = list(set(sequences3))
a3 = [a3[i][mut_region_start:mut_region_start + mut_region_length] for i in range(0,len(a3))]



    

   
batchvec0 = [0 for i in range(0,len(a))]
batchvec1 = [1 for i in range(0,len(a1))]
batchvec2 = [2 for i in range(0,len(a2))]
batchvec3 = [3 for i in range(0,len(a3))]

batch_vec_temp = np.array(batchvec0 + batchvec1 + batchvec2 + batchvec3)

seqs = a + a1 + a2 + a3 
'''
#batch_vec_temp = [batch_vec_temp[i] for i in range(0,len(seqs)) if seqs[i].count('A') > 2 and len(seqs[i]) == mut_region_length]
print len(batch_vec_temp)
#seqs = [seqs[i] for i in range(0,len(seqs)) if seqs[i].count('A') > 2 and len(seqs[i]) == mut_region_length]
print len(seqs)

seq_mat_temp = np.empty([4,len(seqs[1]),len(seqs)])


for i, line in enumerate(seqs):
    seq_mat_temp[:,:,i] = MCMC_utils.seq2mat(line)


# Run matrix on only section of data



# shuffle the elements of seq_mat and batch_vec. This will prevent
# spuriously high mutual information values
print len
index_shuf = range(len(batch_vec_temp))
sp.random.shuffle(index_shuf)
seq_mat = sp.zeros([4,len(seq_mat_temp[0,:,0]),len(seq_mat_temp[0,0,:])],dtype = 'int')
batch_vec = sp.zeros_like(batch_vec_temp)
for i, i_s in enumerate(index_shuf):
    seq_mat[:,:,i] = seq_mat_temp[:,:,i_s]
    batch_vec[i] = batch_vec_temp[i_s]

###############################################################################
############### Stochastics go here ###########################################
###############################################################################
@pymc.stochastic(observed=True,dtype=int)
def sequences(value=seq_mat):
    return 0

@pymc.stochastic(observed=True,dtype=int)
def batches(value=batch_vec):
    return 0

@pymc.stochastic(dtype=float)
def emat(s=sequences,b=batches,value=emat_0):
    n_seqs = len(b)
    MI, f_reg = MCMC_utils_mscs.compute_MI(s,b,value)
    return n_seqs*MI
