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
pairnum = 2
config = ConfigParser.RawConfigParser()
# everything that needs to be changed goes here ########

##################
config.read(os.path.expanduser(mscscfg.cfg_fn))

mut_region_start = config.getint('Input','mut_region_start')
mut_region_length = config.getint('Input','mut_region_length')
seq_start = config.getint('Input','seqstart')
seq_end = config.getint('Input','seqend')
datab0 = config.get('Input','data_fn1')
datab01 = config.get('Input','data_fn2')
datab02 = config.get('Input','data_fn3')
datab03 = config.get('Input','data_fn4')

# initialize random energy matrix
emat_0 = MCMC_utils_mscs.fix_matrix_gauge(sp.randn(16,mut_region_length-1))
# load in the data


numseq = [[] for i in range(0,4)]
seq_mat = [[] for i in range(0,4)]

sequences0 = readcollatedremoveoligos.collatedmat(datab0)
sequences1 = readcollatedremoveoligos.collatedmat(datab01)
sequences2 = readcollatedremoveoligos.collatedmat(datab02)
sequences3 = readcollatedremoveoligos.collatedmat(datab03)

print 'sequences loaded'
print len(sequences0)
    
sequences0 = [sequences0[i][seq_start:seq_end] for i in range(0,len(sequences0))]

sequences1 = [sequences1[i][seq_start:seq_end] for i in range(0,len(sequences1))]

sequences2 = [sequences2[i][seq_start:seq_end] for i in range(0,len(sequences2))]

sequences3 = [sequences3[i][seq_start:seq_end] for i in range(0,len(sequences3))]


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
batch_vec_temp = batchvec0 + batchvec1 + batchvec2 + batchvec3

seqs = a + a1 + a2 + a3

#batch_vec_temp = [batch_vec_temp[i] for i in range(0,len(seqs)) if seqs[i].count('A') > 2 and len(seqs[i]) == mut_region_length]
print len(batch_vec_temp)
#seqs = [seqs[i] for i in range(0,len(seqs)) if seqs[i].count('A') > 2 and len(seqs[i]) == mut_region_length]
print len(seqs)

seq_mat_temp = np.empty([4,len(seqs[1]),len(seqs)])
seq_mat_temp2 = np.zeros([16,len(seqs[1])-1,len(seqs)])
seq_dict = {'A':0,'C':1,'G':2,'T':3}
for i, line in enumerate(seqs):
    seq_mat_temp[:,:,i] = MCMC_utils.seq2mat(line)
    for q in range(0,len(line)-1):
        pos = seq_dict[line[q]]*4 + seq_dict[line[np.mod(q+pairnum,20)]]
        seq_mat_temp2[pos,q,i] = 1
    


        

# Run matrix on only section of data



# shuffle the elements of seq_mat and batch_vec. This will prevent
# spuriously high mutual information values
print len
index_shuf = range(len(batch_vec_temp))
sp.random.shuffle(index_shuf)
seq_mat2 = sp.zeros([16,len(seq_mat_temp2[0,:,0]),len(seq_mat_temp2[0,0,:])],dtype = 'int')
batch_vec = sp.zeros_like(batch_vec_temp)
for i, i_s in enumerate(index_shuf):
    seq_mat2[:,:,i] = seq_mat_temp2[:,:,i_s]
    batch_vec[i] = batch_vec_temp[i_s]

###############################################################################
############### Stochastics go here ###########################################
###############################################################################
@pymc.stochastic(observed=True,dtype=int)
def sequences(value=seq_mat2):
    return 0

@pymc.stochastic(observed=True,dtype=int)
def batches(value=batch_vec):
    return 0

@pymc.stochastic(dtype=float)
def emat(s=sequences,b=batches,value=emat_0):
    n_seqs = len(b)
    MI, f_reg = MCMC_utils_mscs.compute_MI(s,b,value)
    return n_seqs*MI
