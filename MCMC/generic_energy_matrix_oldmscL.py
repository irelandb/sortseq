import sys
sys.path.append('/home/wireland/')
sys.path.append('/home/wireland/sortseq/utils/')
import scipy as sp
import pymc
import scipy.ndimage
import glob
import MCMC_utils_mscs
import MCMC_utils
import ConfigParser
import glob
import mscscfg
import os
import numpy as np
import csv

import readuniqueseqssingleend



numbins = 4
config = ConfigParser.RawConfigParser()
# everything that needs to be changed goes here ########

config.read('/home/wireland/mscS4-8-15/data/oldmscL.cfg') #location of config file, change this to run on different datasets
mut_region_start = config.getint('Input','mut_region_start') #not base pair number, distance from start of mut region
mut_region_length = config.getint('Input','mut_region_length')
data_fn = config.get('Input','data_fn')
seq_start = config.getint('Input','seq_start')
seq_end = config.getint('Input','seq_end')
#data_fnbase = config.get('Input','data_fnbase')
#expname = config.get('Input','expname') #ex MscS mut1, describes the experiment without the different batch numbers.
#fnnames = glob.glob(data_fnbase + expname + '*.fasta')
#fnnames.sort()

#this section determines from the barcode file where the mutated region starts for this oligo
'''
barcodefn = config.get('Input','barcodefn')

barcode_dict = {}
reverse_dict = {}
csvfile = open(barcodefn,'r')
reader = csv.DictReader(csvfile)
for row in reader:
    barcode_dict[row['experiment_name']] = row['fwd_barcode']
    reverse_dict[row['experiment_name']] = row['rev_barcode']
'''

f = open(data_fn)
# read lines into one big list and transform into a set. This
# automatically gets rid of duplicates
# lines with region of interest selected
roi_list = [(line.split(',')[0][mut_region_start:mut_region_start+mut_region_length], line.split(',')[1].strip()) for line in f if line.strip()]
f.close()
N = len(roi_list)
index_shuf = range(N)
batch_vec_temp = np.array([float(roi_list[z][1]) for z in index_shuf])
batch_vec_temp = batch_vec_temp - batch_vec_temp.min()

seqs = [roi_list[z][0] for z in index_shuf]
seqs = [seqs[z][seq_start:seq_end] for z in range(0,len(seqs))]
seqs = [seqs[z][mut_region_start:mut_region_start + mut_region_length] for z in range(0,len(seqs))]
seqstemp = []
batch_vec2 = []
for i in range(numbins):
        indexes = np.nonzero(batch_vec_temp == i)[0]
	s = []
	for u in indexes:
	     s.append(seqs[u])
	seqstemp = seqstemp + list(set(s))
	batch_vec2 = batch_vec2 + [i for z in range(len(seqstemp))]
seqs = seqstemp
batch_vec_temp = batch_vec2


print len(batch_vec_temp)

print len(seqs)

seq_mat_temp = np.empty([4,len(seqs[1]),len(seqs)])


for i, line in enumerate(seqs):
    seq_mat_temp[:,:,i] = MCMC_utils.seq2mat(line)

#initial energy matrix
emat_0 = MCMC_utils.fix_matrix_gauge(sp.randn(4,mut_region_length))




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
