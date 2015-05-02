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


f = open('/home/wireland/mscS4-8-15/runsdetails/numtestinfo.txt','r')
raw = f.read().split('-')
condbase = int(raw[0])
condident = raw[1]

numbins = 4
config = ConfigParser.RawConfigParser()
# everything that needs to be changed goes here ########

config.read('/home/wireland/mscS4-8-15/data/mscL.cfg') #location of config file, change this to run on different datasets
mut_region_start = config.getint('Input','mut_region_start') #not base pair number, distance from start of mut region
mut_region_length = config.getint('Input','mut_region_length')
data_fnbase = config.get('Input','data_fnbase')
expname = config.get('Input','expname') #ex MscS mut1, describes the experiment without the different batch numbers.
fnnames = glob.glob(data_fnbase + expname + '*.fasta')
fnnames.sort()

#this section determines from the barcode file where the mutated region starts for this oligo
barcodefn = config.get('Input','barcodefn')

barcode_dict = {}
reverse_dict = {}
csvfile = open(barcodefn,'r')
reader = csv.DictReader(csvfile)
for row in reader:
    barcode_dict[row['experiment_name']] = row['fwd_barcode']
    reverse_dict[row['experiment_name']] = row['rev_barcode']


numseq = [[] for i in range(0,4)]
seq_mat = [[] for i in range(0,4)]

#read in sequences from files
sequences = [readuniqueseqssingleend.collatedmat(fn) for fn in fnnames]
seq_start = len(barcode_dict[expname])
seq_end = sequences[0][0].find(reverse_dict[expname][0:5])
print 'sequences loaded'
print len(sequences[0])
for i in range(0,numbins):
    sequences[i] = [sequences[i][z][seq_start:seq_end] for z in range(0,len(sequences[i]))]

batch_vec_temp = []
seqs = []
#filter sequences for only unique sequences in the target range.
for i in range(0,numbins):
    tempseqs = sequences[i]
    s2 = list(set([tempseqs[z][mut_region_start:mut_region_start + mut_region_length] for z in range(0,len(tempseqs)) if tempseqs[z][condbase] == condident]))
    seqs = seqs + s2
    batch_vec_temp = batch_vec_temp + [i for z in range(len(s2))]

batch_vec_temp = np.array(batch_vec_temp)


print len(batch_vec_temp)

print len(seqs)

seq_mat_temp = np.empty([4,len(seqs[1]),len(seqs)])

if condbase >= mut_region_start and condbase < mut_region_start + mut_region_length:
	for i, line in enumerate(seqs):
    		seq_mat_temp[:,:,i] = MCMC_utils.seq2mat(line)
    		seq_mat_temp[:,condbase-mut_region_start,i] = MCMC_utils.seq2mat(np.random.choice(['A','C','G','T'])).transpose()
else:
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
