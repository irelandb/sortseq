from __future__ import division
import MCMC_utils
import scipy as sp
import pymc
import scipy.ndimage
import ConfigParser
import os
import numpy as np
import matplotlib.pyplot as plt
numbins = 4
config = ConfigParser.RawConfigParser()
config.read('/home/bill/Documents/energymatrix/danieltest/runs_analysis/oldtest35_full-cAMP/35_full-cAMP.cfg')

mut_region_start = config.getint('Input','mut_region_start')
mut_region_length = config.getint('Input','mut_region_length')
data_fn = config.get('Input','data_fn')
seq_mat_temp, batch_vec_temp = MCMC_utils.load_unique_seqs_batches(data_fn,mut_region_start,mut_region_length)

ifoot = np.zeros([len(seq_mat_temp[0,:,0])])
ifootrenorm = np.zeros([len(seq_mat_temp[0,:,0])])  
pbatch = np.array([float(batch_vec_temp.tolist().count(i))/float(len(batch_vec_temp)) for i in range(0,numbins)])
wtseq = ''
seqdict = {0:'A', 1:'C', 2:'G', 3:'T'}
Anumseq2 = [[] for i in range(0,len(seq_mat_temp[0,:,0]))]
Cnumseq2 = [[] for i in range(0,len(seq_mat_temp[0,:,0]))]
Gnumseq2 = [[] for i in range(0,len(seq_mat_temp[0,:,0]))]
Tnumseq2 = [[] for i in range(0,len(seq_mat_temp[0,:,0]))]
mutrate = np.zeros(mut_region_length)
for z in range(0,len(seq_mat_temp[0,:,0])):
#for z in range(0,21):
    Apos = np.nonzero(seq_mat_temp[0,z,:])[0]
    Cpos = np.nonzero(seq_mat_temp[1,z,:])[0]
    Gpos = np.nonzero(seq_mat_temp[2,z,:])[0]
    Tpos = np.nonzero(seq_mat_temp[3,z,:])[0]
    Anumseq = [batch_vec_temp[Apos].tolist().count(i) for i in range(0,numbins)]
    Cnumseq = [batch_vec_temp[Cpos].tolist().count(i) for i in range(0,numbins)]
    Gnumseq = [batch_vec_temp[Gpos].tolist().count(i) for i in range(0,numbins)]
    Tnumseq = [batch_vec_temp[Tpos].tolist().count(i) for i in range(0,numbins)]
    Anumseq2[z] = [batch_vec_temp[Apos].tolist().count(i) for i in range(0,numbins)]
    Cnumseq2[z] = [batch_vec_temp[Cpos].tolist().count(i) for i in range(0,numbins)]
    Gnumseq2[z] = [batch_vec_temp[Gpos].tolist().count(i) for i in range(0,numbins)]
    Tnumseq2[z] = [batch_vec_temp[Tpos].tolist().count(i) for i in range(0,numbins)]
    Asum = np.array(np.sum(Anumseq))
    Csum = np.array(np.sum(Cnumseq))
    Gsum = np.array(np.sum(Gnumseq))
    Tsum = np.array(np.sum(Tnumseq))
    mutrate[z] = 1.0 - np.max([Asum,Csum,Gsum,Tsum])/np.sum([Asum,Csum,Gsum,Tsum])
    wtseq = wtseq + seqdict[np.argmax(np.array([Asum,Csum,Gsum,Tsum]))]
    try:
        condA = np.array([float(Anumseq[i])/float(Asum) for i in range(0,numbins)])
        condC = np.array([float(Cnumseq[i])/float(Csum) for i in range(0,numbins)])
        condG = np.array([float(Gnumseq[i])/float(Gsum) for i in range(0,numbins)])
        condT = np.array([float(Tnumseq[i])/float(Tsum) for i in range(0,numbins)])
        totsum = np.sum(np.array([Asum,Csum,Gsum,Tsum]))
        #print np.sum(-pbatch*np.log2(pbatch)) + np.sum((Asum/totsum)*(condA)*np.log2(condA) + Csum/totsum*(condC)*np.log2(condC) + Gsum/totsum*(condG)*np.log2(condG) +  Tsum/totsum*condT*np.log2(condT))
        ifoot[z] =  np.sum(-pbatch*np.log2(pbatch)) + np.sum((Asum/totsum)*(condA)*np.log2(condA) + Csum/totsum*(condC)*np.log2(condC) + Gsum/totsum*(condG)*np.log2(condG) +  Tsum/totsum*condT*np.log2(condT))
        ifootrenorm[z] = np.sum(-pbatch*np.log2(pbatch)) + .25*np.sum((condA)*np.log2(condA) + (condC)*np.log2(condC) + (condG)*np.log2(condG) +  condT*np.log2(condT))
    except:
        print z
        pass
fig = plt.figure()
ax = fig.add_subplot(111)

x = range(0, 0 + len(ifoot))
x[:] = [a - .45 for a in x]
numBins = 165
ax.bar(x,ifoot,color = 'b')
plt.xlabel('Base Position',fontsize = 20)
plt.ylabel('Mutual Information (bits)', fontsize = 20)
figtitle = 'Information Footprint mscL (Original Method)'
plt.title(figtitle, fontsize = 22)   
plt.show() 