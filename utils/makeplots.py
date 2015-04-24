# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 15:09:50 2015

@author: bill
"""
from __future__ import division
import sys, os
sys.path.append('/home/bill/Documents/energymatrix')
sys.path.append(os.path.expanduser('~/sortseq/utils'))
import scipy as sp
import matplotlib.pyplot as plt
#import seaborn as sns
import MCMC_utils
import ConfigParser
import glob
import numpy as np
import csv
import readuniqueseqssingleend
import os
import plottingutils
import glob
import pymc
from matplotlib.mlab import find
from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost


def clean_up_xticklabels(start, length):
    # convenience function to make nice xticklabels. We want the
    # labels on only multiples of 5
    axis_coords = range(length)
    TSS_coords = range(start,start+length)
    indices = find([sp.remainder(coord,5)==0 for coord in TSS_coords])
    xtick_labels = [str(TSS_coords[ind]) for ind in indices]

    return indices, xtick_labels
    
numbins = 4
config = ConfigParser.RawConfigParser()
cfgname = os.path.expanduser('~/sortseq/cfgfiles/mscL.cfg')
config.read(cfgname)
#datafnbase = config.get('Input','data_fnbase')
#barcodefn = config.get('Input','barcodefn')
datafnbase = '/home/bill/Documents/energymatrix/mscS4815/'
barcodefn = '/home/bill/Documents/energymatrix/mscS4815/barcodes.csv'
maindir = os.path.expanduser('~/Documents/energymatrix/mscS4815')
resultsfn = glob.glob(maindir + '/results/*.sql')
infofn = glob.glob(maindir + '/runsdetails/*.txt')
resultsfn.sort()
infofn.sort()

#Super Crude Error Checking!!
for i, name in enumerate(resultsfn):
    fn = name.split('/results/')[1].split('.sql')[0]
    if fn not in infofn[i]:
        sys.exit('details do not match files!')
        
wtseqs = [plottingutils.getwtseq(infofn[i],barcodefn,datafnbase) for i in range(0,len(resultsfn))]


masteroutputdir = os.path.expanduser('~/Documents/energymatrix/mscS4815/results/analyzeddata/')

start_dict = {}
csvfile = open(barcodefn,'r')
reader = csv.DictReader(csvfile)
for row in reader:
    start_dict[row['experiment_name']] = row['startseq']


# file with optimized matrix
for i, fn in enumerate(resultsfn):
    info_dict = plottingutils.readinfo(infofn[i])
    fnname = fn.split('/results/')[1].split('.sql')[0]
    # load the matrix, set each column to zero
    outputdir = os.path.join(masteroutputdir,fnname)
    try:
        os.makedirs(outputdir)
    except:
        pass
    burn_in = 1000
    db = pymc.database.sqlite.load(fn)
    # emat_mean = db.emat.stats()['mean']
    emat_mean = sp.mean(db.trace('emat')[burn_in:],axis=0)
    
    # change the sign of emat_mean if necessary. We want a negative
    # correlation between energy and batch number (because the higher
    # the energy, the lower the expression)
    seq_mat, batch_vec = plottingutils.loadseqsunique(infofn[i],barcodefn,datafnbase)
    energies = sp.zeros(len(batch_vec))
    for i in range(len(batch_vec)):
        energies[i] = sp.sum(seq_mat[:,:,i]*emat_mean)
    r = sp.stats.pearsonr(energies,batch_vec)[0]
    if r>0:
        emat_mean = MCMC_utils.fix_matrix_gauge(-emat_mean)
    else:
        emat_mean = MCMC_utils.fix_matrix_gauge(emat_mean)
    sp.savetxt(os.path.join(outputdir,fnname+'_emat_mean.txt'),emat_mean)
    
    
    MI,f_reg = MCMC_utils.compute_MI(seq_mat,batch_vec,emat_mean)
    MI_f = open(os.path.join(outputdir,fnname+'_MI.txt'),'w')
    MI_f.write(str(MI))
    
    plt.imshow(f_reg,interpolation='nearest',aspect='auto')
    plt.title('Joint regularized pdf, ' + fnname + ', MI: %.5f' % MI)
    plt.xlabel('Rank order')
    plt.ylabel('Batch number')
    plt.savefig(os.path.join(outputdir,fnname+'_regpdf.png'))
    
    skip = 100
    thinned_trace = db.trace('emat')[::skip]
    n_samples = len(thinned_trace)
    MI_vec = sp.zeros(n_samples)
    for i in range(n_samples):
        MI_vec[i],foo = MCMC_utils.compute_MI(seq_mat,batch_vec,thinned_trace[i])
    plt.clf()
    plt.plot(MI_vec)
    plt.xlabel('MCMC iteration x '+str(skip))
    plt.ylabel('Mutual information')
    plt.title('MI, ' + fnname + ' MI: %.5f' % MI )
    plt.savefig(os.path.join(outputdir,fnname+'_MItrace.png'))
    
    emat0 = emat_mean
    emat = emat0 - emat0.min(axis=0)
    

    # now make the plot
    site_seq = wtseqs[i]



    fig = plt.figure()
    ax1 = SubplotHost(fig, 1,1,1)
    fig.add_subplot(ax1)


    ax2 = ax1.twin()
    ax1.imshow(emat,interpolation='nearest')

    ax1.set_xlabel('Position w.r.t. transcription start site')
    ax1.set_yticks([0,1,2,3])
    ax1.set_yticklabels(['A','C','G','T'])


    # label positions with respect to transcription start site
    tick_start = int(start_dict[info_dict['exp_name']])+int(info_dict['mut_region_start'])
    tick_end = int(start_dict[info_dict['exp_name']])+int(info_dict['mut_region_start']) + int(info_dict['mut_region_length'])
    indices, xtick_labels = clean_up_xticklabels(tick_start,tick_end-tick_start)
    ax1.set_xticks(indices)
    ax1.set_xticklabels(xtick_labels)

    # put the sequence above it
    ax2.set_yticks([])
    ax2.set_yticklabels([])
    ax2.set_xticks(range(20))
    ax2.set_xticklabels([bp for bp in site_seq])
    #plt.title('$lac$ promoter, -35 region')
    plt.text(0.5,1.25,'Energy Matrix from ' + str(tick_start) + 'to' + str(tick_end),horizontalalignment='center',
         fontsize=16, transform = ax2.transAxes)

    plt.show()


