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

plt.rcdefaults()
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
olddatafn = '/home/bill/cfg/mscL_MG1655_seqsbins.csv'
maindir = os.path.expanduser('~/Documents/energymatrix/mscS4815')
resultsfn = glob.glob(maindir + '/results/*.sql')
infofn = glob.glob(maindir + '/runsdetails/*.txt')
resultsfn.sort()
infofn.sort()

#Super Crude Error Checking!!
namedict = {}
for i, name in enumerate(resultsfn):
    fn = name.split('/results/')[1].split('.sql')[0]
    try:
        ind = next(i for i in range(0,len(infofn)) if fn in infofn[i])
        namedict[name] = ind
    except:
        namedict[name] = 'Fail'
        print name + ': there is no info!'
'''
    if fn not in infofn[i]:
        sys.exit(str(i) + 'details do not match files!')
'''        



masteroutputdir = os.path.expanduser('~/Documents/energymatrix/mscS4815/results/analyzeddata/')

start_dict = {}
csvfile = open(barcodefn,'r')
reader = csv.DictReader(csvfile)
for row in reader:
    start_dict[row['experiment_name']] = row['startseq']


# file with optimized matrix
for i, fn in enumerate(resultsfn):
    try:
        info_dict = plottingutils.readinfo(infofn[namedict[fn]])
    except:
        continue
    '''
    if 'conditional' in infofn[namedict[fn]]:
            continue
    '''    
    fnname = fn.split('/results/')[1].split('.sql')[0]
    # load the matrix, set each column to zero
    outputdir = os.path.join(masteroutputdir,fnname)
    burn_in = 1000
    try:
        db = pymc.database.sqlite.load(fn)
    except:
        continue
    
    #if len(db.trace('emat')[:]) != 3000:
    #    print 'Database is not complete'
    #    continue
    try:
        os.makedirs(outputdir)
    except:
        print 'already done'
        continue
    burn_in = 1000
    db = pymc.database.sqlite.load(fn)
    # emat_mean = db.emat.stats()['mean']
    try:
        emat_mean = sp.mean(db.trace('emat')[burn_in:],axis=0)
    
        # change the sign of emat_mean if necessary. We want a negative
        # correlation between energy and batch number (because the higher
        # the energy, the lower the expression)
        if 'old' in fnname:
            seq_mat,batch_vec = MCMC_utils.load_unique_seqs_batches(olddatafn,18 + int(info_dict['mut_region_start']),int(info_dict['mut_region_length']))
        else:
            try:
                unique = info_dict['unique']
                print 'is unique'
                seq_mat, batch_vec = plottingutils.loadseqsuniquemut(infofn[namedict[fn]],barcodefn,datafnbase)
            except:
                try:
                    cond = info_dict['condbase']
                    print 'is conditional'
                    seq_mat,batch_vec = plottingutils.loadcond(infofn[namedict[fn]],barcodefn,datafnbase)
                except:
                        print 'is not unique'
                        seq_mat, batch_vec = plottingutils.loadseqsuniquemut(infofn[namedict[fn]],barcodefn,datafnbase)
          
        energies = sp.zeros(len(batch_vec))
        for u in range(len(batch_vec)):
            energies[u] = sp.sum(seq_mat[:,:,u]*emat_mean)
        r = sp.stats.pearsonr(energies,batch_vec)[0]
        if r>0:
            emat_mean = MCMC_utils.fix_matrix_gauge(-emat_mean)
        else:
            emat_mean = MCMC_utils.fix_matrix_gauge(emat_mean)
        emat_mean = emat_mean*-1
        sp.savetxt(os.path.join(outputdir,fnname+'_emat_mean.txt'),emat_mean)
        
    
        MI,f_reg = MCMC_utils.compute_MI(seq_mat,batch_vec,emat_mean)
        MI_f = open(os.path.join(outputdir,fnname+'_MI.txt'),'w')
        MI_f.write(str(MI))
        
        plt.imshow(f_reg,interpolation='nearest',aspect='auto')
        plt.title('Joint regularized pdf, ' + fnname + ', MI: %.5f' % MI)
        plt.xlabel('Rank order')
        plt.ylabel('Batch number')
        plt.savefig(os.path.join(outputdir,fnname+'_regpdf.png'))
        plt.close()
        skip = 100
        thinned_trace = db.trace('emat')[::skip]
        n_samples = len(thinned_trace)
        MI_vec = sp.zeros(n_samples)
        for u in range(n_samples):
            MI_vec[u],foo = MCMC_utils.compute_MI(seq_mat,batch_vec,thinned_trace[u])
        plt.clf()
        plt.plot(MI_vec)
        plt.xlabel('MCMC iteration x '+str(skip))
        plt.ylabel('Mutual information')
        plt.title('MI, ' + fnname + ' MI: %.5f' % MI )
        plt.savefig(os.path.join(outputdir,fnname+'_MItrace.png'))
        plt.close()
        emat0 = emat_mean
        emat = emat0 - emat0.min(axis=0)
        

        # now make the plot
        site_seq = plottingutils.getwtseq(infofn[namedict[fn]],barcodefn,datafnbase)


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
        if 'conditional' in infofn[namedict[fn]]:
            condbase = str(int(start_dict[info_dict['exp_name']]) + int(info_dict['condbase']))
            plt.text(0.5,1.25,'Conditional Matrix: ' + info_dict['condident'] + ' at ' + condbase, horizontalalignment='center', fontsize=12, transform = ax2.transAxes)
        else:
            plt.text(0.5,1.25,'Energy Matrix: ' + str(tick_start) + ' to ' + str(tick_end),horizontalalignment='center', fontsize=12, transform = ax2.transAxes)
        plt.show()
        plt.savefig(os.path.join(outputdir,fnname + 'emat' + str(tick_start) + ' to ' + str(tick_end) + '.pdf'))
        plt.close()
    except:
        print fn + 'failed'
        continue

