# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 15:21:48 2015

@author: bill
"""
from __future__ import division
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import MCMC_utils
import ConfigParser
import glob
import numpy as np
import csv
import readuniqueseqssingleend


def readinfo(infoname):
    infodict = {}
    if 'conditional' in infoname:
        f = open(infoname)
        raw = f.readlines()
        infodict['condbase'] = raw[0].split('-')[0]
        infodict['condident'] = raw[0].split('-')[1]
        for line in raw[1:]:
            infodict[line.split(' = ')[0]] = line.split(' = ')[1].strip('\n') 
    else:
        f = open(infoname)
        raw = f.readlines()
        for line in raw:
            infodict[line.split(' = ')[0]] = line.split(' = ')[1].strip('\n')
    
    return infodict

def loadseqsunique(infoname,barcodefn,datafnbase):
    infodict = readinfo(infoname)
    numbins = 4
    
    mut_region_start = int(infodict['mut_region_start'])
    mut_region_length = int(infodict['mut_region_length'])
    expname = infodict['exp_name']
    fnnames = glob.glob(datafnbase + expname + '*.fasta')
    fnnames.sort()

    barcode_dict = {}
    reverse_dict = {}
    csvfile = open(barcodefn,'r')
    reader = csv.DictReader(csvfile)
    for row in reader:
        barcode_dict[row['experiment_name']] = row['fwd_barcode']
        reverse_dict[row['experiment_name']] = row['rev_barcode']

    sequences = [readuniqueseqssingleend.collatedmat(fn) for fn in fnnames]
    seq_start = len(barcode_dict[expname])
    seq_end = sequences[0][0].find(reverse_dict[expname][0:5])
    print 'sequences loaded'
    print len(sequences[0])
    for i in range(0,numbins):
        sequences[i] = [sequences[i][z][seq_start:seq_end] for z in range(0,len(sequences[i]))]

    batch_vec_temp = []
    seqs = []
    for i in range(0,numbins):
        tempseqs = list(set(sequences[i]))
        seqs = seqs + [tempseqs[z][mut_region_start:mut_region_start + mut_region_length] for z in range(0,len(tempseqs))]
        batch_vec_temp = batch_vec_temp + [i for z in range(0,len(tempseqs))]

    batch_vec_temp = np.array(batch_vec_temp)

    #batch_vec_temp = [batch_vec_temp[i] for i in range(0,len(seqs)) if seqs[i].count('A') > 2 and len(seqs[i]) == mut_region_length]
    print len(batch_vec_temp)
    #seqs = [seqs[i] for i in range(0,len(seqs)) if seqs[i].count('A') > 2 and len(seqs[i]) == mut_region_length]
    print len(seqs)

    seq_mat_temp = np.empty([4,len(seqs[1]),len(seqs)])

    for i, line in enumerate(seqs):
        seq_mat_temp[:,:,i] = MCMC_utils.seq2mat(line)
    return seq_mat_temp, batch_vec_temp
    
def loadcond(infoname,barcodefn,datafnbase):
    infodict = readinfo(infoname)
    numbins = 4
    
    condbase = int(infodict['condbase'])
    condident = infodict['condident']
    
    mut_region_start = int(infodict['mut_region_start'])
    mut_region_length = int(infodict['mut_region_length'])
    expname = infodict['exp_name']
    fnnames = glob.glob(datafnbase + expname + '*.fasta')
    fnnames.sort()

    barcode_dict = {}
    reverse_dict = {}
    csvfile = open(barcodefn,'r')
    reader = csv.DictReader(csvfile)
    for row in reader:
        barcode_dict[row['experiment_name']] = row['fwd_barcode']
        reverse_dict[row['experiment_name']] = row['rev_barcode']

    sequences = [readuniqueseqssingleend.collatedmat(fn) for fn in fnnames]
    seq_start = len(barcode_dict[expname])
    seq_end = sequences[0][0].find(reverse_dict[expname][0:5])
    print 'sequences loaded'
    print len(sequences[0])
    for i in range(0,numbins):
        sequences[i] = [sequences[i][z][seq_start:seq_end] for z in range(0,len(sequences[i]))]

    batch_vec_temp = []
    seqs = []
    for i in range(0,numbins):
        tempseqs = list(set(sequences[i]))
        seqs = seqs + [tempseqs[z][mut_region_start:mut_region_start + mut_region_length] for z in range(0,len(tempseqs)) if tempseqs[z][condbase] == condident]
        batch_vec_temp = batch_vec_temp + [i for z in range(0,len(tempseqs)) if tempseqs[z][condbase] == condident]

    batch_vec_temp = np.array(batch_vec_temp)

    #batch_vec_temp = [batch_vec_temp[i] for i in range(0,len(seqs)) if seqs[i].count('A') > 2 and len(seqs[i]) == mut_region_length]
    print len(batch_vec_temp)
    #seqs = [seqs[i] for i in range(0,len(seqs)) if seqs[i].count('A') > 2 and len(seqs[i]) == mut_region_length]
    print len(seqs)

    seq_mat_temp = np.empty([4,len(seqs[1]),len(seqs)])

    if condbase >= mut_region_start and condbase < mut_region_start + mut_region_length:
	for i, line in enumerate(seqs):
    		seq_mat_temp[:,:,i] = MCMC_utils.seq2mat(line)
    		seq_mat_temp[:,condbase-mut_region_start,i] = MCMC_utils.seq2mat(np.random.choice(['A','C','G','T'])).transpose()
    else:
	for i, line in enumerate(seqs):
    		seq_mat_temp[:,:,i] = MCMC_utils.seq2mat(line)
    return seq_mat_temp, batch_vec_temp

def loadseqsuniquemut(infoname,barcodefn,datafnbase):
    infodict = readinfo(infoname)
    numbins = 4
    
    mut_region_start = int(infodict['mut_region_start'])
    mut_region_length = int(infodict['mut_region_length'])
    expname = infodict['exp_name']
    fnnames = glob.glob(datafnbase + expname + '*.fasta')
    fnnames.sort()

    barcode_dict = {}
    reverse_dict = {}
    csvfile = open(barcodefn,'r')
    reader = csv.DictReader(csvfile)
    for row in reader:
        barcode_dict[row['experiment_name']] = row['fwd_barcode']
        reverse_dict[row['experiment_name']] = row['rev_barcode']

    sequences = [readuniqueseqssingleend.collatedmat(fn) for fn in fnnames]
    seq_start = len(barcode_dict[expname])
    seq_end = sequences[0][0].find(reverse_dict[expname][0:5])
    print 'sequences loaded'
    print len(sequences[0])
    for i in range(0,numbins):
        sequences[i] = [sequences[i][z][seq_start:seq_end] for z in range(0,len(sequences[i]))]

    batch_vec_temp = []
    seqs = []
    for i in range(0,numbins):
        #tempseqs = list(set(sequences[i])) I instead added this to the next line to make sure only sequences with unique mutated regions are counted.
        tempseqs = sequences[i]
        s2 = list(set([tempseqs[z][mut_region_start:mut_region_start + mut_region_length] for z in range(0,len(tempseqs))]))
        seqs = seqs + s2
        batch_vec_temp = batch_vec_temp + [i for z in range(len(s2))]

    batch_vec_temp = np.array(batch_vec_temp)

    #batch_vec_temp = [batch_vec_temp[i] for i in range(0,len(seqs)) if seqs[i].count('A') > 2 and len(seqs[i]) == mut_region_length]
    print len(batch_vec_temp)
    #seqs = [seqs[i] for i in range(0,len(seqs)) if seqs[i].count('A') > 2 and len(seqs[i]) == mut_region_length]
    print len(seqs)

    seq_mat_temp = np.empty([4,len(seqs[1]),len(seqs)])

    for i, line in enumerate(seqs):
        seq_mat_temp[:,:,i] = MCMC_utils.seq2mat(line)
    return seq_mat_temp, batch_vec_temp
    
def getwtseq(infoname,barcodefn,datafnbase):
    
    seq_mat_temp, bat = loadseqsunique(infoname,barcodefn,datafnbase)
    seqdict = {0:'A', 1:'C', 2:'G', 3:'T'}
    wtseq = ''
    for z in range(0,len(seq_mat_temp[0,:,0])):
        #for z in range(0,21):
        Apos = np.nonzero(seq_mat_temp[0,z,:])[0]
        Cpos = np.nonzero(seq_mat_temp[1,z,:])[0]
        Gpos = np.nonzero(seq_mat_temp[2,z,:])[0]
        Tpos = np.nonzero(seq_mat_temp[3,z,:])[0]
        Asum = len(Apos)
        Csum = len(Cpos)
        Gsum = len(Gpos)
        Tsum = len(Tpos)
        wtseq = wtseq + seqdict[np.argmax(np.array([Asum,Csum,Gsum,Tsum]))]
        
    return wtseq
            
            
        