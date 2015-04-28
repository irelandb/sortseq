from __future__ import division
import sys
sys.path.append('/home/bill/Documents/energymatrix/mscsstuff/')
sys.path.append('/home/bill/Documents/energymatrix/lacstardata/')
sys.path.append('/home/bill/Documents/DiversityData/')
import MCMC_utils
import ConfigParser
import glob
import numpy as np
import csv

import matplotlib.pyplot as plt
import readuniqueseqssingleend

''' This script generates information footprints'''
numbins = 4
config = ConfigParser.RawConfigParser()
# everything that needs to be changed goes here ########
allexp = ['MscS mut1', 'MscS mut3','Mut2','Mut3','MscL mut1 new(sorted)']
##################
#config.read(os.path.expanduser(mscscfg.cfg_fn))
config.read('/home/bill/Documents/energymatrix/mscS4815/mscL.cfg')
mut_region_start = config.getint('Input','mut_region_start')
mut_region_length = config.getint('Input','mut_region_length')
data_fnbase = config.get('Input','data_fnbase')
#expname = config.get('Input','expname')

barcodefn = config.get('Input','barcodefn')

barcode_dict = {}
reverse_dict = {}
csvfile = open(barcodefn,'r')
reader = csv.DictReader(csvfile)
for row in reader:
    barcode_dict[row['experiment_name']] = row['fwd_barcode']
    reverse_dict[row['experiment_name']] = row['rev_barcode']
for expname in allexp:
    fnnames = glob.glob(data_fnbase + expname + '*.fasta')
    fnnames.sort()
    numseq = [[] for i in range(0,4)]
    seq_mat = [[] for i in range(0,4)]

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

    ifoot = np.zeros([len(seqs[0]),len(seqs[0])*4])    
    ifootrenorm = np.zeros([len(seqs[0]),len(seqs[0])*4])  
    for condbase in range(0,len(seqs[0])):
        for p, condident in enumerate(['A','C','G','T']):
            #filter for only sequences equal to our condbase
            seqscond = []
            batchveccond = []
            for i in range(0,len(seqs)):
                if seqs[i][condbase] == condident: 
                    seqscond.append(seqs[i])
                    batchveccond.append(batch_vec_temp[i])
        
            batchveccond = np.array(batchveccond)

            #batch_vec_temp = [batch_vec_temp[i] for i in range(0,len(seqs)) if seqs[i].count('A') > 2 and len(seqs[i]) == mut_region_length]
            print len(batchveccond)
            #seqs = [seqs[i] for i in range(0,len(seqs)) if seqs[i].count('A') > 2 and len(seqs[i]) == mut_region_length]
            print len(seqscond)

            seq_mat_temp = np.empty([4,len(seqscond[1]),len(seqscond)])


            for i, line in enumerate(seqscond):
                seq_mat_temp[:,:,i] = MCMC_utils.seq2mat(line)

            pbatch = np.array([float(batchveccond.tolist().count(i))/float(len(batchveccond)) for i in range(0,numbins)])
            wtseq = ''
            seqdict = {0:'A', 1:'C', 2:'G', 3:'T'}
            Anumseq2 = [[] for i in range(0,len(seq_mat_temp[0,:,0]))]
            Cnumseq2 = [[] for i in range(0,len(seq_mat_temp[0,:,0]))]
            Gnumseq2 = [[] for i in range(0,len(seq_mat_temp[0,:,0]))]
            Tnumseq2 = [[] for i in range(0,len(seq_mat_temp[0,:,0]))]
            for z in range(0,len(seq_mat_temp[0,:,0])):
                #for z in range(0,21):
                Apos = np.nonzero(seq_mat_temp[0,z,:])[0]
                Cpos = np.nonzero(seq_mat_temp[1,z,:])[0]
                Gpos = np.nonzero(seq_mat_temp[2,z,:])[0]
                Tpos = np.nonzero(seq_mat_temp[3,z,:])[0]
                Anumseq = [batchveccond[Apos].tolist().count(i) for i in range(0,numbins)]
                Cnumseq = [batchveccond[Cpos].tolist().count(i) for i in range(0,numbins)]
                Gnumseq = [batchveccond[Gpos].tolist().count(i) for i in range(0,numbins)]
                Tnumseq = [batchveccond[Tpos].tolist().count(i) for i in range(0,numbins)]
                Anumseq2[z] = [batchveccond[Apos].tolist().count(i) for i in range(0,numbins)]
                Cnumseq2[z] = [batchveccond[Cpos].tolist().count(i) for i in range(0,numbins)]
                Gnumseq2[z] = [batchveccond[Gpos].tolist().count(i) for i in range(0,numbins)]
                Tnumseq2[z] = [batchveccond[Tpos].tolist().count(i) for i in range(0,numbins)]
                Asum = np.array(np.sum(Anumseq))
                Csum = np.array(np.sum(Cnumseq))
                Gsum = np.array(np.sum(Gnumseq))
                Tsum = np.array(np.sum(Tnumseq))
                #wtseq = wtseq + seqdict[np.argmax(np.array([Asum,Csum,Gsum,Tsum]))]
                try:
                    condA = np.array([float(Anumseq[i])/float(Asum) for i in range(0,numbins)])
                    condC = np.array([float(Cnumseq[i])/float(Csum) for i in range(0,numbins)])
                    condG = np.array([float(Gnumseq[i])/float(Gsum) for i in range(0,numbins)])
                    condT = np.array([float(Tnumseq[i])/float(Tsum) for i in range(0,numbins)])
                    totsum = np.sum(np.array([Asum,Csum,Gsum,Tsum]))
                    #print np.sum(-pbatch*np.log2(pbatch)) + np.sum((Asum/totsum)*(condA)*np.log2(condA) + Csum/totsum*(condC)*np.log2(condC) + Gsum/totsum*(condG)*np.log2(condG) +  Tsum/totsum*condT*np.log2(condT))
                    ifoot[z,condbase*4 + p] =  np.sum(-pbatch*np.log2(pbatch)) + np.sum((Asum/totsum)*(condA)*np.log2(condA) + Csum/totsum*(condC)*np.log2(condC) + Gsum/totsum*(condG)*np.log2(condG) +  Tsum/totsum*condT*np.log2(condT)) - (len(sequences)-1)*3*np.log2(np.exp(1))/(2*len(seq_mat_temp[0,0,:]))
                    ifootrenorm[z,condbase*4 + p] = np.sum(-pbatch*np.log2(pbatch)) + .25*np.sum((condA)*np.log2(condA) + (condC)*np.log2(condC) + (condG)*np.log2(condG) +  condT*np.log2(condT)) - (len(sequences)-1)*3*np.log2(np.exp(1))/(2*len(seq_mat_temp[0,0,:]))
                except:
                    print z
                    pass
            #print 'Wildtype Seq is ' + wtseq

    np.save('/home/bill/Documents/energymatrix/mscS4815/infofootprints/ifootconditional' + expname + '.npy', ifoot)
    np.save('/home/bill/Documents/energymatrix/mscS4815/infofootprints/ifootrenormconditional' + expname + '.npy', ifootrenorm)




'''
fig = plt.figure()
ax = fig.add_subplot(111)

x = range(0, 0 + len(ifoot))
x[:] = [a - .45 for a in x]
numBins = 165
ax.bar(x,ifoot,color = 'b')
plt.xlabel('Base Position',fontsize = 20)
plt.ylabel('Mutual Information (bits)', fontsize = 20)
figtitle = 'Information Footprint ' + expname
plt.title(figtitle, fontsize = 22)   
plt.show()
plt.savefig(data_fnbase + expname + '.pdf')
plt.close()

fig = plt.figure()
ax = fig.add_subplot(111)

x = range(0, 0 + len(ifoot))
x[:] = [a - .45 for a in x]
numBins = 165
ax.bar(x,ifootrenorm,color = 'b')
plt.xlabel('Base Position',fontsize = 20)
plt.ylabel('Mutual Information (bits)', fontsize = 20)
figtitle = 'Information Footprint Renorm ' + expname
plt.title(figtitle, fontsize = 22)   
plt.show()
plt.savefig(data_fnbase + expname + 'renorm.pdf')
plt.close()
'''