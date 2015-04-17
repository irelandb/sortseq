#!/opt/hpc/pkg/python-2.7/bin/python

import os,ConfigParser
import glob
import matplotlib.pyplot as plt
import pymc
import sys
sys.path.append('/home/bill/Documents/energymatrix/mscsstuff/')
sys.path.append('/home/bill/Documents/energymatrix/lacstardata/')
sys.path.append('/home/bill/Documents/energymatrix/31515mscs/')
import scipy as sp
import scipy.ndimage, scipy.stats
import MCMC_utils
import MCMC_utils_mscs


config = ConfigParser.RawConfigParser()
# stuff to be changed goes here
#cfg_fn = '/data/kinney/data/illumina_sequencing/11.09.25_ecoli/pymc/Sample_mscL1/runs/mscL_MG1655_site_53_72/site_53_72.cfg'
cfg_fn = '/home/wireland/mscS4-8-15/data/mscL.cfg'
# end stuff to be changedm
config.read(cfg_fn)
n_runs = config.getint('Output','num_runs')
runs_dir = config.get('Output','output_dir')
print runs_dir
mut_region_start = config.getint('Input','mut_region_start')
mut_region_length = config.getint('Input','mut_region_length')
data_fnbase = config.get('Input','data_fnbase')
output_basefn = config.get('Output','output_basename')
seq_start = config.getint('Input','seqstart')
seq_end = config.getint('Input','seqend')

# name of database file containing analysis run
# runs_dir = '/bluearc/home/kinney/djones/Code/pymc'
# fn = '/data/kinney/data/illumina_sequencing/12.08.17_lac/12_inducibility_shared/pymc_runs/x61_batches34.sql'
# data_dir = '/data/kinney/data/simulations/sortseq/sim_002'
# data_fn = '/data/kinney/data/illumina_sequencing/12.08.17_lac/BaseCalls/processed_by_script.py/pymc/x61_seqsbins34.csv'
# data_fn = '/data/kinney/data/simulations/sortseq/sim_002/sim_sortseq_data_002.txt'
# mut_region_start = 0
# mut_region_length = 14


# now loop through all the sql files
plt.figure()
for fn in glob.glob(os.path.join(runs_dir,'*.sql')):
    print fn
    
    # make a new directory corresponding to the run where we'll save everything
    run_name = os.path.basename(fn).split('.')[0]
    output_dir = os.path.join(os.path.dirname(fn),run_name)
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    # data_name = '_'.join(run_name.split('_')[:2]) + '.csv'
    # data_name = data_name.replace('batches','seqsbins')
    # print data_name
    #data_fn = os.path.join(data_dir,data_name)


    # load the database and save and plot the mean value of the emat
    burn_in = 1000
    db = pymc.database.sqlite.load(fn)
    # emat_mean = db.emat.stats()['mean']
    emat_mean = sp.mean(db.trace('emat')[burn_in:],axis=0)
    
    # change the sign of emat_mean if necessary. We want a negative
    # correlation between energy and batch number (because the higher
    # the energy, the lower the expression)
    seq_mat, batch_vec = MCMC_utils_mscs.load_unique_seqs_batches(data_fnbase,4,seq_start,seq_end,mut_region_start,mut_region_length)
    energies = sp.zeros(len(batch_vec))
    for i in range(len(batch_vec)):
        energies[i] = sp.sum(seq_mat[:,:,i]*emat_mean)
    r = scipy.stats.pearsonr(energies,batch_vec)[0]
    if r>0:
        emat_mean = MCMC_utils.fix_matrix_gauge(-emat_mean)
    else:
        emat_mean = MCMC_utils.fix_matrix_gauge(emat_mean)
    sp.savetxt(os.path.join(output_dir,run_name+'_emat_mean.txt'),emat_mean)

    # compute mutual information and joint pdf
    MI,f_reg = MCMC_utils.compute_MI(seq_mat,batch_vec,emat_mean)
    MI_f = open(os.path.join(output_dir,run_name+'_MI.txt'),'w')
    MI_f.write(str(MI))

    print MI

    # make heat map of mean matrix
    plt.clf()
    plt.imshow(MCMC_utils.zero_matrix(emat_mean),interpolation='nearest')
    plt.title('Energy Matrix, ' + run_name + ', MI: %.5f' % MI)
    plt.savefig(os.path.join(output_dir,run_name+'.png'))


    plt.imshow(f_reg,interpolation='nearest',aspect='auto')
    plt.title('Joint regularized pdf, ' + run_name + ', MI: %.5f' % MI)
    plt.xlabel('Rank order')
    plt.ylabel('Batch number')
    plt.savefig(os.path.join(output_dir,run_name+'_regpdf.png'))

    # now, plot mutual information over time
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
    plt.title('MI, ' + run_name + ' MI: %.5f' % MI )
    plt.savefig(os.path.join(output_dir,run_name+'_MItrace.png'))
