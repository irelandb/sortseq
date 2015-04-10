#!/usr/bin/python
import csv, gzip, subprocess, sys, glob
from Bio import SeqIO, Seq, AlignIO
from Bio.Align.Applications import MuscleCommandline
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import scipy as sp
print 'finished import'
noalign = 0
# These will be replaced with a config file
read_length = 100
# fn_raw_F = ['lane1_Undetermined_L001_R1_001.fastq.gz']
# fn_raw_R = ['lane1_Undetermined_L001_R2_001.fastq.gz']
fn_raw_F = glob.glob('lane1_Undetermined_*.fastq.gz')

fn_raw_F.sort()

qualfilt_fn = 'all.qcfilt.fastq.gz'
qualfilt_fh = gzip.open(qualfilt_fn,'w')
print 'Created File'
# First, filter the sequences on the basis of quality
# Filter on both minimum and median quality

min_qual = 10
median_qual = 20
n_ok_qual = 0

print 'about to start loop'
for i in range(len(fn_raw_F)):
    print 'Processing ', fn_raw_F[i]
    forward_iterator = SeqIO.parse(gzip.open(fn_raw_F[i]),'fastq')
    try:
        while True:
            forward_record = forward_iterator.next()
            all_quals = (forward_record.letter_annotations['phred_quality'])
            min_qual_all = min(all_quals)
            median_qual_all = sp.median(all_quals)
            if min_qual_all >= min_qual and median_qual_all >= median_qual:
                SeqIO.write(forward_record,qualfilt_fh,'fastq')
                # should really check that the id's are the same....
                n_ok_qual = n_ok_qual+1
    except StopIteration:
        pass

# close the quality filter file handles
qualfilt_fh.close()
print 'Number passing quality filter: ', n_ok_qual
'''
seq_iterator = SeqIO.parse(gzip.open('all.qcfilt.fastq.gz'),'fastq')
     wtseqs = []
    ...: dstpaseqs = []
    ...: print 'Beginning Barcode Filtering'
    ...: try:
    ...:     while True:
    ...:         a = seq_iterator.next()
    ...:         if str(a.seq[0:4]) == 'CAAA':
    ...:             wtseqs.append(str(a.seq))
    ...:         if str(a.seq[0:4]) == 'CACC':
    ...:             dstpaseqs.append(str(a.seq))
    ...: except StopIteration:
    ...:     pass
'''
# next, split everything up according to bar code. we'll still keep
# forward and reverse separate at this point
barcode_fn = '/home/bill/Documents/DiversityData/April8/barcodes.csv'
barcode_dict = {} # instantiate barcode dictionary
bcfiles_dict = {}
reverse_dict = {}
wt_seq_dict = {}
csvfile = open(barcode_fn,'r')
reader = csv.DictReader(csvfile)
for row in reader:
    barcode_dict[row['fwd_barcode']] = row['experiment_name'] + '_B' + row['batch_number']
    reverse_dict[row['fwd_barcode']] = row['rev_barcode']
    bcfiles_dict[row['experiment_name'] + '_B' + row['batch_number']] = gzip.open(row['experiment_name'] + '_B' + row['batch_number']+'.qcfilt.fastq.gz','w')
    wt_seq_dict[row['experiment_name'] + '_B' + row['batch_number']] = row['wt_seq']
print barcode_dict
print bcfiles_dict
n = 0
seq_iterator = SeqIO.parse(gzip.open(qualfilt_fn),'fastq')
print 'Beginning barcode filtering'
try:
    while True:
        forward_record = seq_iterator.next()
        for i, code in enumerate(barcode_dict.keys()):
            cat_barcode = str(forward_record.seq[0:len(code)]) + str(forward_record.seq[len(code) + 60:len(code) + 60 + len(reverse_dict[code])])
            if cat_barcode in code + reverse_dict[code]:
                #print cat_barcodes
                SeqIO.write(forward_record,bcfiles_dict[barcode_dict[code]],'fastq')
except StopIteration:
    pass
# close the files now that they're done
for bcfile in bcfiles_dict.values():
    bcfile.close()
for batch in bcfiles_dict.keys():
    print batch
    # open the appropriate files
    bcfilt_file = gzip.open(batch + '.qcfilt.fastq.gz', 'r')
    aln_file = open(batch + '.qcfilt.aln.fasta', 'w')
    seq_iterator = SeqIO.parse(bcfilt_file,'fastq')
    try:
        while True:
            forward_record = seq_iterator.next()
            SeqIO.write(forward_record,aln_file,'fasta')
    except StopIteration:
        pass
            
