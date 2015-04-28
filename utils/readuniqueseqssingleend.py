#The purpose of this script is to read in collated full libraries and output seq_mat and a copy number of the sequences
import csv, gzip, subprocess, sys, glob
from Bio import SeqIO, Seq, AlignIO
from Bio.Align.Applications import MuscleCommandline
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def collatedmat(filename):
    f = open(filename,'r')
    seqs = []
    for record in SeqIO.parse(f,'fasta'):
         seqs.append(str(record.seq))

    #There seems to alwasy be a trailing empty selection so...
   
           
    
    return seqs
    