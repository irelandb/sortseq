#!/home/wireland/bin/python2.7
import argparse, os
import sys
sys.path.append('/home/wireland/sortseq/MCMC/')
import pymc
import ConfigParser
config = ConfigParser.RawConfigParser()
# everything that needs to be changed goes here ########

config.read('/home/wireland/mscS4-8-15/data/mscL.cfg') #location of config file, change this to run on different datasets
mut_region_start = config.getint('Input','mut_region_start') #not base pair number, distance from start of mut region
mut_region_length = config.getint('Input','mut_region_length')
expname = config.get('Input','expname') #ex MscS mut1, describes the experiment without the different batch numbers.

# file to save trace to is passed as command line argument
# also, the config file should be passed as cli. This enables some 
# rudimentary error checking that everything is consistent.
parser = argparse.ArgumentParser()
parser.add_argument('runnum',help='Filename of database file to save traces to.')
parser.add_argument('savefn')
parser.add_argument('condbase')
parser.add_argument('condident')
#parser.add_argument('cfg_fn',help='Filename pymc config file for the run.')

args = parser.parse_args()

# make sure the generic_energy_matrix module and this script are
# looking at the same data
# cli_cfg = os.path.abspath(args.cfg_fn)
# module_cfg = os.path.abspath(generic_energy_matrix.cfg_fn)
# if cli_cfg != module_cfg:
#     raise ValueError("generic_energy_matrix module %s and script cli %s are not consistent" % (module_cfg,cli_cfg))
f = open('/home/wireland/mscS4-8-15/runsdetails/condinfo.txt','w')
f.writelines([str(args.condbase) + '-' + str(args.condident) + '- = condbase + condident'])
f.close()
import generic_energy_matrix_conditional
import stepper
#print args.db_fn
fulldbname = '/home/wireland/mscS4-8-15/results/' + str(args.savefn) + str(args.runnum) + '.sql'
#fulldbname = '/home/wireland/lassoresults/MCMC/MCMCtest3_' + str(args.runnum) + '.sql'
M = pymc.MCMC(generic_energy_matrix_conditional,db='sqlite',dbname=fulldbname)
M.use_step_method(stepper.GaugePreservingStepper,generic_energy_matrix_conditional.emat)
f = open('/home/wireland/mscS4-8-15/runsdetails/' + str(args.savefn) + str(args.runnum) + '.txt','w')
f.writelines([str(args.condbase) + '-' + str(args.condident) + '- = condbase + condident \n', 'dbname = ' + fulldbname + '\n', 'mut_region_start = ' + str(mut_region_start) + '\n', 'mut_region_length = ' + str(mut_region_length) + '\n', 'exp_name = ' + str(expname)])
f.close()
M.sample(30000,thin=10)

