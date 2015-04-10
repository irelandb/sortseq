#!/home/wireland/bin/python2.7
import argparse, os
import pymc
import generic_energy_matrix_mscsvariedbin
import stepper
import sys
sys.path.append('/home/wireland/sortseq/')
sys.path.append('/home/wireland/')

# file to save trace to is passed as command line argument
# also, the config file should be passed as cli. This enables some 
# rudimentary error checking that everything is consistent.
parser = argparse.ArgumentParser()
parser.add_argument('runnum',help='Filename of database file to save traces to.')
parser.add_argument('savefn')
#parser.add_argument('cfg_fn',help='Filename pymc config file for the run.')

args = parser.parse_args()

# make sure the generic_energy_matrix module and this script are
# looking at the same data
# cli_cfg = os.path.abspath(args.cfg_fn)
# module_cfg = os.path.abspath(generic_energy_matrix.cfg_fn)
# if cli_cfg != module_cfg:
#     raise ValueError("generic_energy_matrix module %s and script cli %s are not consistent" % (module_cfg,cli_cfg))

#print args.db_fn
fulldbname = '/home/wireland/' + str(args.savefn) + str(args.runnum) + '.sql'
#fulldbname = '/home/wireland/lassoresults/MCMC/MCMCtest3_' + str(args.runnum) + '.sql'
M = pymc.MCMC(generic_energy_matrix_mscsvariedbin,db='sqlite',dbname=fulldbname)
M.use_step_method(stepper.GaugePreservingStepper,generic_energy_matrix_mscsvariedbin.emat)

M.sample(30000,thin=10)

