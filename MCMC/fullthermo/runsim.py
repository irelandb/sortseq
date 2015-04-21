# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 17:31:25 2015

@author: bill
"""

import argparse, os
import sys
sys.path.append(os.path.expanduser('~'))
sys.path.append(os.path.expanduser('~/sortseq/MCMC/fullthermo'))
sys.path.append(os.path.expanduser('~/sortseq/utils'))
import pymc
import ThermoSim
import stepper
# everything that needs to be changed goes here ########


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
fulldbname = '/home/bill/Documents/energymatrix/SimResults/' + str(args.savefn) + str(args.runnum) + '.sql'
#fulldbname = '/home/wireland/lassoresults/MCMC/MCMCtest3_' + str(args.runnum) + '.sql'
M = pymc.MCMC(ThermoSim,db='sqlite',dbname=fulldbname)
M.use_step_method(stepper.GaugePreservingStepper,ThermoSim.ematQ)
M.use_step_method(stepper.GaugePreservingStepper,ThermoSim.ematR)
M.use_step_method(pymc.Metropolis,ThermoSim.gamma)
M.use_step_method(pymc.Metropolis,ThermoSim.sR)
M.use_step_method(pymc.Metropolis,ThermoSim.sQ)
M.use_step_method(pymc.Metropolis.ThermoSim.R_0)
 
M.sample(30000,thin=10)