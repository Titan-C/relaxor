# -*- coding: utf-8 -*-
"""
@author: Óscar Nájera

Functions to perform experiments and data analysis for simulated relaxors
"""
from __future__ import division
import numpy as np
import libmat as lm
import re
import copy
from pylab import plot, errorbar
from scipy.integrate import simps
from os import path, remove
from multiprocessing import Pool, cpu_count

def steps_calibrator(setup):
    """Stablishes proper amount of MCS during simulation
       edits setup for later labeling and returns steps interval"""
    periods = setup.Eiter // setup.tau
    if periods > setup.osc:
        setup.Eiter = periods * setup.tau
    else:
        setup.Eiter = setup.osc * setup.tau

    if setup.Qiter > setup.Eiter:
        setup.Qiter = setup.Eiter

    return np.arange(0, 2*np.pi*periods, 2*np.pi / setup.tau)

def gen_label(setup):
    """Generate string to label simulation files"""

    label  = setup.name + '_N' + str(setup.L**3)
    label += '_p'+str(setup.rho) +'_E'+str(setup.field) +'_t'+str(setup.tau)
    label += '_Ti'+str(setup.Ti) +'_Tf'+str(setup.Tf) +'_dT'+str(setup.dT)
    label += '_X' + str(setup.Eiter) + '_Q' + str(setup.Qiter)
    return label

def get_label(filename):
    """Extract from filename simulation setup"""
    pat  = r'_N(?P<N>\d+)'
    pat += r'_p(?P<rho>[\d\.]+)_E(?P<E>[\d\.]+)_t(?P<tau>\d+)'
    pat += r'_Ti(?P<Ti>[\d\.]+)_Tf(?P<Tf>[\d\.]+)_dT(?P<dT>[-+\d\.]+)'
    pat += r'_X(?P<Eiter>\d+)_Q(?P<Qiter>\d+)'
    experiment = re.search(pat, filename).groupdict()
    for key in ['N', 'Eiter', 'Qiter']:
        experiment[key] = int(experiment[key])
    for key in ['rho', 'E', 'tau', 'Ti', 'Tf', 'dT']:
        experiment[key] = float(experiment[key])

    return experiment

def build_temp(experiment):
    """Generate numpy array for temperatures"""
    temp = np.arange( experiment['Ti'], experiment['Tf'], experiment['dT'] )
    return temp

def eval_frozen(filename, frozen = None):
    """Evaluate frozen dipoles proportion for determinate
       experiment(filename) given the list of thresholds
       to consider a dipole frozen
       eval_frozen(filename,[1,0.95])"""
    if frozen is None:
        frozen = [1, 0.9]
    experiment = get_label(filename)
    thermostat = build_temp(experiment)
    with open(filename, 'rb') as datf:
        sigma = np.fromfile(datf, dtype=np.int32) / experiment['Eiter']
    sig = np.abs(sigma).reshape(-1, thermostat.size, experiment['N'])

    slow_dipoles = []
    for bound in frozen:
        slow_dipoles.append( np.sum(sig>=bound, axis=2) / experiment['N'] )

    for i in xrange(len(slow_dipoles)):
        errorbar( thermostat, slow_dipoles[i].mean(axis=0),
                 slow_dipoles[i].std(axis=0) )

    return slow_dipoles, thermostat

def get_pol_log(filename):
    """Get polarization data and experiment setup from file"""
    experiment = get_label(filename)
    thermostat =  build_temp(experiment)
    with open(filename, 'rb') as datf:
        pol_log = np.fromfile(datf, dtype=np.double)
        pol_log = pol_log.reshape( -1, thermostat.size, experiment['Eiter'] )
    return pol_log, thermostat, experiment

def eval_pol(filename):
    """Calculate polarization"""
    pol_log, thermostat = get_pol_log(filename)[:2]
    pol  = pol_log.mean(axis=2)
    pol_err = pol_log.std(axis=2)

    plot(thermostat, pol.mean(axis=0), 'k*-',
         thermostat, np.abs(pol).mean(axis=0))

    return pol, pol_err, thermostat

def eval_sus(filename):
    """Calculate dielectric susceptibility according to Liu paper"""
    pol_log, thermostat, experiment = get_pol_log(filename)

    periods = experiment['Eiter'] // experiment['tau']
    steps = np.arange(0, 2*np.pi*periods, 2*np.pi / experiment['tau'] )
    cos_wave = np.cos(steps)
    sin_wave = np.sin(steps)

    sus_re = simps(cos_wave*pol_log)/experiment['Eiter']/experiment['E']
    sus_im = simps(sin_wave*pol_log)/experiment['Eiter']/experiment['E']

    plot(thermostat, sus_re.mean(axis=0), thermostat, sus_im.mean(axis=0))

    return sus_re, sus_im, thermostat

def susGUI(filename, frozen = None):
    """Calculate dielectric susceptibility X=(1-p)/T
       where p:proportion of slow dipoles & T:temperature"""
    if frozen is None:
        frozen = [0.9, 0.8, 0.6]
    slow_dipoles, thermostat = eval_frozen(filename, frozen)

    for i in range(len(slow_dipoles)):
        plot(thermostat,(1-slow_dipoles[i].mean(axis=0))/thermostat,'*')

def required_simulations(numexps, experiment_label):
    """Calculate the number of needed remaining simulations in order
       to reach the required number"""
    filename = 'log_pol_' + experiment_label + '.dat'

    try:
        size = path.getsize(filename)
    except OSError:
        return numexps

    experiment = get_label(filename)
    thermostat = build_temp(experiment)
    count = size/8/experiment['Eiter']/thermostat.size
    if int(count) == count:
        if numexps - count > 0:
            return numexps - int(count)
        else:
            return 0
    else:
        remove(filename)
        return numexps

def setup_experiments(instructions):
    """Arrange multiple experiments accordind to instructions"""
    experiment = copy.deepcopy(instructions)
    experiment_list = []

    for rho in instructions.rho:
        for field in instructions.field:
            for tau in instructions.tau:
                experiment.rho = rho
                experiment.field = field
                experiment.tau = tau
                experiment_list.append(copy.deepcopy(experiment))

    pool = Pool(cpu_count())
    return pool.map(do_experiment, experiment_list)

def do_experiment(setup):
    """Perform experiment accoding to setup"""
    exp_label = gen_label(setup)

    num = required_simulations(setup.numexps, exp_label)
    if num > 0:
        relaxor = lm.Material(setup.L, setup.rho, exp_label,
                            setup.pol, setup.logH, setup.logS)
        steps = steps_calibrator(setup)
        field = lm.double_vector()
        field[:] = setup.field * np.cos(steps)
        temperatures = lm.double_vector()
        temperatures[:] = build_temp(vars(setup))
        print exp_label, relaxor.oven(num, setup.Qiter, temperatures,
                           field, setup.pol)
