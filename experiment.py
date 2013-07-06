# -*- coding: utf-8 -*-
"""
@author: Óscar Nájera

Functions to perform simulations of relaxors
"""
from __future__ import division
import numpy as np
import libmat as lm
import copy
from os import path, remove
from multiprocessing import Pool, cpu_count
from labels import get_label, gen_label

def build_temp(experiment):
    """Generate numpy array for temperatures"""
    temp = np.arange( experiment['Ti'], experiment['Tf'], experiment['dT'] )
    return temp

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
