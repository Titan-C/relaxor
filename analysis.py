# -*- coding: utf-8 -*-
"""
@author: Óscar Nájera

Funtions for analyzing saved data from simulations
"""
from __future__ import division
import numpy as np
from experiment import build_temp
from labels import get_label
from scipy.integrate import simps

def get_sigmas_log(filename):
    """Get individual dipole history log and experiment setup from file"""
    experiment = get_label(filename)
    thermostat = build_temp(experiment)
    with open(filename, 'rb') as datf:
        sigma = np.fromfile(datf, dtype=np.int32) / experiment['Eiter']
        sigma = sigma.reshape(-1, thermostat.size, experiment['N'])

    return sigma, thermostat, experiment

def get_pol_log(filename):
    """Get polarization data and experiment setup from file"""
    experiment = get_label(filename)
    thermostat =  build_temp(experiment)
    with open(filename, 'rb') as datf:
        pol_log = np.fromfile(datf, dtype=np.double)
        pol_log = pol_log.reshape( -1, thermostat.size, experiment['Eiter'] )
    return pol_log, thermostat, experiment

def eval_frozen(filename, frozen = None):
    """Evaluate frozen dipoles proportion for determinate
       experiment(filename) given the list of thresholds
       to consider a dipole frozen
       eval_frozen(filename,[1,0.95])"""
    sigma, thermostat, experiment = get_sigmas_log(filename)
    sigma_abs = np.abs(sigma)

    slow = []
    if frozen is None:
        frozen = [1, 0.9]
    for bound in frozen:
        slow.append(np.sum(sigma_abs>=bound, axis=2) / experiment['N'])

    return slow, thermostat, experiment

def eval_pol(filename):
    """Calculate polarization"""
    pol_log, thermostat, experiment = get_pol_log(filename)
    pol  = pol_log.mean(axis=2)
    pol_err = pol_log.std(axis=2)

    return pol, pol_err, thermostat, experiment

def eval_sus(filename):
    """Calculate dielectric susceptibility according to Liu paper"""
    pol_log, thermostat, experiment = get_pol_log(filename)

    periods = experiment['Eiter'] // experiment['tau']
    steps = np.arange(0, 2*np.pi*periods, 2*np.pi / experiment['tau'] )
    cos_wave = np.cos(steps)
    sin_wave = np.sin(steps)

    sus_re = simps(cos_wave*pol_log)/experiment['Eiter']/experiment['E']
    sus_im = simps(sin_wave*pol_log)/experiment['Eiter']/experiment['E']

    return sus_re, sus_im, thermostat, experiment
