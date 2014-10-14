# -*- coding: utf-8 -*-
"""
@author: Óscar Nájera

Plotting scripts
"""
import numpy as np
from analysis import get_sigmas_log
from labels import get_label, labeler
from pylab import hist, plot, xlabel, ylabel, title

def sigmahist(filename):
    """Plots histogram of dipoles permanence in a certain orientation"""
    sigma, thermostat, setup = get_sigmas_log(filename)
    sigma = np.rollaxis(sigma, 1).reshape(thermostat.size, -1)
    for i in range(thermostat.size):
        hist(sigma[i], 64 , normed=True)
    return sigma, thermostat, setup

def h_log_plot(filename):
    """Plots the energy evolution of the experiment"""
    experiment = get_label(filename)
    with open(filename, 'rb') as datf:
        hlog = np.fromfile(datf, dtype=np.double)
        
    plot(hlog,'+')
    xlabel('Iteration Count $[MCS/dipole]$')
    ylabel(r'System Energy $[\Delta J]$')
    title('Energy evolution ' + labeler(experiment) )