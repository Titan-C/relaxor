# -*- coding: utf-8 -*-
"""
@author: Óscar Nájera

Plotting scripts
"""
import numpy as np
from analysis import get_sigmas_log
from pylab import hist

def sigmahist(filename):
    """Plots histogram of dipoles permanence i a certain orientation"""
    sigma, thermostat, setup = get_sigmas_log(filename)
    sigma = np.rollaxis(sigma, 1).reshape(thermostat.size, -1)
    for i in range(thermostat.size):
        hist(sigma[i], 64 , normed=True)

