# -*- coding: utf-8 -*-
"""
@author: Óscar Nájera

Simulations of reference papers
"""
import main
import numpy as np
import experiment as ee
import analysis as an
import plots as pt
from pylab import plot, figure, show
#Gui, H., Gu, B., & Zhang, X. (1995). Dynamics of the freezing process in
#relaxor ferroelectrics.pdf. Physical Review B, 52(5), 3135.
def susGUI(filename, frozen = None):
    """Calculate dielectric susceptibility X=(1-p)/T
       where p:proportion of slow dipoles & T:temperature"""
    if frozen is None:
        frozen = [0.9, 0.8, 0.6]
    slow_dipoles, thermostat, experiment = an.eval_frozen(filename, frozen)

    for i in range(len(slow_dipoles)):
        plot(thermostat,(1-slow_dipoles[i].mean(axis=0))/thermostat,'*')

#Gui95 setup experiment 1, figs 2 & 3
setup = main.parser.parse_args('-L 16 -p 0. -E 0. -t 1 -Ti 2.5 -Tf 0 -dT -0.1 -Eiter 400 --logH --logS -name Gui95fix'.split())
ee.setup_experiments(setup)
label = 'Gui95fix_N4096_p0.0_E0.0_t1_Ti2.5_Tf0.0_dT-0.1_X400_Q200.dat'
figure()
pt.h_log_plot('log_H_'+label)
show()
figure()
susGUI('log_sigma_'+label)
show()
##Gui95 setup experiment 2, figs 4 & 5
setup = main.parser.parse_args('-L 16 -p 0. -E 0. 0.8 2.0 10. -t 1 -Ti 2.5 -Tf 0 -dT -0.1 -Eiter 400 --logS -name Gui95fix'.split())
ee.setup_experiments(setup)

#Lui00

