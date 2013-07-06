# -*- coding: utf-8 -*-
"""
@author: Óscar Nájera

Simulations of reference papers
"""
import main
import numpy as np
import experiment as ee
from pylab import plot,figure
#Gui, H., Gu, B., & Zhang, X. (1995). Dynamics of the freezing process in
#relaxor ferroelectrics.pdf. Physical Review B, 52(5), 3135.
def susGUI(filename, frozen = None):
    """Calculate dielectric susceptibility X=(1-p)/T
       where p:proportion of slow dipoles & T:temperature"""
    if frozen is None:
        frozen = [0.9, 0.8, 0.6]
    slow_dipoles, thermostat = eval_frozen(filename, frozen)

    for i in range(len(slow_dipoles)):
        plot(thermostat,(1-slow_dipoles[i].mean(axis=0))/thermostat,'*')

#Gui95 setup experiment 1, figs 2 & 3
setup = main.parser.parse_args('-L 16 -p 0. -E 0. -t 1 -Ti 2.5 -Tf 0 -dT -0.1 -Eiter 400 --logH --logS -name Gui95fix'.split())
ee.setup_experiments(setup)
label = 'Gui95fix_N4096_p0.0_E0.0_t1_Ti2.5_Tf0.0_dT-0.1_X400_Q200.dat'
with open('log_H_'+label, 'rb') as f:
    H = np.fromfile(f, dtype=np.double)
figure()
plot(H)
figure()
ee.eval_frozen('log_sigma_'+label)
figure()
ee.susGUI('log_sigma_'+label)

##Gui95 setup experiment 2, figs 4 & 5
setup = main.parser.parse_args('-L 16 -p 0. -E 0. 0.8 2.0 10. -t 1 -Ti 2.5 -Tf 0 -dT -0.1 -Eiter 400 --logS -name Gui95fix'.split())
ee.setup_experiments(setup)

#Lui00

