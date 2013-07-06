# -*- coding: utf-8 -*-
"""
@author: Óscar Nájera

Funtions for labeling files & plots
"""
import re

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