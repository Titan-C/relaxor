# -*- coding: utf-8 -*-
"""
@author: Óscar Nájera

Argparser for relaxor simulations
"""
import argparse
import experiment as ee

import copy

parser = argparse.ArgumentParser(description="Relaxor Simulator")

parser.add_argument('-L', metavar='L', type=int, default=12,
                    help='Side length of simulation box')
parser.add_argument('-n', '--numexps', metavar='N', type=int, default=4,
                    help='Number of sampling experiments')
parser.add_argument('-p', '--rho', type=float, nargs='*',
                    default=[0], help='Mean Ferroelectricity')
parser.add_argument('-E', '--field', metavar='E0', type=float, nargs='*',
		    default=[0.1], help='External Electric Field amplitude')
parser.add_argument('-t', '--tau', metavar='t', type=float, nargs='*',
		    default=[100], help='External Electric Field period')
parser.add_argument('-Ti', metavar='T', type=float, default=7.15,
		    help='Starting temperature of experiment')
parser.add_argument('-Tf', metavar='T', type=float, default=0,
		    help='Finishing temperature of experiment')
parser.add_argument('-dT', metavar='T', type=float, default=-0.5,
		    help='Temperature step during experiment')
parser.add_argument('-Eiter', metavar='N', type=int, default=3000,
		    help='Time steps for Experiment')
parser.add_argument('-osc', metavar='N', type=int, default=2,
            help='Minimum count of field oscilation during Experiment')
parser.add_argument('-Qiter', metavar='N', type=int, default=200,
		    help='Time steps for Equilibration')
parser.add_argument('-name', type=str, default='',
		    help='Additional experiment label')
parser.add_argument('--pol', action='store_true', default=False,
                    help='Initiate System in a fully polarized state')
parser.add_argument('--logH', action='store_true', default=False,
                    help='Keep a log of systems total energy')
parser.add_argument('--logS', action='store_true', default=False,
                    help='Keep a log of systems individual dipoles state')



if __name__ == "__main__":
    setup = parser.parse_args()
    experiment = copy.deepcopy(setup)

    for p in setup.rho:
        for E in setup.field:
            for t in setup.tau:
                experiment.rho = p
                experiment.field = E
                experiment.tau = t
                ee.do_experiment(experiment)
