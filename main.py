import argparse

parser = argparse.ArgumentParser(description="Relaxor Simulator")

parser.add_argument('-L', '--length', metavar='L',
		    type=int, nargs=1, default=12,
		    help='Side length of simulation box')
parser.add_argument('-n', '--numexps', metavar='N',
		    type=int, nargs=1, default=4,
		    help='Number of sampling experiments')
parser.add_argument('-p', '--rho', type=float, nargs='*',
		    default=0, help='Mean Ferroelectricity')
parser.add_argument('-T', '--temperatures', metavar='T', type=float,
		    nargs=3, default=[7,0.2,1],
		    help='Temperature range of experiment. Corresponds to [Ti dT Tf]')
parser.add_argument('-E', '--fields', metavar='E0', type=float, nargs='*',
		    default=0.1, help='External Electric Field amplitude')
parser.add_argument('-t', '--tau', metavar='t', type=float, nargs='*',
		    default=100, help='External Electric Field period')