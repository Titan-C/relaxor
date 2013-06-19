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
parser.add_argument('-E', '--field', metavar='E0', type=float, nargs='*',
		    default=0.1, help='External Electric Field amplitude')
parser.add_argument('-t', '--tau', metavar='t', type=float, nargs='*',
		    default=100, help='External Electric Field period')
parser.add_argument('-Ti', metavar='T', type=float,
		    nargs=1, default=7.15,
		    help='Starting temperature of experiment')
parser.add_argument('-Tf', metavar='T', type=float,
		    nargs=1, default=0.1,
		    help='Finishing temperature of experiment')
parser.add_argument('-dT', metavar='T', type=float,
		    nargs=1, default=-0.5,
		    help='Temperature step during experiment')
parser.add_argument('--Nexp', metavar='N', type=int, nargs=2,
		    default=[3000,100], help='Time steps for Experiment and Equilibration')


if __name__ == "__main__":
  args=parser.parse_args()
  #main(args)