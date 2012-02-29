# -*- coding: utf-8 -*-
from numpy import *
from pylab import *
import glob

def plot_pol():
  return None


def plot_sus(file):
  rho,E,tau = prop(file)
  data = genfromtxt(file)
  T = data[:,0]
  X = data[:,1]
  Xe = data[:,2]
  I = data[:,3]
  Ie = data[:,4]
  plot(T,X,'o-',label='E='+E+' $\\tau=$'+tau)
  #errorbar(T,X,yerr=Xe,label='E='+E+' $\\tau=$'+tau)

def prop(file):
  rho = file.find('_p')
  E = file.find('_E')
  tau = file.find('_t')
  dot = file.find('.d')

  rho = file[rho+2:E]
  E = file[E+2:tau]
  tau = file[tau+2:dot]

  return rho,E,tau

if __name__ == "__main__":
  files = glob.glob('sus*.dat')
  for file in files:
    plot_sus(file)


  xlabel('Temperature [$\\Delta J /k_B$]')
  ylabel('Electric susceptibility $\\chi$')
  title('Simulation Data')
  legend()
  show()