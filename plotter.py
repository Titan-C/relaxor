# -*- coding: utf-8 -*-
from pylab import *
from sys import argv
from glob import glob

def simpolPlot(file):
  rho,E,tau = identifier(file)
  data = genfromtxt(file)
  T = data[:,0]
  p = data[:,1]
  pe = data[:,2]
  errorbar(T,p,yerr=pe,label='E='+E+' $\\tau=$'+tau)
  return rho,E,tau

def simsusPlot(file):
  rho,E,tau = identifier(file)
  data = genfromtxt(file)
  T = data[:,0]
  X = data[:,1]
  Xe = data[:,2]
  I = data[:,3]
  Ie = data[:,4]
  #plot(T,X,'o-',label='E='+E+' $\\tau=$'+tau)
  errorbar(T,X,yerr=Xe,label='E='+E+' $\\tau=$'+tau)
  return rho,E,tau
  
def expsusPlot(file):
  data = genfromtxt(file)
  T = data[:,0]
  X = data[:,1]
  Xe = data[:,2]*X
  plot(T,X,T,Xe)

def fittedPlot(equation):
  T = equation.dataCache.allDataCacheDictionary['IndependentData'][0]
  E = equation.dataCache.allDataCacheDictionary['DependentData']
  a,b,c,d = equation.solvedCoefficients
  x=arange(T.min()*0.8,T.max()*1.1,T.max()/150)
  F=a*(x-b)/(x*x+c*x+d)
  plot(T,E,'o')
  plot(x,F)

def susLabel():
  xlabel('Temperature [$\\Delta J /k_B$]')
  ylabel('Electric susceptibility $\\chi$')
  title('Simulation Data & fits')
  legend(('DATA','FIT'))
  show()


def identifier(file):
  rho = file.find('_p')
  E = file.find('_E')
  tau = file.find('_t')
  dot = file.find('.d')

  rho = file[rho+2:E]
  E = file[E+2:tau]
  tau = file[tau+2:dot]

  return rho,E,tau

def generator(path):
  files = glob(path)
  for file in files:
    rho = simsusPlot(file)[0]
  xlabel('Temperature [$\\Delta J /k_B$]')
  ylabel('Electric susceptibility $\\chi$')
  title('Simulation Data for $\\rho=$'+rho)
  legend()
  show()

if __name__ == "__main__":
  generator('sus*'+argv[1]+'*.dat')