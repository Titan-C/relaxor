# -*- coding: utf-8 -*-
from pylab import *
from sys import argv
from glob import glob

def elogPlot(file,fmt='+'):
  H=genfromtxt(file)
  plot(H,fmt)
  xlabel('Iteraciones $[MCS/dipolo]$')
  ylabel(u'Energía del  sistema $[\Delta J]$')
  title(u'Evolución de la energía en un proceso de enfriamiento')
  show()

def simpolPlot(file):
  lb=legendSet(file)
  data = genfromtxt(file)
  T = data[:,0]
  p = data[:,1]
  pe = data[:,2]
  errorbar(T,p,yerr=pe,label=lb)
  return rho,E,tau

def susPlot(file,fmt='-o', error=False):
  lb=legendSet(file)
  data = genfromtxt(file)
  T = data[:,0]
  X = data[:,1]
  if not error:
    plot(T,X,fmt,label=lb)
  else:
    std = data[:,3]
    errorbar(T,X,yerr=std,label=lb)

def fittedPlot(equation,file):
  T = equation.dataCache.allDataCacheDictionary['IndependentData'][0]
  E = equation.dataCache.allDataCacheDictionary['DependentData']
  a,b,c,d = equation.solvedCoefficients
  x=arange(T.min()*0.8,T.max()*1.1,T.max()/150)
  F=a*(x-b)/(x*x+c*x+d)
  plot(T,E,'o',x,F,label=legendSet(file))

def scalefittedPlot(equation,coefs,file):
  T = equation.dataCache.allDataCacheDictionary['IndependentData'][0]
  E = equation.dataCache.allDataCacheDictionary['DependentData']
  a,b,c,d = [float(x) for x in coefs]
  j,u = equation.solvedCoefficients
  x=arange(T.min()*0.8,T.max()*1.1,T.max()/150)
  F=u/j*a*(j*x-b)/(j**2*x*x+c*j*x+d)
  plot(T,E,'o',x,F,label=legendSet(file))

def axisLabel():
  xlabel('Temperatura [$\\Delta J /k_B$]')
  ylabel('Electric susceptibility $\\chi$')
  title('Data & fits')

def simIdentifier(file):
  rho = file.find('_p')
  E = file.find('_E')
  tau = file.find('_t')
  L = file.find('_L')

  rho = file[rho+2:E]
  E = file[E+2:tau]
  tau = file[tau+2:L]

  return rho,E,tau

def legendSet(file):
  if file.find('_p') > 0:
    rho,E,tau = simIdentifier(file)
    return '$\\rho=$'+rho+' E='+E+' $\\tau=$'+tau
    
  else:
    mat = file.find('P')
    return file[mat:mat+5]+' '+file[mat+5:-4]

def filesusPlot(path,fmt='-o',error=False):
  files = sort(glob(path))
  for file in files:
    susPlot(file,fmt,error)
  axisLabel()
  legend()
  show()

if __name__ == "__main__":
  filesusPlot(argv[1])
