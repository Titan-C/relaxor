# -*- coding: utf-8 -*-
from pylab import *
from sys import argv
from glob import glob
from utils import zoom_effect02

def HlogPlot(file,stl='+'):
  H=genfromtxt(file)
  plot(H,stl,label=legendSet(file))

def polPlot(file,stl='-o', error=False):
  lb=legendSet(file)
  data = genfromtxt(file)
  T = data[:,0]
  p = data[:,1]
  if not error:
    plot(T,p,stl,label=u'Polarización')
  else:
    pe = data[:,2]
    errorbar(T,p,yerr=pe,fmt=stl,label=lb)

def susPlot(file,stl='-o', error=False, IM=False):
  lb=legendSet(file)
  data = genfromtxt(file)
  T = data[:,0]
  X = data[:,1]
  if not error:
    plot(T,X,stl,label='$\\chi$ $\\mathcal{R}eal$')
  else:
    std = data[:,3]
    errorbar(T,X,yerr=std,fmt=stl, label='$\\mathcal{R}eal$')
  if IM:
    Xi = data[:,4]
    plot(T,Xi,stl,label='$\\chi$ $\\mathcal{I}maginaria$')

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

def susLabel():
  xlabel('Temperatura [$\\Delta J /k_B$]')
  ylabel(u'Susceptibilidad dieléctrica $\\chi$')
  title(u'Susceptibilidad en un proceso de enfriamiento\n $\\rho=0.8$, $E_0=0.4$ y $\\tau=1500$')
  
def dieLabel():
  xlabel('Temperatura [$^{\circ}K$]')
  ylabel(u'Constante dieléctrica $\\varepsilon_r$')
  title(u'Curvas de ajuste a los datos experimentales')

def polLabel():
  xlabel('Temperatura [$\\Delta J /k_B$]')
  ylabel(u'Polarización normada del sistema $[\\overline{\\mu}/N]$')
  title(u'Polarización global en un proceso de enfriamiento\nCampo externo nulo')

def HlogLabel():
  xlabel('Iteraciones $[MCS/dipolo]$')
  ylabel(u'Energía del sistema $[\Delta J]$')
  suptitle(u'Evolución de la energía en un proceso de enfriamiento\nCampo externo nulo')

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
    return '$E=$'+E
    
  else:
    mat = file.find('P')
    return file[mat:mat+5]+' '+file[mat+5:-4]

def filesusPlot(path,stl='-o',error=False,IM=False):
  files = sort(glob(path))
  for file in files:
    susPlot(file,stl,error,IM)
  susLabel()
  legend()
  show()

def filepolPlot(path,stl='-o',error=False):
  files = sort(glob(path))
  for file in files:
  #  if (float(simIdentifier(file)[0]) > 0.4):
      polPlot(file,stl,error)
  polLabel()
  legend()
  show()

def fileHlogPlot(path,stl='+',error=False):
  '''Grafica el historial de la energía del sistema y hace
     un zoom a un intervalo predefinido'''
  files = sort(glob(path))
  ax1 = plt.subplot(211)
  ax1.set_xlim(5.8*3350, 6.3*3350)
  ax1.set_ylim(-11200, -8400)
  ax1.annotate(u'$Equilibración$', xy=(6*3350+50,-9000),xytext=(6*3350+50, -9000))
  plt.axvspan(6*3350, 6*3350+350, facecolor='r', alpha=0.25)
  for file in files:
    HlogPlot(file,stl)
  ax2 = plt.subplot(212)
  for file in files:
    HlogPlot(file,stl)
  zoom_effect02(ax1, ax2)
  HlogLabel()
  legend()
  show()

if __name__ == "__main__":
  filesusPlot(argv[1])
