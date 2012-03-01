
'''Primera función que ajusta un archivo al mi función tomada'''
import os, sys, inspect

# ensure pyeq2 can be imported
import pyeq2
from numpy import *
from pylab import *
import glob

def data_fitter(file):
  f=open(file,'r')
  data = f.read()
  #Busca los coeficientes
  equation = pyeq2.Models_2D.BioScience.MembraneTransport('SSQABS')
  pyeq2.dataConvertorService().ConvertAndSortColumnarASCII(data,equation,False)
  equation.estimatedCoefficients = [9e4,-13,-8e2,2.5e5]
  equation.Solve()
  equation.CalculateCoefficientAndFitStatistics()
  print 'R-squared:',  equation.r2
  f.close()
  return equation
  
def data_fitterOff(file):
  f=open(file,'r')
  data = f.read()
  #Busca los coeficientes
  equation = pyeq2.Models_2D.BioScience.MembraneTransport('SSQABS','Offset')
  pyeq2.dataConvertorService().ConvertAndSortColumnarASCII(data,equation,False)
  equation.estimatedCoefficients = [9e5,-13,-9e2,2.5e5,200]
  equation.Solve()
  equation.CalculateCoefficientAndFitStatistics()
  print 'R-squared:',  equation.r2
  f.close()
  return equation

def plotter(equation):
  T = equation.dataCache.allDataCacheDictionary['IndependentData'][0]
  E = equation.dataCache.allDataCacheDictionary['DependentData']
  a,b,c,d = equation.solvedCoefficients
  x=arange(-350,T.max()*1.1,T.max()/150)
  F=a*(x-b)/(x*x+c*x+d)
  plot(T,E,'o')
  plot(x,F)

if __name__ == "__main__":
  files=glob.glob('data/P*K.dat')
  for file in files:
    eq = data_fitter(file)
    print file, eq.solvedCoefficients
    plotter(eq)
  xlabel('Temperature [$\\Delta J /k_B$]')
  ylabel('Electric susceptibility $\\chi$')
  title('Simulation Data & fits')
  legend(('DATA','FIT'))
  show()
