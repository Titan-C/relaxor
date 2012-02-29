
'''Primera función que ajusta un archivo al mi función tomada'''
import os, sys, inspect

# ensure pyeq2 can be imported
if os.path.join(sys.path[0][:sys.path[0].rfind(os.sep)], '../..') not in sys.path:
    sys.path.append(os.path.join(sys.path[0][:sys.path[0].rfind(os.sep)], '../..'))
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
#  equation.estimatedCoefficients = [9e4,-1.3e2,-8e2,2.5e5]
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
  equation.estimatedCoefficients = [9e5,-1.3e2,-9e2,2.5e5,200]
  equation.Solve()
  equation.CalculateCoefficientAndFitStatistics()
  print 'R-squared:',  equation.r2
  f.close()
  return equation

def plotter(equation):
  T = equation.dataCache.allDataCacheDictionary['IndependentData'][0]
  E = equation.dataCache.allDataCacheDictionary['DependentData']
  a,b,c,d = equation.solvedCoefficients
  x=arange(0,T.max()*1.1,T.max()/150)
  F=a*(x-b)/(x*x+c*x+d)
  plot(T,E,'o')
  plot(x,F)

if __name__ == "__main__":
  files=glob.glob('/home/oscar/workplace/Relaxor/ICTPdat/sus*.dat')
  for file in files:
    eq = data_fitter(file)
    print eq.solvedCoefficients
    #plotter(eq)
    #eq0 = data_fitterOff(file)
    #print eq0.solvedCoefficients
    plotter(eq)
  xlabel('Temperature [$\\Delta J /k_B$]')
  ylabel('Electric susceptibility $\\chi$')
  title('Simulation Data & fits')
  legend(('DATA','FIT'))
  show()
  
#Fitted Parameters:
    #a = 1.2370050415987713E+05
    #b = 2.6765810051540325E+02
    #c = -9.9125857400687005E+02
    #d = 2.5241678263901567E+05
#R-squared: 0.987034132705
#Root Mean Squared Error (RMSE): 143.790205244 
