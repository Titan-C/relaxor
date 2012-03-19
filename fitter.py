'''Primera función que ajusta un archivo al mi función tomada'''
import pyeq2
from plotter import *
from pylab import *
from glob import glob

def dataFitter(file):
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
  
def dataFitterOff(file):
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

def scaleFitter(fit_eq,path):
  record = 0.0
  files=glob(path)
  for file in files:
    S=open(file)
    data = S.read()
    equation = pyeq2.Models_2D.UserDefinedFunction.UserDefinedFunction('SSQABS','Default', fit_eq)
    pyeq2.dataConvertorService().ConvertAndSortColumnarASCII(data, equation, False)
    equation.estimatedCoefficients = [100,0.008]
    pyeq2.solverService().SolveUsingSimplex(equation)
    
    #imprimir resultados
    equation.CalculateCoefficientAndFitStatistics()
    if equation.r2 > record:
      record = equation.r2
      print "Material", file
      print 'R-squared:',  equation.r2, equation.solvedCoefficients
    
  
  
def recorder(equation,file):
  f=open('Expdata','a')
  f.write(file+'\t')
  f.write(str(equation.r2)+'\t')
  for coef in equation.solvedCoefficients:
    f.write(str(coef)+'\t')
  f.write('\n')
  f.close()

def experimentFit(plot=False):
  files=glob('data/P*K.dat')
  
  for file in files:
    eq = dataFitter(file)
    recorder(eq,file)
    print file, eq.solvedCoefficients
    if plot: fittedPlot(eq)
  
  if plot: susLabel()


def equationGenerator(info):
  if float(info[3]) < 0:
    fit_eq = '+'+ info[3][1:]
  else:
    fit_eq = '-'+ info[3]
  fit_eq = 'u/j*'+info[2]+'*(j*X'+fit_eq+')/(j**2*X*X+'+info[4]+'*j*X+'+info[5]+')'
  return fit_eq
  
if __name__ == "__main__":
  Exper=open('Expdata')
  for exper in Exper:
    print '\n Para el material ', exper.split()[0]
    fit_eq = equationGenerator( exper.split() )
    scaleFitter(fit_eq, argv[1])
  Exper.close()
    
