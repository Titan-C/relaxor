'''Primera función que ajusta un archivo al mi función tomada'''
import pyeq2
from plotter import *
from pylab import *
from glob import glob
from sys import argv

def dataFitter(path, estimated):
  f=open(path,'r')
  data = f.read()
  #Crea el objeto de ecuación de ajuste y lo llena de datos
  equation = pyeq2.Models_2D.BioScience.MembraneTransport('SSQABS')
  pyeq2.dataConvertorService().ConvertAndSortColumnarASCII(data,equation,False)
  f.close()

  #Calcula el ajuste
  equation.estimatedCoefficients = estimated
  equation.Solve()
  equation.CalculateCoefficientAndFitStatistics()
  print 'R-squared:',  equation.r2
  return equation
  
def dataFitterOff(file):
  f=open(file,'r')
  data = f.read()
  #Busca los coeficientes
  equation = pyeq2.Models_2D.BioScience.MembraneTransport('SSQABS','Offset')
  pyeq2.dataConvertorService().ConvertAndSortColumnarASCII(data,equation,False)
  f.close()

  equation.estimatedCoefficients = [9e5,-13,-9e2,2.5e5,200]
  equation.Solve()
  equation.CalculateCoefficientAndFitStatistics()
  print 'R-squared:',  equation.r2
  return equation

def scaleFitter(fit_eq, path, estimated=[100,0.008]):
  f=open(path)
  data = f.read()
  #Crea el objeto de función a escala y lo llena de datos
  equation = pyeq2.Models_2D.UserDefinedFunction.UserDefinedFunction('SSQABS','Default', fit_eq)
  pyeq2.dataConvertorService().ConvertAndSortColumnarASCII(data, equation, False)
  f.close()

  #Calcula el ajuste
  equation.estimatedCoefficients = estimated
  pyeq2.solverService().SolveUsingSimplex(equation)
  equation.CalculateCoefficientAndFitStatistics()
  return equation

def fitRecorder(equation, datafile, targetfile):
  f=open(targetfile,'a')
  f.write(datafile+'\t')
  f.write(str(equation.r2)+'\t')
  for coef in equation.solvedCoefficients:
    f.write(str(coef)+'\t')
  f.write('\n')
  f.close()

def filesFit(path, writefile, estimated=[9e4,-13,-8e2,2.5e5], plot=False):
  files=glob(path)
  
  for file in files:
    eq = dataFitter(file,estimated)
    fitRecorder(eq,file,writefile)
    print file, eq.solvedCoefficients
    if plot: fittedPlot(eq)
  
  if plot: susLabel()

def getfitdata(material):
  materialFitdata = open('Expdata')
  fit_eq = ''
  for data in materialFitdata:
    if data.find(material) > 0:
      fit_eq = scale_eqGenerator(data.split())
  materialFitdata.close()
  return fit_eq

def filesScaleFit(path, material):
  fit_eq = fetfitdata(material)
  files=glob(path)
  for file in files:
    eq = scaleFitter(fit_eq, file)
    fitRecorder(eq,file,material+'Fits')

def scale_eqGenerator(info):
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
    
