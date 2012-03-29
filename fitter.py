'''Primera función que ajusta un archivo al mi función tomada'''
import pyeq2
from plotter import *
from glob import glob
from sys import argv

def dataFitter(file, estimated, upBound, weight):
  f=open(file,'r')
  data = f.read()
  #Crea el objeto de ecuación de ajuste y lo llena de datos
  equation = pyeq2.Models_2D.BioScience.MembraneTransport('SSQABS')
  pyeq2.dataConvertorService().ConvertAndSortColumnarASCII(data,equation,weight)
  f.close()

  #Calcula el ajuste
  equation.estimatedCoefficients = estimated
  equation.upperCoefficientBounds = upBound
  equation.Solve()
  equation.CalculateCoefficientAndFitStatistics()
  print 'R-squared:',  equation.r2
  return equation

def scaleFitter(fit_eq, file, estimated,weight):
  f=open(file)
  data = f.read()
  #Crea el objeto de función a escala y lo llena de datos
  equation = pyeq2.Models_2D.UserDefinedFunction.UserDefinedFunction('SSQABS','Default', fit_eq)
  pyeq2.dataConvertorService().ConvertAndSortColumnarASCII(data, equation,weight)
  f.close()

  #Calcula el ajuste
  equation.estimatedCoefficients = estimated
  pyeq2.solverService().SolveUsingSimplex(equation)
  equation.CalculateCoefficientAndFitStatistics()
  return equation

def filesFit(path, writefile, estimated=[9e4,-200,-8e2,2.5e5], upBound =[None,0,0,None], weight=False, plot=True):
  files=glob(path)
  
  for file in files:
    eq = dataFitter(file, estimated, upBound, weight)
    fitRecorder(eq,file,writefile)
    print file, eq.solvedCoefficients
    if plot: fittedPlot(eq,file)
  
  if plot: axisLabel()

def filesScaleFit(path, material,estimated=[100,0.008],weight=True):
  fit_eq = getfitdata(material)
  files=glob(path)
  for file in files:
    eq = scaleFitter(fit_eq, file, estimated, weight)
    fitRecorder(eq,file,material+'Fits')

def getfitdata(material):
  materialFitdata = open('Expdata')
  fit_eq = ''
  for data in materialFitdata:
    if data.find(material) > 0:
      fit_eq = scale_eqGenerator(data.split())
  materialFitdata.close()
  return fit_eq

def fitRecorder(equation, datafile, writefile):
  f=open(writefile,'a')
  f.write(datafile+'\t')
  f.write(str(equation.r2)+'\t')
  for coef in equation.solvedCoefficients:
    f.write(str(coef)+'\t')
  f.write('\n')
  f.close()

def scale_eqGenerator(info):
  if float(info[3]) < 0:
    fit_eq = '+'+ info[3][1:]
  else:
    fit_eq = '-'+ info[3]
  fit_eq = 'u/j*'+info[2]+'*(j*X'+fit_eq+')/(j**2*X*X+'+info[4]+'*j*X+'+info[5]+')'
  return fit_eq
  
if __name__ == "__main__":
  filesFit(argv[1], argv[2], estimated=[9e4,-200,-8e2,2.5e5], weight=False, plot=True)
