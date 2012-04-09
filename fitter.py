'''Primera función que ajusta un archivo al mi función tomada'''
import pyeq2
from plotter import *

def dataFitter(file, estimated, upBound, lowBound, weight):
  '''Ajusta los datos de archivo file dadas las listas,
     de coeficientes estimados, cotas superior e inferior,
     más una bandera de datos pesados'''
  f=open(file,'r')
  data = f.read()
  #Crea el objeto de ecuación de ajuste y lo llena de datos
  equation = pyeq2.Models_2D.BioScience.MembraneTransport('SSQABS')
  pyeq2.dataConvertorService().ConvertAndSortColumnarASCII(data,equation,weight)
  f.close()

  #Calcula el ajuste
  equation.estimatedCoefficients = estimated
  equation.upperCoefficientBounds = upBound
  equation.lowerCoefficientBounds = lowBound
  equation.Solve()
  equation.CalculateCoefficientAndFitStatistics()
  print 'R-squared:',  equation.r2
  return equation

def scaleFitter(fit_eq, file, estimated,weight):
  '''Ajusta datos del archivo para mi función arbitraria
     que es un escalamiento de la función Membrane Transport'''
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

def filesFit(path, writefile, estimated=[5e4,-3e2,-9e2,2.2e5], upBound =[72000,0,0,None],lowBound =[None,None,-1040,None], weight=False, plot=True):
  '''Permite llamar a la función de ajuste para todos los archivos que se encuentran en path'''
  files=sort(glob(path))
  for file in files:
    eq = dataFitter(file, estimated, upBound, lowBound, weight)
    fitRecorder(eq,file,writefile)
    print file, eq.solvedCoefficients
    if plot: fittedPlot(eq,file)
  if plot:
      legend()
      axisLabel()

def filesScaleFit(path, material,estimated=[100,0.008],weight=True, plot=True):
  '''Permite llamar a la función de ajuste de escala para todos los archivos que se encuentran en path, para el material deseado'''
  fit_eq = getfitdata(material)
  files=sort(glob(path))
  for file in files:
    eq = scaleFitter(fit_eq[0], file, estimated, weight)
    fitRecorder(eq,file,material+'Fits.csv')
    print file, eq.r2, eq.solvedCoefficients
    if plot: scalefittedPlot(eq,fit_eq[1],file)
  if plot:
      legend()
      axisLabel()

def getfitdata(material):
  '''Requiere el nombre del material y devuelve
      la ecuación de ajuste
      Ej: eq=getfitdata('P2BIT1K') '''
  materialFitdata = open('Expdata')
  fit_eq = ''
  for data in materialFitdata:
    if data.find(material) > 0:
      fit_eq = scale_eqGenerator(data.split()[2:])
  materialFitdata.close()
  return fit_eq,data.split()[2:]

def fitRecorder(equation, datafile, writefile):
  '''Graba los resultados del ajuste, requiero
	el objeto de la ecuación ajustada, archivo de datos de origen
	y archivo donde guardar'''
  f=open(writefile,'a')
  f.write(datafile+'\t')
  f.write(str(equation.r2)+'\t')
  for coef in equation.solvedCoefficients:
    f.write(str(coef)+'\t')
  f.write('\n')
  f.close()

def scale_eqGenerator(coefs):
  '''Devuelve un string con la ecuación a ajustar por escala'''
  if float(coefs[1]) < 0:
    fit_eq = 'u/j*'+coefs[0]+'*(j*X+'+ coefs[1][1:]+')'
  else:
    fit_eq = coefs[0]+'*u*X'
  fit_eq = fit_eq+'/(j**2*X**2'+coefs[2]+'*j*X+'+coefs[3]+')'
  return fit_eq
  
if __name__ == "__main__":
  filesFit(argv[1], argv[2], estimated=[9e4,-200,-8e2,2.5e5], weight=False, plot=True)
