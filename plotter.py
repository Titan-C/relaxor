# -*- coding: utf-8 -*-
from pylab import *
from sys import argv
from glob import glob
from utils import zoom_effect02

#Plotters
# Single files plots

def HlogPlot(file,stl='+'):
  H=genfromtxt(file)
  plot(H,stl,label=legendSet(file))

def uPlot(file, lb, stl='-o', error=False, IM=False,col=1):
  data = genfromtxt(file)
  T = data[:,0]
  Y1 = data[:,col]
#  lb='$\\chi$ $\\mathcal{R}eal$'
  if not error:
    return plot(T,Y1,stl,label=lb)
  else:
    Y1e = data[:,2]
    return errorbar(T,Y1,yerr=Y1e,fmt=stl,label=lb)
  if IM:
    Y2 = data[:,4]
    return plot(T,Y2,stl,label='$\\chi$ $\\mathcal{I}maginaria$')

def fittedPlot(equation,file):
  T = equation.dataCache.allDataCacheDictionary['IndependentData'][0]
  E = equation.dataCache.allDataCacheDictionary['DependentData']
  a,b,c,d = equation.solvedCoefficients
  x=arange(T.min()*0.8,T.max()*1.1,T.max()/150)
  F=a*(x-b)/(x*x+c*x+d)
  plot(T,E,'o',label=legendSet(file)+'Hz')
  plot(x,F,label='Ajuste '+legendSet(file)+'Hz')

def scalefittedPlot(equation,coefs,file):
  T = equation.dataCache.allDataCacheDictionary['IndependentData'][0]
  E = equation.dataCache.allDataCacheDictionary['DependentData']
  a,b,c,d = [float(x) for x in coefs]
  j,u = equation.solvedCoefficients
  x=arange(T.min()*0.8,T.max()*1.1,T.max()/150)
  F=u/j*a*(j*x-b)/(j**2*x*x+c*j*x+d)
  plot(T,E,'o',x,F,label=legendSet(file))

def sigmahist(rho,tempind,sigmas = None,tr=0.7):
  fig=figure()
  file = glob('sigmas*p'+str(rho)+'_*')[0]
  n=int(file[file.find('_n')+2])
  if sigmas == None:
    sigmas = genfromtxt(file)
  T = genfromtxt('pol'+file[6:])[:,0]
  Sigstat = array([concatenate([sigmas[temps+i*len(T)] for i in range(n)]) for temps in range(len(T))])
  for temp in tempind:
#    n,bins,patches = hist(Sigstat[temp],96,normed=True, label='T='+str(T[temp])+'[$\\Delta J /k_B$]',alpha=tr)
    hist(Sigstat[temp]*rand(len(Sigstat[temp])),96,normed=True,label='T='+str(T[temp])+'[$\\Delta J /k_B$]',alpha=tr)
#    w=bins[:-1]*n
#    plot(bins, normpdf(bins, mean(w),std(w)/2), 'g+-', label='fit', linewidth=3)
  legend()
  rho,E,tau,L = simIdentifier(file)
  dist_sig_avgTitle(rho,E,tau)
  return sigmas
    
#Multifile plots
def filesPlot(files,Frho,FE,Ftau,yAxis,stl,error,IM):
  plots=list()
  for file in files:
    rho,E,tau,L=simIdentifier(file)
#    if float(rho) < 0.55:
    plots.append(uPlot(file,legendSet(file,Frho,FE,Ftau),stl,error,IM)[0])
  ylabel(yAxis)
  return plots

def ufilePlot(procs,rho=None,E=None,tau=None,stl='-o',error=False,IM=False):
  '''Grafica para los procesos *procs*(pol,sus,frozen) los datos que
     correspondan a simulaciones con rho, E y tau fijos'''

  fig=figure()
  ax1 = subplot(111)
  procs = procs.split()
  files = sort(glob(procs[0]+searchstr(rho,E,tau)))
  lin=filesPlot(files,rho,E,tau,axisLabel(procs[0]),stl,error,IM)
  #legend()
  if len(procs)>1:
    ax2 = twinx()
    files = [procs[1]+file[3:] for file in files]
    lin+=filesPlot(files,rho,E,tau,axisLabel(procs[1]),'gp-',error,IM)
  print lin
  ax1.set_xlabel(tempLabel())
  legend(lin,[l.get_label() for l in lin])
  title(procTitle(procs,rho,E,tau))

  show()
  
def matplot(material,stl='-o'):
  '''Grafica Las curvas de constante dieléctrica y pérdidas para el
     material asignado colocar 'P2BIT' o 'P3BIT' '''
  fig=figure()
  ax1 = subplot(111)  
  freqs = ['1K','10K','100K']
  for file in freqs:
    uPlot('data/'+material+file+'.dat',file+'Hz',stl)
  xlabel(tempLabelexp())
  ylabel(dieLabel())
  legend(loc=9)
  ax2 = twinx()
  ax2.set_xlim(290, 880)
  ax2.set_ylim(0, 50)
  for file in freqs:
    uPlot('data/'+material+file+'.dat',file,stl,col=2)
  ylabel(lossLabel())
  arrow(380,36,-30,0,width=1,color='k')
  arrow(430,8,30,0,width=1,color='k')
  
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

# Labelers
def tempLabel():
  return 'Temperatura [$\\Delta J /k_B$]'

def tempLabelexp():
  return 'Temperatura [$^{\circ}K$]'

def sig_avgLabel():
  return '$\\overline{p}$'

def dist_sig_avgLabel():
  return '$\\mathcal{P}(\\overline{p})$'

def polLabel():
  return u'Polarización normada absoluta $[\\overline{\\mu}/N]$'

def frozenLabel():
  return u'Fracción de PNRs congeladas'

def susLabel():
  return u'Susceptibilidad dieléctrica $\\chi$'

def dieLabel():
  return u'Constante dieléctrica $\\varepsilon_r$'

def lossLabel():
  return u'Pérdidas dieléctricas %'

def axisLabel(proc):
  if proc == 'sus':
    return susLabel()
  if proc == 'pol':
    return polLabel()
  if proc == 'fro':
    return frozenLabel()

def procTitle(proc,rho,E,tau):
  setup = u''
  for i in range(len(proc)):
    if proc[i] == 'sus':
      setup += u'Susceptibilidad'
    if proc[i] == 'pol':
      setup += u'Polarización espontánea'
    if proc[i] == 'fro':
      setup += u'PNRs congeladas'
    if i<len(proc)-1:
      setup += u' y '
  setup +=' en función de la tempertatura'

  return setup+fixedCond(rho,E,tau)

def fixedCond(rho,E,tau):
  cond = str()
  if rho!=None:
    cond += '$\\rho = $'+str(rho)
  if E!=None:
    cond += '; $E_0 =$'+str(E)
  if tau!=None and E!=0:
    cond += '; $\\tau =$'+str(tau)
  if len(cond)==0:
    return ''
  if cond[0]==';':
    cond = cond[1:]
  return '\npara: '+cond

def fitTitle():
  title(u'Curvas de ajuste de los datos experimentales')

def HlogLabel(rho,E,tau):
  xlabel('Iteraciones $[MCS/dipolo]$')
  ylabel(u'Energía del sistema $[\Delta J]$')
  suptitle(procTitle(u'Evolución de la energía',rho,E,tau))

def dist_sig_avgTitle(rho,E,tau):
  xlabel(sig_avgLabel())
  ylabel(dist_sig_avgLabel())
  title(u'Densidad de probabilidad de la polarización local'+fixedCond(rho,E,tau))

def simIdentifier(file):
  rho = file.find('_p')
  E = file.find('_E')
  tau = file.find('_t')
  L = file.find('_L')

  rho = file[rho+2:E]
  E = file[E+2:tau]
  tau = file[tau+2:L]
  L = file[L+2:L+4]

  return rho,E,tau,L

def legendSet(file,Frho=None,FE=None,Ftau=None):
  if file.find('_p') > 0:
    rho,E,tau,L = simIdentifier(file)

    legend = str()
    if Frho==None:
      if rho == '0':
        rho = '0.0'
      legend += '$\\rho = $'+rho
    if FE==None:
      legend += '; $E_0 =$'+E
    if Ftau==None and E!=0:
      legend += '; $\\tau =$'+tau
    if len(legend)==0:
      legend = '$N=$'+L.split('_')[0]+'$^3$'
      if file[:3] == 'sus': legend='Susceptibilidad'
      elif file[:3] == 'pol': legend = u'Polarización \n espontánea'
      else: legend = 'PNRs congeladas'
    if legend[0]==';':
      legend = legend[1:]
    return legend
    
  else:
    mat = file.find('P')
    return file[mat+5:-4]

def searchstr(rho,E,tau):
  Sstr='*'
  if rho!=None:
    Sstr+='p'+str(rho)+'_*'
  if E!=None:
    Sstr+='E'+str(E)+'_*'
  if tau!=None:
    Sstr+='t'+str(tau)+'_*'
  return Sstr

def p2anot():
  #sus
  #axvline(x=5.3,ls='--')
  #arrow(5,0.295,-1.5,0,color='k', head_length=0.25,head_width=0.01)
  #annotate('$T_m$',xy=(5.3,0),xytext=(5.5,0.05),arrowprops=dict(arrowstyle="->"))
  #pol
  plot([3.9,4.9],[0.35,0],'g--')
  arrow(4.4,0.2,0.75,0,color='k', head_length=0.25,head_width=0.01)
  annotate('$T_C$',xy=(4.9,0),xytext=(4.2,0.05),arrowprops=dict(arrowstyle="->"))
  #frozen
  axvline(x=4,ls='--')
  arrow(1.6,0.16,-0.55,0,color='k', head_length=0.25,head_width=0.01)
  annotate('$T_f$',xy=(4,0),xytext=(3.5,0.05),arrowprops=dict(arrowstyle="->"))

def p3anot():
  #sus
  #axvline(x=4.3,ls='--')
  #arrow(3,0.04,-0.75,0,color='k', head_length=0.5)
  #annotate('$T_m$',xy=(4.3,0),xytext=(4.5,0.01),arrowprops=dict(arrowstyle="->"))
  #pol
  plot([2.9,3.9],[0.09,0],'g--')
  arrow(3.65,0.025,0.75,0,color='k', head_length=0.5)
  annotate('$T_C$',xy=(3.9,0),xytext=(3.1,0.007),arrowprops=dict(arrowstyle="->"))
  #frozen
  axvline(x=2.9,ls='--')
  arrow(1.3,0.036,-0.55,0,color='k', head_length=0.5)
  annotate('$T_f$',xy=(2.9,0),xytext=(2.5,0.007),arrowprops=dict(arrowstyle="->"))
  


if __name__ == "__main__":
  ufilePlot('sus',0.1,0.3,1000)
