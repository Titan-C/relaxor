import libmat
import numpy as np
import re
from pylab import show,plot,errorbar,figure
from scipy.integrate import simps
from os import path, remove

#def stepCalibrator(Niter,tau,min_periods):
    #periods = Niter/tau
    #if periods > min_periods:
        #return periods*tau
    #else:
        #return min_periods*tau

def genLabel(name, L, rho, field, tau, Ti, Tf, dT, Eiter, Qiter):
  """Generate string to label simulation files"""
  label = name + '_N' + str(L*L*L)
  label+= '_p' + str(rho) + '_E' + str(field) + '_t' + str(tau)
  label+= '_Ti'+ str(Ti) + '_Tf'+ str(Tf) + '_dT'+ str(dT)
  label+= '_X' + str(Eiter) + '_Q' + str(Qiter)
  return label

def getLabel(filename):
  """Extract from filename simulation setup"""
  pat = r'_N(?P<N>\d+)'
  pat+= r'_p(?P<rho>[\d\.]+)_E(?P<E>[\d\.]+)_t(?P<tau>\d+)'
  pat+= r'_Ti(?P<Ti>[\d\.]+)_Tf(?P<Tf>[\d\.]+)_dT(?P<dT>[-+\d\.]+)'
  pat+= r'_X(?P<Eiter>\d+)_Q(?P<Qiter>\d+)'
  return re.search(pat, filename).groupdict()

def buildTemp(experiment):
  """Generate numpy array for temperatures"""
  return np.arange(float(experiment['Ti']), float(experiment['Tf']), float(experiment['dT']))

def eval_frozen(filename, frozen=[1,0.9]):
  """Evaluate frozen dipoles proportion for determinate
     experiment(filename) given the list of thresholds
     to consider a dipole frozen
     eval_frozen(filename,[1,0.95])"""
  experiment = getLabel(filename)
  thermostat = buildTemp(experiment)
  with open(filename, 'rb') as f:
    sigma=np.fromfile(f,dtype=np.int32)/float(experiment['Eiter'])
  sig=np.abs(sigma).reshape(-1, thermostat.size, int(experiment['N']))

  slow_dipoles=[]
  for bound in frozen:
    slow_dipoles.append( np.sum( sig>=bound, axis=2 ) / float( experiment['N'] ) )

  fro=figure('fro')
  for i in xrange(len(slow_dipoles)):
    errorbar( thermostat, slow_dipoles[i].mean(axis=0), slow_dipoles[i].std(axis=0) )

  return slow_dipoles, thermostat

def getPolLog(filename):
  experiment = getLabel(filename)
  thermostat =  buildTemp(experiment)
  with open(filename, 'rb') as f:
    pol_log=np.fromfile(f,dtype=np.double).reshape( -1, thermostat.size, int(info['Eiter']) )
  return pol_log, thermostat, experiment

def eval_pol(filename):
  """Calculate polarization"""
  pol_log, thermostat, experiment = getPolLog(filename)
  pol=pol_log.mean(axis=2)
  pol_err=pol_log.std(axis=2)

  plot(thermostat, pol.mean(axis=0), 'k*-',thermostat, np.abs(pol).mean(axis=0))  

  return pol, pol_err, thermostat

def eval_sus(filename):
  """Calculate dielectric susceptibility according to Liu paper"""
  pol_log, thermostat, experiment = getPolLog(filename)

  periods = int(experiment['Eiter'])/int(experiment['tau'])
  t=np.arange(0, 2*np.pi*periods, 2*np.pi/int(experiment['tau']) )
  cos_wave=np.cos(t)
  sin_wave=np.sin(t)

  susRe=simps(cos_wave*pol_log)/int(experiment['Eiter'])/float(experiment['E'])
  susIm=simps(sin_wave*pol_log)/int(experiment['Eiter'])/float(experiment['E'])

  plot(thermostat, susRe.mean(axis=0), thermostat, susIm.mean(axis=0))

  return susRe, susIm, thermostat
  
def susGUI(filename, frozen=[0.9, 0.8, 0.6]):
  """Calculate dielectric susceptibility X=(1-p)/T
     where p:proportion of slow dipoles & T:temperature"""
  slow_dipoles, thermostat = eval_frozen(filename, frozen)

  for i in range(len(slow_dipoles)):
    plot(thermostat,(1-slow_dipoles[i].mean(axis=0))/thermostat,'*')

def requiredSimulations(numexps, filename):
  """Calculate the number of needed remaining simulations in order
     to reach the required number"""
  try:
    size = path.getsize(filename)
  except:
    return numexps

  info = getLabel(filename)
  thermostat = np.arange(float(info['Ti']), float(info['Tf']), float(info['dT']))
  n = size/8./float(info['Eiter'])/thermostat.size
  if int(n) == n:
    if numexps - n > 0:
      return numexps - n
    else:
      return 0
  else:
    remove(filename)
    return numexps
