import libmat
import numpy as np
import re
from pylab import show,plot,errorbar,figure

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

def eval_frozen(filename, frozen=[1,0.9]):
  """Evaluate frozen dipoles proportion for determinate
     experiment(filename) given the list of thresholds
     to consider a dipole frozen
     eval_frozen(filename,[1,0.95])"""
  info = getLabel(filename)
  with open(filename, 'rb') as f:
    sigma=np.fromfile(f,dtype=np.int32)/float(info['Eiter'])

  thermostat = np.arange(float(info['Ti']), float(info['Tf']), float(info['dT']))
  sig=np.abs(sigma).reshape(-1,thermostat.size,int(info['N']))
  slow_dipoles=[]
  for bound in frozen:
    slow_dipoles.append(np.sum(sig>=bound,axis=2)/float(info['N']))

  fro=figure('fro')
  for i in xrange(len(slow_dipoles)):
    errorbar(thermostat, slow_dipoles[i].mean(axis=0), slow_dipoles[i].std(axis=0))
  show()

  return slow_dipoles, thermostat

def susGUI(filename, frozen=[0.9, 0.8, 0.6]):
  """Calculate dielectric susceptibility X=(1-p)/T
     where p:proportion of slow dipoles & T:temperature"""
  slow_dipoles, thermostat = eval_frozen(filename, frozen)

  sus=figure('sus')
  plot(thermostat, 1/thermostat, '--')
  for i in range(len(slow_dipoles)):
    plot(thermostat,(1-slow_dipoles[i].mean(axis=0))/thermostat,'*')
  show()
