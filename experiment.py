import libmat
import numpy as np
import re
from pylab import show,plot,errorbar,figure

#a=libmat.Material(16,0,"Gui95",False,True,True)
#a.inicio(True)
#print a
#thermostat = np.arange(2.5,0,-0.1)
#field = 0.1*np.cos(np.arange(500)*2*np.pi/250)
#wave = libmat.double_vector()
#wave[:] = field

#for n in xrange(8):
  #for t in thermostat:
    #a.state(t,wave,0,True)

#with open('log_pol_fine.dat', 'rb') as f: pol=np.fromfile(f,dtype=np.double)

#pol = pol.reshape(-1,thermostat.size,500)
#avg = pol.mean(2)
#sd = pol.std(2)
#comp=np.arange(0,avg.size,1/500.)
#base=np.arange(avg.size)
##plot(comp,pol.reshape(-1))
#errorbar(base,avg.reshape(-1),yerr=sd.reshape(-1))
#show()
#plot(thermostat, avg.T)
#plot(thermostat, np.abs(avg.T),'--')
#show()

#def explabel(exp_setup):
    #pass

#def stepCalibrator(Niter,tau,min_periods):
    #periods = Niter/tau
    #if periods > min_periods:
        #return periods*tau
    #else:
        #return min_periods*tau
def setLabel(name, L, rho, field, tau, Ti, Tf, dT, Eiter, Qiter):
  label = name + '_N' + str(L*L*L)
  label+= '_p' + str(rho) + '_E' + str(field) + '_t' + str(tau)
  label+= '_Ti'+ str(Ti) + '_Tf'+ str(Tf) + '_dT'+ str(dT)
  label+= '_X' + str(Eiter) + '_Q' + str(Qiter)
  return label

def getLabel(filename):
  pat = r'_N(?P<N>\d+)'
  pat+= r'_p(?P<rho>[\d\.]+)_E(?P<E>[\d\.]+)_t(?P<tau>\d+)'
  pat+= r'_Ti(?P<Ti>[\d\.]+)_Tf(?P<Tf>[\d\.]+)_dT(?P<dT>[-+\d\.]+)'
  pat+= r'_X(?P<Eiter>\d+)_Q(?P<Qiter>\d+)'
  #import pdb; pdb.set_trace()
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
    slow_dipoles.append(np.sum(sig>=bound,axis=2))
  
  fro=figure('fro')
  
  for i in xrange(len(frozen)):
    errorbar(thermostat, slow_dipoles[i].mean(axis=0), slow_dipoles[i].std(axis=0))
  show()


#Gui95 setup experiment 1
label = setLabel("Gui95err", 16, 0, 0, 1, 2.5, 0, -0.1, 400, 1000)
print label
relaxor=libmat.Material(16,0,label,False,True,True)

print relaxor
thermostat = np.arange(2.5,0,-0.1)
equi = libmat.double_vector()
equi[:] = np.zeros(1000)
obs = libmat.double_vector()
obs[:] = np.zeros(400)
for n in xrange(4):
  relaxor.inicio(False)
  for t in thermostat:
    relaxor.state(t,equi,0,False)
    relaxor.state(t,obs,0,True)

with open('logH'+label+'.dat', 'rb') as f: H=np.fromfile(f,dtype=np.double)
plot(H)
show()
eval_frozen('log_sigma_'+label+'.dat')

  
