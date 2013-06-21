#Simulations for paper
from experiment import *
#Gui, H., Gu, B., & Zhang, X. (1995). Dynamics of the freezing process in relaxor ferroelectrics.pdf. Physical Review B, 52(5), 3135.
#Gui95 setup experiment 1, figs 2 & 3
label = genLabel("Gui95err", 16, 0, 0., 1, 2.5, 0, -0.1, 400, 1000)
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
susGUI('log_sigma_'+label+'.dat')

#Gui95 setup experiment 2, figs 4 & 5
relaxor2=libmat.Material(16,0,"NULL",False,False,True)
obs = libmat.double_vector()
equi = libmat.double_vector()
thermostat = np.arange(2.5,0,-0.1)
for E in [0., 0.8, 2.0, 10.]:
  label = genLabel("Gui95err", 16, 0, E, 1, 2.5, 0, -0.1, 400, 1000)
  print label
  relaxor2.set_ExpID(label)
  equi[:] = E*np.ones(1000)
  obs[:] = E*np.ones(400)
  for n in xrange(4):
    relaxor2.inicio(False)
    for t in thermostat:
      relaxor2.state(t,equi,0,False)
      relaxor2.state(t,obs,0,True)

  susGUI('log_sigma_'+label+'.dat', [0.9])