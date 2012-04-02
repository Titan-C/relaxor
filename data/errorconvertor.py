from numpy import array, genfromtxt, savetxt
from glob import glob

def include_err(file,errpow):
    data=genfromtxt(file)
    T=data[:,0]
    E=data[:,1]
    d=data[:,2]
    Ex=E*d/100
    err=(1/Ex)**errpow
    new=array([T,E,err,Ex,d])
    savetxt('er'+str(errpow)+file,new.T,fmt='%3.6g',delimiter='\t')
    
if __name__ == "__main__":
  files = glob('P*dat')
  for file in files:
    include_err(file,1/4.0)