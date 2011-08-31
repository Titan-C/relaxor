/*
Sistema de 2 dimensiones de un ferroelectrico relaxor

TO DO
Unidades del sistema
umbtener valor de mu según datos de l a PNR
J intercambio debe contener valores del mu
*/
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <gsl/gsl_rng.h>
#include "sistema.h"
#include "impresor.h"

using namespace std;

int main(int argc, char **argv) {
  //iniciar sistema
  time_t start, end;
  time(&start);
  
//  vaciar archivo de datos en cada ejecución
  system("rm *.dat");
  
  gsl_rng * rng = gsl_rng_alloc (gsl_rng_taus);

  unsigned int L=16, numexps = 4, Equi_iter=200, Exp_iter= 3000;
  double T=12,dT = 0.4, DeltaJ = 1;
  
  gsl_rng_set(rng, time(NULL) );

  clock_t cl_start = clock();
  Sistema relaxor(L, rng, DeltaJ);
  vector<double> temperaturas;
  temperaturas = temp_array(DeltaJ, T, dT, false);
  clock_t cl_stop = clock();
  cout<<"Iniciar sistema "<<cl_stop-cl_start<<"\n";
  // weak field
  cl_start = clock();
  double field0[]={0.1};
  vector<double> campos (field0, field0 + sizeof(field0) / sizeof(double));
  ostringstream frec, fieldamp;
  int taus0[] = {100,50,20,10};
  vector<int> tau (taus0, taus0 + sizeof(taus0) / sizeof(int) );
  
  for(unsigned int t=0;t<tau.size();t++){
    for(unsigned int E=0; E<campos.size(); E++){
      frec.str(""); fieldamp.str("");
      frec<<tau[t]; fieldamp<<campos[E];
      for(unsigned int n = 0; n < numexps; n++){
	for(unsigned int T = 0; T<temperaturas.size();T++){
	  relaxor.experimento(temperaturas[T], campos[E], tau[t], Equi_iter, false, rng, "cool_E"+fieldamp.str()+"_t"+frec.str());
	  relaxor.experimento(temperaturas[T], campos[E], tau[t], Exp_iter, true, rng, "cool_E"+fieldamp.str()+"_t"+frec.str());
	}
      }
      calc_sus(numexps,tau[t],Exp_iter,L,1,temperaturas,campos[E], "cool_E"+fieldamp.str()+"_t"+frec.str() );
    }
  }
  cl_stop = clock();
  cout<<"Experimeto "<<cl_stop-cl_start<<"\n";
  // multi field
  cl_start = clock();
  double field1[]={0.5,1.0,1.5,2.0};
  campos.assign (field1, field1 + sizeof(field1) / sizeof(double));
  int taus1[] = {50,10};
  tau.assign (taus1, taus1 + sizeof(taus1) / sizeof(int) );
  
  for(unsigned int t=0;t<tau.size();t++){
    for(unsigned int E=0; E<campos.size(); E++){
      frec.str(""); fieldamp.str("");
      frec<<tau[t]; fieldamp<<campos[E];
      for(unsigned int n = 0; n < numexps; n++){
	for(unsigned int T = 0; T<temperaturas.size();T++){
	  relaxor.experimento(temperaturas[T], campos[E], tau[t], Equi_iter, false, rng, "cool_E"+fieldamp.str()+"_t"+frec.str());
	  relaxor.experimento(temperaturas[T], campos[E], tau[t], Exp_iter, true, rng, "cool_E"+fieldamp.str()+"_t"+frec.str());
	}
      }
      calc_sus(numexps,tau[t],Exp_iter,L,1,temperaturas,campos[E], "cool_E"+fieldamp.str()+"_t"+frec.str() );
    }
  }
  cl_stop = clock();
  cout<<"Experimeto "<<cl_stop-cl_start<<"\n";
  
  
//   system("gnuplot ../plots.p");
  
  time(&end);
  cout<<difftime(end,start)<<endl;
  return 0;
}

