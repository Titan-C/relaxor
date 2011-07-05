/*
Sistema de 2 dimensiones de un ferroelectrico relaxor

TO DO
Unidades del sistema
umbtener valor de mu según datos de l a PNR
J intercambio debe contener valores del mu
*/
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <gsl/gsl_rng.h>
#include "sistema.h"
#include "impresor.h"

using namespace std;

int main(int argc, char **argv) {
  //iniciar sistema
  time_t start, end;
  clock_t cl_start = clock();
  time(&start);
  
  //vaciar archivo de datos en cada ejecución
  system("rm *.dat");
  
  gsl_rng * rng = gsl_rng_alloc (gsl_rng_taus);

  unsigned int L=12, numexps = 8, Equi_iter=1000, Exp_iter= 400;
  double T=5,dT = 0.1, H=2, dH=0.1;
  gsl_rng_set(rng, time(NULL) );

  Sistema relaxor(L, Exp_iter, rng);
//Enfriar
  cout<<"Energía inicial = "<<relaxor.total_E(0)<<endl;
  vector<double> campos, temperaturas;
  temperaturas = temp_array(relaxor.DeltaJ, T, dT, false);
  campos = field_array(relaxor.DeltaJ, H, dH);
  for(unsigned int E=0;E < campos.size(); E++){
    for(unsigned int n = 0; n < numexps; n++){
      for(unsigned int T = 0; T < temperaturas.size(); T++){
	relaxor.experimento(temperaturas[T], 2.5*campos[E], Equi_iter, false , rng,"cool");
	relaxor.experimento(temperaturas[T], 2.5*campos[E], Exp_iter, true , rng,"cool");
      }
    }
  }
  // Analisis de datos
  Sus_proc(numexps, Exp_iter, L, relaxor.DeltaJ, temperaturas, campos,"cool");
  Pol_proc(numexps, Exp_iter, L, relaxor.DeltaJ, temperaturas, campos,"cool");
  plot_data_sus(relaxor.DeltaJ, temperaturas, campos, "cool");

//Calentar desde pol
  Equi_iter = 50;
  temperaturas = temp_array(relaxor.DeltaJ, 3, dT, true);
  campos.assign(3,0);
  campos[1]=0.2*relaxor.DeltaJ;
  campos[2]=0.6*relaxor.DeltaJ;
  for(unsigned int E=0;E < campos.size(); E++){
    for(unsigned int n = 0; n < numexps; n++){
      relaxor.init_pol(rng,true);
      for(unsigned int T = 0; T < temperaturas.size(); T++){
	relaxor.experimento(temperaturas[T], 2.5*campos[E], Equi_iter, false , rng,"heatpol");
	relaxor.experimento(temperaturas[T], 2.5*campos[E], Exp_iter, true , rng,"heatpol");
      }
    }
  }
  // Analisis de datos
  Pol_proc(numexps, Exp_iter, L, relaxor.DeltaJ, temperaturas, campos,"heatpol");
//Calentar desde unpol
  for(unsigned int E=0;E < campos.size(); E++){
    for(unsigned int n = 0; n < numexps; n++){
      relaxor.init_pol(rng,false);
      for(unsigned int T = 0; T < temperaturas.size(); T++){
	relaxor.experimento(temperaturas[T], 2.5*campos[E], Equi_iter, false , rng,"heatunpol");
	relaxor.experimento(temperaturas[T], 2.5*campos[E], Exp_iter, true , rng,"heatunpol");
      }
    }
  }
  // Analisis de datos
  Pol_proc(numexps, Exp_iter, L, relaxor.DeltaJ, temperaturas, campos,"heatunpol");
  
//Enfriar T<Tf
  Equi_iter = 400; string file="coolpol"; campos.assign(1,0);
  for(unsigned int i=0;i<4;i++){
    temperaturas = temp_array(relaxor.DeltaJ, 1-i*0.2, 0.05, false);
    for(unsigned int n = 0; n < numexps; n++){
      relaxor.init_pol(rng,true);
      for(unsigned int T = 0; T < temperaturas.size(); T++){
	relaxor.experimento(temperaturas[T], 0, Equi_iter, false , rng,file.c_str());
	relaxor.experimento(temperaturas[T], 0, Exp_iter, true , rng,file.c_str());
      }
    }
    Pol_proc(numexps, Exp_iter, L, relaxor.DeltaJ, temperaturas, campos ,file.c_str());
    file+="_t";
  }
  time(&end);
  cout<<difftime(end,start)<<endl;
  cout<<(double)(clock()-cl_start)/CLOCKS_PER_SEC;
  return 0;
}

