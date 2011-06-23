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

  unsigned int L=8, numexps = 4, Equi_iter=1000, Exp_iter= 400;
  double T=2.5,dT = 0.1, H=2, dH=0.1;
  gsl_rng_set(rng, time(NULL) );

  Sistema relaxor(L, Exp_iter, rng);
//Enfriar
  cout<<"Energía inicial = "<<relaxor.total_E(0)<<endl;
  vector<double> campos, temperaturas;
  temp_array(temperaturas, relaxor.DeltaJ, T, dT, false);
  field_array(campos, relaxor.DeltaJ, H, dH);
  for(unsigned int E=0;E < campos.size(); E++){
    for(unsigned int n = 0; n < numexps; n++){
      for(unsigned int T = 0; T < temperaturas.size(); T++){
	relaxor.experimento(temperaturas[T], 2.5*campos[E], Equi_iter, false , rng);
	relaxor.experimento(temperaturas[T], 2.5*campos[E], Exp_iter, true , rng);
      }
    }
  }

  // Analisis de datos
  procesar(numexps, Exp_iter, L, relaxor.DeltaJ, temperaturas, campos);
  //Graficos
  graficos(relaxor.DeltaJ, temperaturas, campos);

  time(&end);
  cout<<difftime(end,start)<<endl;
  cout<<(double)(clock()-cl_start)/CLOCKS_PER_SEC;
  system("gnuplot ../plotsf.p");
  system("rename .dat _1.dat *.dat");
  
  
  
  return 0;
}
