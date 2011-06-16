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
  
  gsl_rng * rng = gsl_rng_alloc (gsl_rng_taus);

  unsigned int L=8, numexps = 10, Equi_iter=1000, Exp_iter= 400;
  double T=2.5,dT = 0.1;
  gsl_rng_set(rng, time(NULL) );

  Sistema relaxor(L, Exp_iter, rng);

  //vaciar archivo de datos en cada ejecución
  file_wipe("log");
  file_wipe("sum_sigma_time.dat");
  file_wipe("sum_sigma_conf.dat");
  file_wipe("energy_log.dat");
  cout<<"Energía inicial = "<<relaxor.total_E(0)<<endl;
  vector<double> campos, temperaturas;
  temp_array(temperaturas, relaxor.DeltaJ, T, dT);
  field_array(campos, relaxor.DeltaJ);

  for(unsigned int i=0;i < campos.size(); i++){
    for(unsigned int j = 0; j < numexps; j++){
      for(unsigned int k = 0; k < temperaturas.size(); k++){
	relaxor.experimento(temperaturas[k], campos[i], Equi_iter, false , rng);
	relaxor.experimento(temperaturas[k], campos[i], Exp_iter, true , rng);
      }
    }
  }

  // Analisis de datos
  procesar(numexps, Exp_iter, L, relaxor.DeltaJ, temperaturas, campos);

  time(&end);
  cout<<difftime(end,start)<<endl;
  cout<<(double)(clock()-cl_start)/CLOCKS_PER_SEC;
  system("gnuplot ../plots.p");
  
  return 0;
}
