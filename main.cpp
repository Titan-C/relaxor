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
  time(&start);

  gsl_rng * rng = gsl_rng_alloc (gsl_rng_taus);
  
  unsigned int L=8, Niter=400, numexps = 1;
  double T=2.5,dT = 0.1;/*
  cin>>L;
  cin>> Niter;
  cin>>T;*/
  gsl_rng_set(rng, time(NULL) );
  
  Sistema relaxor(L, Niter, rng);
  cout<<"desviación standard = "<<relaxor.DeltaJ<<endl;

  //vaciar archivo de datos en cada ejecución
  file_wipe("log");
  file_wipe("sum_sigma_time.dat");
  file_wipe("sum_sigma_conf.dat");
  file_wipe("energy_log.dat");
  cout<<"Energía inicial = "<<relaxor.total_E(0)<<endl;
  vector<double> campos, temperaturas;
  temp_array(temperaturas, T, dT);
  field_array(campos);
  
  for(unsigned int i=0;i < campos.size(); i++){
    for(unsigned int j = 0; j < numexps; j++){
      for(unsigned int k = 0; k < temperaturas.size(); k++){
	for(unsigned int l = 0; l < 3;l++)
	  relaxor.experimento(temperaturas[k], campos[i], Niter, false , rng);
	relaxor.experimento(temperaturas[k], campos[i], Niter, true , rng);
      }
    }
  }
  
  // Analisis de datos
  procesar(Niter, L, temperaturas);

  time(&end);
  cout<<difftime(end,start)<<endl;
  cout<<(double)clock()/CLOCKS_PER_SEC;
  system("gnuplot ../plots.p");
  
  return 0;
}
