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
//   system("rm *.dat");
  
  gsl_rng * rng = gsl_rng_alloc (gsl_rng_taus);
  
  unsigned int L=16, numexps = 8, Equi_iter=350, Exp_iter= 3000;
  double DeltaJ = 1, rho=0.3;
  
  gsl_rng_set(rng, time(NULL) );
  
  clock_t cl_start = clock();
  Sistema relaxor(L, rng);
  vector<double> temperaturas, campos, tau;
  clock_t cl_stop = clock();
  cout<<"Iniciar sistema "<<cl_stop-cl_start<<"\n";
  
  // weak field
  temperaturas.clear();
  temperaturas=step2vec(DeltaJ,9,0.1,0.5,temperaturas);
  campos = str2vec(DeltaJ, "0.5 1.7 1 2.4 3.5 4");
  tau = str2vec(1,"100 50 30 10 5");
  relaxor.Gen_exp(temperaturas,campos,tau,numexps,DeltaJ,rho,Equi_iter,Exp_iter,"cool",rng);
  
  //Multifield temperature steps
  temperaturas = str2vec(DeltaJ,"0.5 1 1.5 2 3 4 5");
  campos.clear();
  campos = step2vec(DeltaJ, 0.2,9,0.5,campos);
  tau = str2vec(1,"10 50");
  relaxor.Gen_exp(temperaturas,campos,tau,numexps,DeltaJ,rho,Equi_iter,Exp_iter,"riseE",rng);
  
  //Single temp histeresis
  temperaturas = str2vec(DeltaJ, "1 2.5");
  campos = loop2vec(DeltaJ,9,10);
  tau = str2vec(1,"1");
  relaxor.Gen_exp(temperaturas,campos,tau,numexps,DeltaJ,rho,Equi_iter,Exp_iter,"hist_loop",rng);

  
  time(&end);
  cout<<difftime(end,start)<<endl;
  return 0;
}

