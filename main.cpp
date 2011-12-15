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
  
  unsigned int L=8, numexps = 4, Equi_iter=350, Exp_iter= 1000;
  double rho = 0.3;
  gsl_rng_set(rng, time(NULL) );
  
  vector<double> tau;
  // Cooling process
  vector<double> DeltaJ;
  DeltaJ = str2vec(1," 0.2 0.5 1 1.5 3");
  
  tau = str2vec(1,"50");
  
  for(unsigned int DJ=0; DJ<DeltaJ.size(); DJ++){
  vector<double> temperaturas, campos;
  temperaturas=step2vec(DeltaJ[DJ],60,0.1,2,temperaturas);
  campos = str2vec(DeltaJ[DJ], "0.3");
    Gen_exp(temperaturas,campos,tau,numexps,DeltaJ[DJ],rho,L,Equi_iter,Exp_iter,"cool",rng);
  }
//   //Multifield temperature steps
//   temperaturas = str2vec(DeltaJ,"0.5 1 1.5 2 3 4 5");
//   campos.clear();
//   campos = step2vec(DeltaJ, 0.2,9,0.5,campos);
//   tau = str2vec(1,"10 50");
//   relaxor.Gen_exp(temperaturas,campos,tau,numexps,DeltaJ,rho,Equi_iter,Exp_iter,"riseE",rng);
//   
//   //Single temp histeresis
//   temperaturas = str2vec(DeltaJ, "1 2.5");
//   campos = loop2vec(DeltaJ,9,10);
//   tau = str2vec(1,"1");
//   relaxor.Gen_exp(temperaturas,campos,tau,numexps,DeltaJ,rho,Equi_iter,Exp_iter,"hist_loop",rng);

  
  time(&end);
  cout<<difftime(end,start)<<endl;
  return 0;
}

