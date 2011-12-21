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
  
  unsigned int L=16, numexps = 10, Equi_iter=350, Exp_iter= 3000;
  gsl_rng_set(rng, time(NULL) );
  
  vector<double> tau;
  // Cooling process
  vector<double> DeltaJ, rho;
  DeltaJ = step2vec(1,0.2,1.5,0.1,DeltaJ);
  rho = step2vec(1,0,0.2,0.01,rho);
  vector<double> temperaturas, campos;
  tau = str2vec(1,"2000 1000 500 300 100 50");
  
  for(unsigned int DJ=0; DJ<DeltaJ.size(); DJ++){
    temperaturas.clear();
    temperaturas=step2vec(DeltaJ[DJ],20,0,0.15,temperaturas);
    campos = str2vec(DeltaJ[DJ], "0.35");
    for(unsigned int p=0; p<rho.size(); p++)
      Gen_exp(temperaturas,campos,tau,numexps,DeltaJ[DJ],rho[p],L,Equi_iter,Exp_iter,"cool",rng);
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

