/*
 * Sistema de 2 dimensiones de un ferroelectrico relaxor
 * 
 * TO DO
 * Unidades del sistema
 * umbtener valor de mu según datos de l a PNR
 * J intercambio debe contener valores del mu
 */
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include "sistema.h"
#include "impresor.h"

using namespace std;

int main(int argc, char **argv) {
  //iniciar sistema
  time_t start, end;
  time(&start);
  
  unsigned int L=atoi(argv[1]), numexps = atoi(argv[2]), Equi_iter=atoi(argv[3]), Exp_iter= atoi(argv[4]);
  cout<<argv[0];
  vector<double> tau;
  // Cooling process
  vector<double> rho;
  rho = step2vec(0,0.2,0.01,rho);
  vector<double> temperaturas, campos;
  tau = str2vec("100");
  
  temperaturas.clear();
  temperaturas=step2vec(12,0.1,0.15,temperaturas);
  campos = str2vec("0.35");
  //  for(unsigned int p=0; p<rho.size(); p++)
  Gen_exp(temperaturas,campos,tau,numexps,0.6,L,Equi_iter,Exp_iter,"cool");
  
  time(&end);
  cout<<difftime(end,start)<<endl;
  return 0;
}

