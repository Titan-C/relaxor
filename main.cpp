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
  
  
  //Parametros de entrada
  unsigned int L=	atoi(argv[1]);
  unsigned int numexps=	atoi(argv[2]);
  vector<double> rho = 	str2vec(argv[3]);
  vector<double> Temp=	str2vec(argv[4]);
  vector<double> campos=str2vec(argv[5]);
  vector<double> tau =	str2vec(argv[6]);
  
  for(unsigned int p=0; p<rho.size(); p++)
    Gen_exp(L,numexps,rho[p],Temp,campos,tau,"cool");
  
  time(&end);
  cout<<difftime(end,start)<<endl;
  return 0;
}

