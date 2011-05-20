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
#include <fstream>
#include <cmath>
#include <vector>
#include <gsl/gsl_rng.h>
#include "sistema.h"
#include "impresor.h"

#define _pi atan(1)*4
#define r_max 2

using namespace std;

int main(int argc, char **argv) {

  gsl_rng * rng = gsl_rng_alloc (gsl_rng_taus);
  
  unsigned int L=8, Niter=400, numexps = 1;
  double T=2.5,dT = 0.1,std;/*
  cin>>L;
  cin>> Niter;
  cin>>T;*/
  gsl_rng_set(rng, time(NULL) );
  
  Sistema relaxor(L, Niter, rng);
  std = relaxor.DeltaJ;
  cout<<"desviación standard = "<<std<<endl;

  
  //iniciar sistema
  time_t start, end;
  time(&start);
  //vaciar archivo de datos en cada ejecución
  file_wipe("log");
  file_wipe("energy_log.dat");
  cout<<"Energía inicial = "<<relaxor.total_E(0)<<endl;
  bool grabar;
  
  file_wipe("sum_sigma_time.dat");
  file_wipe("sum_sigma_conf.dat");
  do{
    grabar = (numexps % 4 == 0);
    //ejecutar una corrida en el tiempo correspondiente a Niter  
    numexps += relaxor.experimento(T, 0, Niter, grabar , rng);
    if (grabar)
      T-=dT;
  }while(T>0);
  
  // Analisis de datos

  //procesar datos polarizacion congelada
  int CongMed=0;
  vector< vector<double> > Polarizacion;
  vector<double> temp;
  unsigned int mediciones = numexps/4;
  
    //Arreglo de temperaturas
  temp.resize(mediciones);
  for(unsigned int i = 0; i<temp.size(); i++)
    temp[i] = T + (mediciones-i)*dT;
  
  CongMed = relaxor.eval_congelamiento_susceptibilidad("Congelamiento.dat",temp,Niter,mediciones);

  
  time(&end);
  cout<<difftime(end,start);
  
  


  return 0;
}
