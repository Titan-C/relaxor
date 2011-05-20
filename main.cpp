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
  file_wipe("sum_sigma_time.dat");
  file_wipe("sum_sigma_conf.dat");
  file_wipe("energy_log.dat");
  cout<<"Energía inicial = "<<relaxor.total_E(0)<<endl;
  bool grabar;

  do{
    grabar = (numexps % 4 == 0);
    //ejecutar una corrida en el tiempo correspondiente a Niter  
    numexps += relaxor.experimento(T, 0, Niter, grabar , rng);
    if (grabar)
      T-=dT;
  }while(T>0);
  
  // Analisis de datos

  //procesar datos polarizacion congelada
  vector< vector<double> > sigmas_time, sigmas_conf;
  vector< vector<double> > S_frozen, Susceptibilidad, Polarizacion;
  vector<double> temp;
  unsigned int lentos=4, mediciones = numexps/4, L3=L*L*L;
  
    //Arreglo de temperaturas
  temp.resize(mediciones);
  for(unsigned int i = 0; i<temp.size(); i++)
    temp[i] = T + (mediciones-i)*dT;
  
  S_frozen.resize(mediciones);
  for(unsigned int i=0; i< S_frozen.size(); i++){
    S_frozen[i].resize(lentos+1);
    S_frozen[i][lentos] = temp[i];
    for(unsigned int j=0; j< lentos; j++)
      S_frozen[i][j] = 0;
  }

  import_data(sigmas_time, "sum_sigma_time.dat", mediciones, L3);
  
  for(unsigned int i=0 ; i<sigmas_time.size(); i++){
    for(unsigned int j=0; j<sigmas_time[i].size(); j++){
      sigmas_time[i][j] = abs(sigmas_time[i][j]/Niter);
      for(unsigned int k=0;k<lentos;k++){
	if (sigmas_time[i][j] >= (1-0.1*k) )
	  S_frozen[i][k]+=(double) 1/L3;
      }
    }    
  }
  array_print(S_frozen, "Congelamiento.dat");
  

  
  //Encontrar la susceptibilidad del material
  Susceptibilidad.resize(mediciones);
  for(unsigned int i=0; i<Susceptibilidad.size(); i++){
    Susceptibilidad[i].resize(lentos+1);
    Susceptibilidad[i][lentos] = temp[i];
    for(unsigned int j=0; j<lentos; j++)
      Susceptibilidad[i][j] = (1 - S_frozen[i][j])/temp[i];
  }
  array_print(Susceptibilidad, "Susceptibilidad.dat");
  
  time(&end);
  cout<<difftime(end,start);
  
  


  return 0;
}
