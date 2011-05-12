#include "sistema.h"
#include <cmath>

/*Constructor:
Dimensiona y encera a los vectores del sistema. Llena sus datos iniciales */
Sistema::Sistema(unsigned int lado, unsigned int dimension, bool polarizar)
{
  // Dimensionado de arreglos caraterísticos del sistema
  sigma.resize(pow(lado,dimension));
  
  sum_sigma.resize(sigma.size());
  for(unsigned int i = 0; i < sum_sigma.size(); i++)
    sum_sigma[i]=0;
  
  G.resize(sigma.size());
  J.resize(sigma.size());
  unsigned int vecinos = 2*dimension;
  for(unsigned int i = 0; i < J.size(); i++){
    G[i].resize(vecinos);
    J[i].resize(vecinos);
  }
  
  // Inicializa al sistema, llenado de datos
  //Crear arreglos auxiliares
  
  
  
  
}

/*Destructor:
libera la memoria asignada a los vectores del sistema*/
Sistema::~Sistema()
{
  J.clear();
  sum_sigma.clear();
  sigma.clear();
}
