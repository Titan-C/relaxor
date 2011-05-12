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

#define _pi atan(1)*4
#define r_max 2

using namespace std;
// Aplica las condiciones de borde toroidales
void condborde ( vector <double>& R, int L){
  double bL=L/2.0;
  for(unsigned int i=0; i<R.size();i++){
    if (R[i]>bL)	R[i]-=L;
    else if (R[i]<-bL)	R[i]+=L;
  }
}
//Producto interno
double dot(const vector<double>& a, const vector<double>& b){
  double A=0;
  for(unsigned int i=0;i<a.size();i++)
    A+=a[i]*b[i];
  return A;
}
// Imprime datos de variables double
void out(double value, string ARCHIVO){
  fstream file;
  file.open(ARCHIVO.c_str(), fstream::out | fstream::app);
  file<<value<<endl;
  file.close();
}
// Imprime datos de los arreglos vectoriales
void array_print(const vector<int>& V, string ARCHIVO) {
  fstream file;
  file.open (ARCHIVO.c_str(), fstream::out | fstream::app);
  for(unsigned int i = 0; i<V.size(); i++)
    file<<V[i]<<"\t";
  file<<endl;
  file.close();
}
// Imprime datos de los arreglos matricales
void array_print(const vector< vector<int> >& M, string ARCHIVO, bool p=true) {
  fstream file;
  file.open (ARCHIVO.c_str(), fstream::out);
  for (unsigned int i=0;i<M.size();i++) {
    for (unsigned int j=0;j<M[i].size();j++) {
      file<<M[i][j]<<"\t";
    }
    if (p); file<<endl;
  }
  file.close();
}
// Imprime datos de los arreglos matricales
void array_print(const vector< vector<double> >& M, string ARCHIVO, bool p=true) {
  fstream file;
  file.open (ARCHIVO.c_str(), fstream::out);
  for (unsigned int i=0;i<M.size();i++) {
    for (unsigned int j=0;j<M[i].size();j++) {
      file<<M[i][j]<<"\t";
    }
    if (p); file<<endl;
  }
  file.close();
}
//limpia los archivos
void file_wipe(string ARCHIVO){
  fstream file;
  file.open (ARCHIVO.c_str(), fstream::out);
  file.close();
}
void import_data(vector < vector < double > >& M, string ARCHIVO, unsigned int filas, unsigned int columnas) {
  fstream file;
  file.open (ARCHIVO.c_str(), fstream::in);
  M.resize(filas);
  for(unsigned int i=0;i<filas;i++){
    M[i].resize(columnas);
    for(unsigned int j=0;j<columnas;j++)
      file>>M[i][j];
  }
  file.close();
}
//Inicializa el sistema
void init ( vector <int>& sigma, vector< vector<double> >& J, vector< vector<int> >& vecinos,
	    int L, gsl_rng * rng){
  unsigned int PNR=sigma.size();
  vector< vector<double> > MAT, mu, R;
  MAT.resize(PNR);
  mu.resize(PNR);
  R.resize(PNR);
  
  //Creacion de las Nanoregiones Polares(PNR)
  for(unsigned int i=0; i<MAT.size();i++){
    MAT[i].resize(3);
    
    MAT[i][0] = r_max*gsl_rng_uniform(rng);// 
    MAT[i][1] = _pi/2*gsl_rng_uniform(rng);//Angulo theta caida desde z
    MAT[i][2] = 2*_pi*gsl_rng_uniform(rng);//Angulo phi plano xy
  }
  //imprimir datos de PNR
  array_print(MAT,"PNR_Material.dat");
  
  /*Obtención de los momentos dipolares eléctricos de cada PNR,
    la matriz de spines dipolares, posiciones de cada PNR*/
  int ind_xy;
  int L2 = L*L;
  for(unsigned int i=0; i<PNR;i++){
    mu[i].resize(3);
    mu[i][0] = sin(MAT[i][1])*cos(MAT[i][2]);
    mu[i][1] = sin(MAT[i][1])*sin(MAT[i][2]);
    mu[i][2] = cos(MAT[i][1]);
    
    sigma[i] = ( mu[i][2] > 0 ) ? 1 : -1;
    ind_xy = i%(L2);
    //Vector posición i-ésima PNR
    R[i].resize(3);
    R[i][0] = ind_xy % L; R[i][1] = ind_xy / L; R[i][2] = i/(L2);
    
    /*Encontrar índices de los primeros vecinos.
      solo existen 6: arriba y abajo(+z, -z), derecha e izquierda(+y, -y), adelante y atraz(+x, -x).
      También debo aplicar las condiciones de borde en este caso */
    vecinos[i][0] = (R[i][2] == L-1 )	?i - (L-1)*L2	:i + L2;//arriba
    vecinos[i][1] = (R[i][2] == 0 )	?i + (L-1)*L2	:i - L2;//abajo
    vecinos[i][2] = (R[i][1] == L-1 )	?i - (L-1)*L	:i + L;//derecha
    vecinos[i][3] = (R[i][1] == 0)	?i + (L-1)*L	:i - L;//izquierda
    vecinos[i][4] = (R[i][0] == L-1 )	?i - L+1	:i + 1;//adelante
    vecinos[i][5] = (R[i][0] == 0 )	?i + L-1	:i - 1;//atraz
  }
  
  //imprimir datos de mu
  array_print(mu, "mu_PNR.dat");
  array_print(R, "pos.dat");
  array_print(vecinos, "vecinos.dat");
  
  //Cálculo de las energías de intercambio de las PNR
  vector <double> Delta_R;
  vector < vector<double> > Jinter;
  double r, Jex;
  Delta_R.resize(3);
  Jinter.resize(PNR);
  
  for(unsigned int i=0; i<PNR; i++){
    Jinter[i].resize(PNR);
    for(unsigned int j=i+1; j<PNR; j++){
      //Calcular el vector diferencia entre 2 celdas y aplicar las condiciones de borde toroidales
      Jex=1.0;
      for(unsigned int k = 0; k<3; k++)
	Delta_R[k] = R[j][k] - R[i][k];
      condborde(Delta_R,L);
      r=sqrt(dot(Delta_R,Delta_R));
      
      Jex*=(3*dot(mu[i],Delta_R)*dot(mu[j],Delta_R)/(r*r*r*r*r) - dot(mu[i],mu[j])/(r*r*r));
      Jinter[i][j] = Jex/2;
    }
  }
    
  //Copias los datos de la matriz triangular superior a la parte inferior
  for(unsigned int i=0; i<PNR; i++){
    for(unsigned int j=i+1; j<PNR; j++)
      Jinter[j][i] = Jinter[i][j];
  }
  //imprime el arreglo de la matriz de interacción
  array_print(Jinter, "Matriz de intercambio.dat");
  
  /*Elaboración del arreglo de Energías de interación para los primeros vecinos
    solo existen 6: arriba y abajo(+z, -z), derecha e izquierda(+y, -y), adelante y atraz(+x, -x).
    También debo aplicar las condiciones de borde en este caso */
  for(unsigned int i=0; i<PNR; i++){
    for(unsigned int j=0; j<J[i].size(); j++)
      J[i][j] = Jinter[i][vecinos[i][j]];
  }
  //imprime el arreglo de energia de intercambio entre primeros vecinos
  array_print(J, "Jvecinos.dat");
}
//Calcula la Energía total del sistema
double total_e( const vector <int>& sigma, const vector< vector<double> >& J, const vector< vector<int> >& vecinos, double E=0.1){
  double EJ = 0;
  for(unsigned int i = 0; i<sigma.size(); i++){
    for(unsigned int j = 0; j<J[i].size(); j++){
      EJ-=J[i][j]*sigma[i]*sigma[vecinos[i][j]];
    }
  }
  return EJ;
}
//Calcula la variación en la energía del sistema debído al cambio de un spin dipolar
double delta_e(const vector<int>& sigma, const vector < vector <double> >& J,
	       const vector< vector<int> >& vecinos, int idflip){
  double dE = 0;
  for(unsigned int i = 0; i<J[idflip].size(); i++){
/*    dE+=J[idflip][i]*sigma[idflip]*sigma[i];
la variación de energía es con signo positivo, pero si sigma[idflip] cambia de orientación, si
no realizo ese cambio, entonces debo incluirlo como un factor de -1 en el cálculo o como un resta */
      dE+=J[idflip][i]*sigma[idflip]*sigma[vecinos[idflip][i]];
  }
  //cout<<dE<<"\t";
  return 2*dE;
}
//realiza el cambio del spin dipolar en una ubicación
void flip(vector<int>& sigma, const vector< vector <double> >& J,
	  const vector< vector<int> >& vecinos, int idflip, double T, gsl_rng * rng){
  double dE = delta_e(sigma, J,vecinos, idflip);
  if ( dE < 0){
    sigma[idflip]*=-1;
    //cout<<"flip directo"<<endl;
  }
  else if ( exp(-dE/T) >= gsl_rng_uniform(rng) ){
    sigma[idflip]*=-1;
    //cout<<"flip metropolis"<<endl;
  }
}
  

int main(int argc, char **argv) {
  vector <int> sigma; //arreglo de spines dipolares
  vector <int> avr_sigma; //suma de spine para calular el promedio de los spines congelados
  vector< vector<int> > vecinos;
  vector< vector <double> > J; //Arreglo de Energía de intercambio
  gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);
  
  unsigned int L=8, Niter=400, numexps = 1;
  double T=7;/*
  cin>>L;
  cin>> Niter;
  cin>>T;*/
  srand( time(NULL) );
  gsl_rng_set(r, time(NULL) );
  sigma.resize(L*L*L);
  avr_sigma.resize(sigma.size());
  for(unsigned int i = 0; i < avr_sigma.size(); i++)
    avr_sigma[i]=0;
  
  vecinos.resize(sigma.size());
  J.resize(sigma.size());
  for(unsigned int i=0; i < J.size() ;i++){
    vecinos[i].resize(6);
    J[i].resize(6);
  }
  
  //iniciar sistema
  time_t start, end;
  time(&start);
  init(sigma, J, vecinos, L, r);
  array_print(sigma, "sigmas.dat");
  //vaciar archivo de datos en cada ejecución
  file_wipe("avr_sigmas.dat");
  file_wipe("energy_log.dat");
  cout<<total_e(sigma,J,vecinos)<<endl;
  unsigned int mediciones = 1;

  do{
    //ejecutar una corrida en el tiempo correspondiente a Niter  
    for(unsigned int i = 0; i < Niter; i++){
      out(total_e(sigma,J,vecinos), "energy_log.dat");
      //1MCS probar cambiar la orientacion de todos los spines del arreglo
      for(unsigned int idflip = 0; idflip < sigma.size(); idflip++)
	flip(sigma, J, vecinos, idflip, T, r);
      
      //tomar datos despues de 3 periodos de equilibración
      if ( numexps % 4 == 0 ){
	for(unsigned int k = 0; k<sigma.size(); k++)
	  avr_sigma[k]+=sigma[k];
      }      
    }
    //Escribir datos del experimento
    if ( numexps % 4 == 0 ){
      array_print(avr_sigma,"avr_sigmas.dat");
      T-=0.1;
      mediciones++;
      //reiniciar tomador de datos
      for(unsigned int i = 0; i<sigma.size(); i++)
	avr_sigma[i]=0;
    }
    
    numexps++;      
  }while(T>0);

  //procesar datos polarizacion congelada
  vector< vector<double> > promedio_sigmas;
  vector< vector<double> > S_frozen, Susceptibilidad;
  vector<double> temp;
  unsigned int lentos=4;
  S_frozen.resize(mediciones-1);
  for(unsigned int i=0; i< S_frozen.size(); i++){
    S_frozen[i].resize(lentos);
    for(unsigned int j=0; j< lentos; j++)
      S_frozen[i][j] = 0;
  }

  import_data(promedio_sigmas, "avr_sigmas.dat", (mediciones-1), sigma.size());
  
  for(unsigned int i=0 ; i<promedio_sigmas.size(); i++){
    for(unsigned int j=0; j<promedio_sigmas[i].size(); j++){
      promedio_sigmas[i][j] = abs(promedio_sigmas[i][j]/Niter);
      for(unsigned int k=0;k<lentos;k++){
	if (promedio_sigmas[i][j] >= (1-0.1*k) )
	  S_frozen[i][k]+=(double) 1/sigma.size();
      }
    }    
  }
  array_print(S_frozen, "Congelamiento.dat");
  
  //Arreglo de temperaturas
  temp.resize(mediciones);
  for(unsigned int i = 0; i<mediciones; i++)
    temp[i] = T + (mediciones-i)*0.1;
  
  //Encontrar la susceptibilidad del material
  Susceptibilidad.resize(mediciones);
  for(unsigned int i=0; i<S_frozen.size(); i++){
    Susceptibilidad[i].resize(lentos);
    for(unsigned int j=0; j<lentos; j++)
      Susceptibilidad[i][j] = (1 - S_frozen[i][j])/temp[i];
  }
  array_print(Susceptibilidad, "Susceptibilidad.dat");
  
  time(&end);
  cout<<difftime(end,start);
  
  


  return 0;
}
