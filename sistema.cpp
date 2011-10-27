#include <cmath>
#include <ctime>
#include <fstream>
#include <sstream>
#include "sistema.h"
#include "impresor.h"
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_randist.h>

#define _2pi 8*atan(1)

/*Constructor:
Dimensiona y encera a los vectores del sistema. Llena sus datos iniciales */
Sistema::Sistema(unsigned int lado,
		 gsl_rng* rng,
		 double Delta_J,
		 unsigned int dim,
		 bool polarizar){
  dimension = dim;
  L = lado;
  // Dimensionado de arreglos caraterísticos del sistema
  sigma.resize(pow(lado,dimension));
  mu_E.resize(sigma.size());

  J.resize(sigma.size());
  G.resize(sigma.size());
  unsigned int vecinos = 2*dimension;
  for(unsigned int i = 0; i < J.size(); i++){
    J[i].resize(vecinos);
    G[i].resize(vecinos);
  }
  //Generar configuración espacial de PNR
  set_space_config();

  // Inicializa al sistema, llenado de datos
  DeltaJ = init(rng, Delta_J, polarizar);  
}

/*Destructor:
libera la memoria asignada a los vectores del sistema*/
Sistema::~Sistema(){
  for(unsigned int i = 0; i < J.size(); i++){
    J[i].clear();
    G[i].clear();
  }
  J.clear();
  G.clear();
  mu_E.clear();
  sigma.clear();
}

double Sistema::set_pol(gsl_rng* rng, bool polarizar){
  if (polarizar)
    sigma.assign(sigma.size(),1);
  else{
    for(unsigned int i=0; i<sigma.size(); i++)
      sigma[i] = (gsl_rng_uniform(rng)-0.5 > 0)? 1:-1;
  }
  //Calcular la polarización inicial
  return norm_pol();
}
void Sistema::set_mu(gsl_rng* rng){
  for(unsigned int i=0; i<sigma.size(); i++)
    mu_E[i]  = gsl_rng_uniform(rng);
}

void Sistema::set_space_config(){
  unsigned int ind_xy, L2=L*L;
  std::vector< std::vector<unsigned int> > R;
  R.resize(sigma.size());

  for(unsigned int i=0; i<R.size(); i++){
    // Coeficientes vector posición i-ésima PNR
    ind_xy = i % L2;
    R[i].resize(dimension);
    R[i][0] = ind_xy % L;
    R[i][1] = ind_xy / L;
    R[i][2] = i / L2;

    /*Encontrar índices de los primeros vecinos.
     *solo existen 6: arriba y abajo(+z, -z), derecha e izquierda(+y, -y), adelante y atraz(+x, -x).
     *También debo aplicar las condiciones de borde en este caso*/
    G[i][0] = (R[i][2] == L-1 )	?i - (L-1)*L2	:i + L2;//arriba
    G[i][1] = (R[i][2] == 0 )	?i + (L-1)*L2	:i - L2;//abajo
    G[i][2] = (R[i][1] == L-1 )	?i - (L-1)*L	:i + L;//derecha
    G[i][3] = (R[i][1] == 0)	?i + (L-1)*L	:i - L;//izquierda
    G[i][4] = (R[i][0] == L-1 )	?i - L+1	:i + 1;//adelante
    G[i][5] = (R[i][0] == 0 )	?i + L-1	:i - 1;//atraz    
  }
  //liberar mem
  for(unsigned int i=0; i<R.size();i++)
    R[i].clear();
  R.clear();
}

double Sistema::Jex(gsl_rng* rng, double Delta_J){
  std::vector< std::vector<double> > Jinter;
  Jinter.resize(sigma.size());
  /*Genera las matriz triangular superior de las
   * energías de intecambio segun una distribución
   * Gaussiana*/
  for(unsigned int i = 0; i<Jinter.size(); i++){
    Jinter[i].resize(sigma.size());
    for(unsigned int j = i+1; j<Jinter[i].size(); j++)
      Jinter[i][j] = gsl_ran_gaussian(rng,Delta_J);
  }
  //Completa la parte inferior de la matriz de intercambio
  for(unsigned int i = 0; i<Jinter.size(); i++){
    for(unsigned int j = i+1; j<Jinter.size(); j++)
      Jinter[j][i] = Jinter[i][j];
  }
  // Elabora el arreglo de interacción de primeros vecinos
  for(unsigned int i=0; i<J.size(); i++){
    for(unsigned int j=0; j<J[i].size(); j++)
      J[i][j] = Jinter[i][G[i][j]];
  }
  //liberar mem
  for(unsigned int i=0; i<Jinter.size();i++)
    Jinter[i].clear();
  Jinter.clear();

  return Delta_J;
}

double Sistema::init(gsl_rng* rng, double Delta_J, bool polarizar){
  // Genera las energías de intercambio de las PNR
  std::cout<<"Des stan Total= "<<Jex(rng,Delta_J)<<"\n";
  
  //Momentos dipolares y polarización
  set_mu(rng);//Componente del momento dipolar eléctrico en eje principal
  std::cout<<"Polarización inicial="<<set_pol(rng, polarizar)<<"\n";
  
  return Delta_J;
}

double Sistema::total_E(double E){
  double Hamil = 0;
  for(unsigned int i = 0; i < G.size(); i++){
    for(unsigned int j = 0; j < G[i].size(); j++)
      Hamil -= J[i][j]*sigma[i]*sigma[G[i][j]];
    Hamil -= E*mu_E[i]*sigma[i];
  }
  return Hamil;
}

double Sistema::delta_E(unsigned int idflip, double E){
  double dHamil = 0;
  for(unsigned int i = 0; i<G[idflip].size(); i++)
    dHamil += J[idflip][i]*sigma[idflip]*sigma[G[idflip][i]];
  dHamil *=4;
  dHamil += 2*E*mu_E[idflip]*sigma[idflip];
  return dHamil;
}

double Sistema::norm_pol(){
  unsigned int N=sigma.size();
  double p=0;
  for(unsigned int i=0; i<N; i++)
    p += mu_E[i]*sigma[i];

  return (double) p / N;
}

int Sistema::experimento(double T, double E, unsigned int tau, unsigned int Niter,
			 bool grabar, gsl_rng* rng, std::string id_proc){
  //vector historial de polarización por experimento
  std::vector<double> pol_mag;
  pol_mag.resize(Niter);
  
  //vector oscilación del campo alterno para un periodo
  std::vector<double> field;
  field = waves(tau,tau,E,true);
  
  /*Simulación del experimento en el número de iteraciones dadas*/

  unsigned int periods = Niter/tau;
  unsigned int step = 0;
  for(unsigned int i = 0; i<periods; i++){
    for(unsigned int j = 0; j< tau; j++){
      //       out(total_E(E), "energy_log.dat");
      /*Realiza el cambio del spin dipolar en una ubicación dada*/
      for(unsigned int idflip = 0; idflip < sigma.size(); idflip++){
	double dH = delta_E(idflip, field[j]);
	if ( dH < 0)
	  sigma[idflip] *= -1;
	else if ( exp(-dH/T) >= gsl_rng_uniform(rng) )
	  sigma[idflip] *= -1;
      }
      
      if (grabar)
	pol_mag[step] = norm_pol();
      step++;
    }
  }
  /* Guardar los datos de polarización en binario */
  if (grabar)
    array_print_bin(pol_mag,"log_pol_"+id_proc+".dat");
  
  pol_mag.clear();
  field.clear();
  
  return 1;
}

double stan_dev(const std::vector< std::vector<double> >& M){
  //Calcular la desviación estandar de una matriz
  unsigned int celdas, columnas;
  columnas = M[1].size();
  celdas = M.size() * columnas;
  double * Aij;
  Aij = new double [celdas];
  for(unsigned int i = 0 ; i<M.size(); i++){
    for(unsigned int j = 0; j<columnas; j++)
      Aij[i*columnas + j] = M[i][j];
  }
  return gsl_stats_sd (Aij, 1, celdas);
  delete[] Aij;
}

std::vector<double> step2vec(double unidad, double v_top, double dv, bool grow){
  std::vector<double> temp;
  temp.resize( (int) ceil(1/dv*v_top));
  for(unsigned int i = 0; i<temp.size() ;i++){
    if (grow)
      temp[i] = (i+1)*dv*unidad;
    else
      temp[i] = (v_top - i*dv)*unidad;
  }
  
  return temp;
}
std::vector<double> str2vec(double unidad, std::string magnitudes){
  std::istringstream data(magnitudes);
  std::vector<double> data_array;
  double num;
  while(!data.eof()){
    data >> num;
    data_array.push_back(num*unidad);
  }
  return data_array;
}
void calc_sus(unsigned int numexps, unsigned int tau, unsigned int Niter, double unidad,
	      const std::vector<double>& x_array, const std::vector<double>& campo, std::string id_proc){
  //Vectores de información
  std::vector<double> pol_hist;
  pol_hist.resize(Niter);
  std::string name ="log_pol_"+id_proc+".dat";
  std::ifstream file (name.c_str());
  
  /*Generar arreglo de peso sin, cos para la integral*/
  std::vector<double> cos_wave, sin_wave;
  cos_wave = waves(Niter,tau,1.0,true);
  sin_wave = waves(Niter,tau,1.0,false);
  
  double * Freal = new double [numexps*x_array.size()];
  double * Fimag = new double [numexps*x_array.size()];
  bool fieldvec = (campo.size()>1) ? true : false;
  for(unsigned int n = 0; n < numexps; n++){    
    for(unsigned int x=0;x<x_array.size();x++){
      file.read((char * )&pol_hist[0],Niter*sizeof(double));
      
      /*Integración por Simpson*/
      double Int_cos=simpson_int(pol_hist,cos_wave);
      double Int_sin=simpson_int(pol_hist,sin_wave);
      
      int ind=(fieldvec) ? x : 0;
      Freal[n*x_array.size()+x]=Int_cos/Niter/campo[ind];
      Fimag[n*x_array.size()+x]=Int_sin/Niter/campo[ind];
    }
  }
  //liberar memoria
  cos_wave.clear();
  sin_wave.clear();
  pol_hist.clear();
  file.close();
  
  /*Calcular susceptibildad más error*/
  std::vector< std::vector<double> > X_mat;
  X_mat.resize(x_array.size());
  for(unsigned int x=0;x<x_array.size();x++){
    X_mat[x].resize(5);
    X_mat[x][0]=x_array[x]/unidad;
  }
  double * data_arrayr = new double [numexps];
  double * data_arrayi = new double [numexps];
  for(unsigned int x=0;x<x_array.size();x++){
    
    for(unsigned int n=0;n<numexps;n++){
      data_arrayr[n]=Freal[n*x_array.size()+x];
      data_arrayi[n]=Fimag[n*x_array.size()+x];
    }
    X_mat[x][1]=gsl_stats_mean(data_arrayr,1,numexps);
    X_mat[x][2]=gsl_stats_sd_m(data_arrayr,1,numexps,X_mat[x][1]);
    X_mat[x][3]=gsl_stats_mean(data_arrayi,1,numexps);
    X_mat[x][4]=gsl_stats_sd_m(data_arrayi,1,numexps,X_mat[x][3]);
  }
  array_print(X_mat, "susceptibilidad_"+id_proc+".dat");
  
  //liberar memoria
  for(unsigned int i=0; i<X_mat.size();i++)
    X_mat[i].clear();
  X_mat.clear();  
  delete[] data_arrayr;
  delete[] data_arrayi;
  delete[] Freal;
  delete[] Fimag;
}

std::vector<double> waves(unsigned int length, unsigned int tau, double amplitude, bool cossin){
  std::vector<double> wave;
  wave.resize(length);
  if (cossin) {
    for(unsigned int i=0; i<tau; i++)
      wave[i]=amplitude*std::cos(_2pi*i/tau);
  }
  else {
    for(unsigned int i=0; i<tau; i++)
      wave[i]=amplitude*std::sin(_2pi*i/tau);
  }
  
  unsigned int periods = length/tau;
  for(unsigned int i=1; i<periods;i++){
    for(unsigned int j = 0; j<tau ;j++)
      wave[i*tau+j]=wave[j];
  }
  return wave;
}

double simpson_int(const std::vector<double>& f_array, const std::vector<double>& weight){
  unsigned int length=f_array.size();
  
  double Integral = f_array[0]*weight[0];
  
  for(unsigned int i=1; i<length-1; i+=2)
    Integral+=4*f_array[i]*weight[i];
  
  for(unsigned int i=2; i<length; i+=2)
    Integral+=2*f_array[i]*weight[i];
  
  length--;
  Integral+=f_array[length]*weight[length];
  
  return Integral/3;
}

void eval_pol(unsigned int Niter, unsigned int numexps, double unidad, const std::vector<double>& x_array, std::string id_proc, bool absolut) {
  
  std::string name = "log_pol_"+id_proc+".dat";
  std::ifstream file(name.c_str());
  double * pol_hist = new double [Niter];
  /* Calcular la polarización media y desviación estandar para el material en cada experimento y para cada temperatura. */  
  std::vector< std::vector<double> > pol_stats;
  pol_stats.resize(x_array.size());
  for(unsigned int n=0; n<numexps; n++){
    for(unsigned int T=0; T< x_array.size(); T++){
      file.read((char *)&pol_hist[0],Niter*sizeof(double));
      pol_stats[T].resize(2*numexps);
      pol_stats[T][2*n] = gsl_stats_mean(pol_hist,1,Niter);
      pol_stats[T][2*n+1] = gsl_stats_sd_m(pol_hist,1,Niter, pol_stats[T][2*n]);
    }
  }
  delete[] pol_hist;
  file.close();
  array_print(pol_stats,"avg_pol"+id_proc+".dat");
  

  //Polarización, o polarización absoluta y desviación estandar
  double * data_array = new double [numexps];
  std::vector< std::vector<double> > pol_final;
  pol_final.resize(x_array.size());
  for(unsigned int T=0;T< x_array.size(); T++){
    
    pol_final[T].resize(3);
    pol_final[T][0]=x_array[T];
    
    for(unsigned int n=0;n<numexps;n++)
      data_array[n]=(absolut) ? std::abs(pol_stats[T][2*n]) : pol_stats[T][2*n];
    pol_final[T][1]=gsl_stats_mean(data_array,1,numexps);
    
    for(unsigned int n=0;n<numexps;n++)
      pol_final[T][2]+=pol_stats[T][2*n+1]*pol_stats[T][2*n+1];
    pol_final[T][2]=sqrt(pol_final[T][2]/numexps);
  }
  delete[] data_array;
  array_print(pol_final,"pol_"+id_proc+".dat");

  /*Liberar memoria*/
  for(unsigned int i=0;i<x_array.size();i++){
    pol_stats[i].clear();
    pol_final[i].clear();
  }
  pol_stats.clear();
  pol_final.clear();
}

