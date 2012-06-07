#include <cmath>
#include <ctime>
#include <cstdio>
#include <fstream>
#include <sstream>
#include "sistema.h"
#include "impresor.h"
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_randist.h>
#include <sys/types.h>
#include <sys/stat.h>

#define _2pi 8*atan(1)

/*Constructor:
Dimensiona y encera a los vectores del sistema. Llena sus datos iniciales */
Sistema::Sistema(unsigned int lado,
		 unsigned int dim,
		 bool polarizar){
  dimension = dim;
  L = lado;
  PNR = pow(lado,dimension);
  rng = gsl_rng_alloc (gsl_rng_taus);
  // Dimensionado de arreglos caraterísticos del sistema
  sigma = new int8_t [PNR];
  mu_E = new double [PNR];

  vecinos = 2*dimension;
  J = new double *[PNR];
  G = new unsigned int *[PNR];
  for(unsigned int i = 0; i < PNR; i++){
    J[i] = new double [vecinos];
    G[i] = new unsigned int [vecinos];
    for(unsigned int j=0;j<vecinos;j++)
      J[i][j] = -1000;
  }
  //Generar configuración espacial de PNR
  set_space_config();
}

/*Destructor:
libera la memoria asignada a los vectores del sistema*/
Sistema::~Sistema(){
  gsl_rng_free (rng);
  for(unsigned int i=0; i<PNR; i++){
    delete[] J[i];
    delete[] G[i];
  }
  delete[] J;
  delete[] G;
  delete[] mu_E;
  delete[] sigma;
}

void Sistema::set_space_config(){
  unsigned int ind_xy, L2=L*L;
  std::vector< std::vector<unsigned int> > R;
  R.resize(PNR);

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

double Sistema::Jex(){
//Genera directamente el arreglo de energías de intercambio
//para la simulación ahorrando memoria
  for(unsigned int i=0;i<PNR;i++){
    for(unsigned int j=0;j<vecinos;j++){
      if (J[i][j] == -1000){
        for(unsigned int k=0; k<vecinos;k++){
	  if (G[ G[i][j] ][k] == i){
	    if (J [G[i][j] ][k] == -1000)
	      J[ G[i][j] ][k] = gsl_ran_gaussian(rng,1)+rho;
	    J[i][j] = J[ G[i][j] ][k];
	  }
	}}}}
  return stan_dev(J,PNR,vecinos);
}

double Sistema::set_pol(bool polarizar){
  if (polarizar)
    for(unsigned int i=0; i<PNR; i++) sigma[i] = 1;
  else{
    for(unsigned int i=0; i<PNR; i++)
      sigma[i] = (gsl_rng_uniform(rng)-0.5 > 0)? 1:-1;
  }

  return norm_pol();
}

double Sistema::set_mu(bool polarizar){
  for(unsigned int i=0; i<PNR; i++)
    mu_E[i]  = gsl_rng_uniform(rng);
  array_print(mu_E,PNR,"mu.dat");
  return set_pol(polarizar);
}

void Sistema::init(double p, bool polarizar, bool write){
  //Cambiar la semilla del generador de números aleatorios dentro de la clase
  gsl_rng_set(rng, std::time(NULL) );
  //Cambia las propidades del sistema
  rho = p;
  // Genera las energías de intercambio de las PNR
  double Delta_J=Jex();
  //Momentos dipolares y polarización
  double pol=set_mu(polarizar);
  
  if (write){
    std::cout<<"Desvición Estandar Total= "<<Delta_J<<"\n";
    std::cout<<"Polarización inicial="<<pol<<"\n";
  }
}

double Sistema::total_E(double E){
  double Hamil = 0;
  for(unsigned int i = 0; i < PNR; i++){
    for(unsigned int j = 0; j < vecinos; j++)
      Hamil -= J[i][j]*sigma[i]*sigma[G[i][j]];
    Hamil -= E*mu_E[i]*sigma[i];
  }
  return Hamil;
}

double Sistema::delta_E(unsigned int idflip, double E){
  double dHamil = 0;
  for(unsigned int i = 0; i<vecinos; i++)
    dHamil += J[idflip][i]*sigma[idflip]*sigma[G[idflip][i]];
  dHamil *=4;
  dHamil += 2*E*mu_E[idflip]*sigma[idflip];

  return dHamil;
}

double Sistema::norm_pol(){
  double P=0;
  for(unsigned int i=0; i<PNR; i++)
    P += mu_E[i]*sigma[i];

  return (double) P / PNR;
}

void Sistema::experimento(double T, double E, unsigned int tau, unsigned int Niter,
			  bool grabar, std::string id_proc){
  //vector historial de polarización por experimento
  std::vector<double> pol_mag;
  pol_mag.resize(Niter);
  
  //vector oscilación del campo alterno para un periodo
  std::vector<double> field;
  double phase = 0;
  if (!grabar)
    phase = _2pi*Niter/tau;

  field = cosarray(Niter,tau,E,phase);

  /*Simulación del experimento en el número de iteraciones dadas*/
  for(unsigned int i = 0; i< Niter; i++){
    /*Realiza el cambio del spin dipolar en una ubicación dada*/
    for(unsigned int idflip = 0; idflip < PNR; idflip++){
      double dH = delta_E(idflip, field[i]);
      if ( dH < 0 || exp(-dH/T) >= gsl_rng_uniform(rng) )
	sigma[idflip] *= -1;
    }
    if (grabar){
      pol_mag[i] = norm_pol();
      array_print_bin(ret_sigarr(),PNR,"log_sigma_"+id_proc+".dat");
    }
  }
  

  /* Guardar los datos de polarización en binario */
  if (grabar)
    array_print_bin(pol_mag,"log_pol_"+id_proc+".dat");

  pol_mag.clear();
  field.clear();
}

void Gen_exp(unsigned int L, unsigned int numexps, std::vector<double> rho, std::vector<double>& Tdat,
	     std::vector<double>& Fields, std::vector<double> tau, std::string Exp_ID)
{
  Sistema relaxor(L);
  unsigned int Equi_iter=350;
  //   Código para bajar (variar) la temperatura a campos fijos
  if (Exp_ID == "cool" || Exp_ID == "heat"){
    for(unsigned int p=0; p<rho.size(); p++){
      std::vector<double> Thermostat;
      Thermostat= thermostat(rho.size(), p, rho[p], Tdat[0], Tdat[1]);
      for(unsigned int t=0; t< tau.size() ; t++){
	unsigned int Exp_iter = stepEstimator(3000,tau[t],2);
	
	for(unsigned int E=0; E<Fields.size(); E++){
	  std::ostringstream id_proc;
	  id_proc<<Exp_ID<<"_p"<<rho[p]<<"_E"<<Fields[E]<<"_t"<<tau[t]<<"_L"<<L<<"_n"<<numexps;
	  id_proc<<"_Ti"<<Thermostat[0]<<"Tf"<<Tdat[1]<<"dT"<<Tdat[0]<<"_X"<<Exp_iter<<"_Q"<<Equi_iter;
	  
	  clock_t cl_start = clock();
	  unsigned int sim_size = sizeof(double)*Exp_iter*Thermostat.size()*numexps;
	  if (needSimulation(id_proc.str(), sim_size)) {
	    for(unsigned int n=0; n<numexps; n++){
	      relaxor.init(rho[p],false);
	      for(unsigned int T=0; T<Thermostat.size(); T++){
		relaxor.experimento(Thermostat[T],Fields[E],tau[t], Equi_iter,false, id_proc.str());
		relaxor.experimento(Thermostat[T],Fields[E],tau[t], Exp_iter,true, id_proc.str());
	      }}
	  }
	  proces_data(Thermostat,Fields[E],tau[t],numexps,relaxor.return_PNR(), rho[p],Exp_iter,id_proc.str());
	  std::cout<<id_proc.str()<<":"<<clock()-cl_start<<"\n";
	}}
    }
  }
//   Código para variar el campo a temperaturas fijas, desactualizado
//   else {
//     for(unsigned int t=0; t< tau.size() ; t++){
//       for(unsigned int T=0; T<Temps.size(); T++){
// 	std::ostringstream id_proc;
// 	id_proc<<Exp_ID<<"_p"<<p<<"_T"<<Temps[T]<<"_t"<<tau[t];
// 	for(unsigned int n=0; n<numexps; n++){
// 	  relaxor.init(p,false);
// 	  for(unsigned int E=0; E<Fields.size(); E++){
// 	    relaxor.experimento(Temps[T],Fields[E],tau[t], Equi_iter,false, id_proc.str());
// 	    relaxor.experimento(Temps[T],Fields[E],tau[t], Exp_iter,true, id_proc.str());
// 	  }}}}}
}

void proces_data(std::vector< double >& Temps, double Field,
		 unsigned int tau, unsigned int numexps, unsigned int PNR,
		 double p, unsigned int Niter, std::string id_proc){
  //procesar los datos para casos de variación de temperatura
    std::vector<double> pol_stats, pol_int_avg;
    pp_data(pol_stats,pol_int_avg,Temps.size(),numexps,tau,Niter,id_proc);
    std::vector<double> intfield (1,Field);
    eval_pol(pol_stats,numexps,Temps,id_proc,true);
    calc_sus(pol_int_avg,numexps,Temps,intfield,id_proc);
    eval_frozen(PNR,Niter, Temps, numexps, id_proc);
    intfield.clear();
    pol_int_avg.clear();
    pol_stats.clear();

//       else {
// 	id_proc<<Exp_ID<<"_p"<<p<<"_T"<<Temps[T]<<"_t"<<tau[t];
// 	pp_data(pol_stats,pol_int_avg,Fields.size(),numexps,tau[t],Niter,id_proc.str());
// 	eval_pol(pol_stats,numexps,Fields,id_proc.str(),(Exp_ID=="hist_loop") ? false : true);
// 	calc_sus(pol_int_avg,numexps,Fields,Fields,id_proc.str());
//       }

}

void pp_data(std::vector<double>& pol_stats, std::vector<double>& pol_int_avg, unsigned int data_length,
	      unsigned int numexps, unsigned int tau, unsigned int Niter, std::string id_proc){
  //Generar vectores de Datos
  double * pol_hist = new double[Niter];
  unsigned int dat_vec_size = data_length*numexps*2;
  pol_stats.resize(dat_vec_size);
  pol_int_avg.resize(dat_vec_size);
  
  /*Generar arreglo de peso sin, cos para la integral*/
  std::vector<double> cos_wave, sin_wave;
  cos_wave = cosarray(Niter,tau,1.0,0);
  sin_wave = cosarray(Niter,tau,1.0,_2pi/4);
  
  //Abrir Archivo y leer
  std::string name = "log_pol_"+id_proc+".dat";
  std::ifstream file(name.c_str());
  for(unsigned int ind=0; ind<dat_vec_size; ind+=2){
      file.read((char *)&pol_hist[0],Niter*sizeof(double));
      /*Calcular media y desviación estandar de polarización por
       numero y condiciones de experimento*/
      pol_stats[ind]=gsl_stats_mean(pol_hist,1,Niter);
      pol_stats[ind+1]=gsl_stats_sd_m(pol_hist,1,Niter, pol_stats[ind]);
      /*Integración por Simpson, para promedio pesado */
      pol_int_avg[ind]=simpson_int(pol_hist,cos_wave)/Niter;
      pol_int_avg[ind+1]=simpson_int(pol_hist,sin_wave)/Niter;
  }
  //liberar memoria
  cos_wave.clear();
  sin_wave.clear();
  delete[] pol_hist;
  file.close();
}
void eval_frozen(unsigned int PNR, unsigned int Niter, const std::vector<double>& Temps, unsigned int numexps, std::string id_proc){
  int8_t * sigma_hist = new int8_t[PNR];
  std::vector< double > sigmaTemp;
  sigmaTemp.assign(numexps*Temps.size()*PNR,0);
  
  //Abrir Archivo, leer guardar datos
  std::string name = "log_sigma_"+id_proc+".dat";
  std::ifstream file(name.c_str());
  for(unsigned int n=0; n<numexps ; n++){
    for(unsigned int T=0; T<Temps.size(); T++){
      for(unsigned int iter = 0; iter < Niter ; iter++){
	file.read((char *)&sigma_hist[0],PNR*sizeof(int8_t));
	for(unsigned int s = 0; s<PNR ; s++)
	  sigmaTemp[n*Temps.size()*PNR+T*PNR + s] += sigma_hist[s];
  }}}
  array_print(sigmaTemp, "sigmas_"+id_proc+".dat", PNR, Niter);
  delete[] sigma_hist;
  //Evaluar % congelamiento
  std::vector< std::vector<double> > Frozen;
  Frozen.resize(Temps.size());
  for(unsigned int T=0; T<Temps.size(); T++){
    Frozen[T].assign(4,0);
    Frozen[T][0]=Temps[T];
    for(unsigned int s = 0; s<PNR ; s++){
      double avgabssigma= std::abs(sigmaTemp[T*PNR + s]/Niter);
      if (avgabssigma > 0.9){
	Frozen[T][3]++;
	if (avgabssigma > 0.95){
	  Frozen[T][2]++;
	  if (avgabssigma > 0.999999)
	    Frozen[T][1]++;
	}}
    }
    for(unsigned int f=1; f<4;f++)
      Frozen[T][f]/=PNR;
  }
  array_print(Frozen,"frozen_"+id_proc+".dat");
  /*Liberar memoria*/
  sigmaTemp.clear();
  for(unsigned int T=0; T<Temps.size(); T++)
    Frozen[T].clear();
  Frozen.clear();
}
void eval_pol(const std::vector<double>& pol_stats, unsigned int numexps, const std::vector<double>& x_array, std::string id_proc, bool absolut) {

  //Polarización, o polarización absoluta y desviación estandar
  double * data_array = new double [numexps];
  double data_length = x_array.size();
  std::vector< std::vector<double> > pol_final;
  pol_final.resize(data_length);
  for(unsigned int x=0;x< data_length; x++){
    pol_final[x].resize(3);
    pol_final[x][0]=x_array[x];
    
    double pol_std=0;
    for(unsigned int n=0;n<numexps;n++){
      unsigned int ind = 2*(n*data_length+x);
      data_array[n]=(absolut) ? std::abs(pol_stats[ind]) : pol_stats[ind];
      pol_std+=pol_stats[ind+1]*pol_stats[ind+1];
    }
    pol_final[x][1]=gsl_stats_mean(data_array,1,numexps);
    pol_final[x][2]=sqrt(pol_std/numexps);
  }
  delete[] data_array;
  array_print(pol_final,"pol_"+id_proc+".dat");

  /*Liberar memoria*/
  for(unsigned int i=0;i<data_length;i++)
    pol_final[i].clear();
  pol_final.clear();
}

void calc_sus(const std::vector<double>& pol_int_avg, unsigned int numexps,
	      const std::vector<double>& x_array, const std::vector<double>& campo, std::string id_proc){

  /*Calcular susceptibildad más error*/
  double * data_arrayr = new double [numexps];
  double * data_arrayi = new double [numexps];
  bool fieldvec = (campo.size()>1) ? true : false;
  double data_length = x_array.size();
  std::vector< std::vector<double> > X_mat;
  X_mat.resize(data_length);
  for(unsigned int x=0;x<data_length;x++){
    X_mat[x].resize(7);
    X_mat[x][0]=x_array[x];
    
    for(unsigned int n=0;n<numexps;n++){
      unsigned int ind = 2*(n*data_length+x);
      unsigned int field_ind=(fieldvec) ? x : 0;
      data_arrayr[n]=pol_int_avg[ind]/campo[field_ind];
      data_arrayi[n]=pol_int_avg[ind+1]/campo[field_ind];
    }
    X_mat[x][1]=gsl_stats_mean(data_arrayr,1,numexps);
    X_mat[x][3]=gsl_stats_sd_m(data_arrayr,1,numexps,X_mat[x][1]);
    X_mat[x][2]=0.0001/X_mat[x][3]/X_mat[x][3];
    X_mat[x][4]=gsl_stats_mean(data_arrayi,1,numexps);
    X_mat[x][6]=gsl_stats_sd_m(data_arrayi,1,numexps,X_mat[x][4]);
    X_mat[x][5]=0.0001/X_mat[x][6]/X_mat[x][6];
  }
  array_print(X_mat, "sus_"+id_proc+".dat");

  //liberar memoria
  for(unsigned int i=0; i<X_mat.size();i++)
    X_mat[i].clear();
  X_mat.clear();  
  delete[] data_arrayr;
  delete[] data_arrayi;
}
std::vector< double > thermostat(unsigned int n, unsigned int i, double rho, double dT, double Tf){
  double Ti = (rho>0.5) ? 10*rho+2.5 : 8;
  double shift = (double) dT/(i+1);
  Ti+=shift;
  
  std::vector<double> Temparray;
  Temparray.clear();
  Temparray=step2vec(Ti,Tf,dT,Temparray);
  
  return Temparray;
}

std::vector<double> step2vec(double v_start, double v_end, double dv, std::vector<double> last, double unidad){
  while(v_start<=v_end) {
    last.push_back(v_start*unidad);
    v_start+=dv;
  }
  while (v_start>v_end) {
    last.push_back(v_start*unidad);
    v_start-=dv;
  }
  
  return last;
}

std::vector<double> loop2vec(double max, int divs, double unidad){
  std::vector<double> vec;
  double step = 1.0 / divs;
  unsigned int l=0;
  for(double x=step; x<=1 ;x+=step){
    vec.push_back(x*x*max*unidad);
    l++;
  }
  for(unsigned int i=l-1; i>0; i--)
    vec.push_back(vec[i-1]);

  unsigned int vec_size= vec.size();
  for(unsigned int i=0; i<vec_size; i++)
    vec.push_back(-1*vec[i]);

  for(unsigned int i=0; i<=l;i++)
    vec.push_back(vec[i]);

  return vec;
}

std::vector<double> str2vec(std::string magnitudes, double unidad){
  std::istringstream data(magnitudes);
  std::vector<double> data_array;
  double num;
  while(!data.eof()){
    data >> num;
    data_array.push_back(num*unidad);
  }
  return data_array;
}

unsigned int stepEstimator(unsigned int Niter, unsigned int tau, unsigned int min_periods){
  unsigned int periods = Niter/tau;
  if (periods > min_periods)
    return periods*tau;
  else
    return min_periods*tau;
}

std::vector<double> cosarray(unsigned int length, unsigned int tau, double amplitude, double phase){
  std::vector<double> wave;
  wave.resize(length);
  unsigned int wavetop = (tau>=length)? length : tau;
    for(unsigned int i=0; i<wavetop; i++)
      wave[i]=amplitude*cos(_2pi*i/tau-phase);
  
  unsigned int periods = length/tau;
  for(unsigned int i=1; i<periods;i++){
    for(unsigned int j = 0; j<tau ;j++)
      wave[i*tau+j]=wave[j];
  }
  return wave;
}

double simpson_int(const double f_array[], const std::vector<double>& weight){
  unsigned int length=weight.size();
  
  double Integral = f_array[0]*weight[0];
  
  for(unsigned int i=1; i<length-1; i+=2)
    Integral+=4*f_array[i]*weight[i];
  
  for(unsigned int i=2; i<length; i+=2)
    Integral+=2*f_array[i]*weight[i];
  
  length--;
  Integral+=f_array[length]*weight[length];
  
  return Integral/3;
}

//Calcular la desviación estandar de una matriz
double stan_dev(double ** M, unsigned int rows, unsigned int cols){
  unsigned int celdas;
  celdas = rows * cols;
  double * Aij = new double [celdas];
  for(unsigned int i = 0 ; i<rows; i++){
    for(unsigned int j = 0; j<cols; j++)
      Aij[i*cols + j] = M[i][j];
  }
  double sd = gsl_stats_sd (Aij, 1, celdas);
  delete[] Aij;
  return sd;
}

bool needSimulation(std::string id_proc, unsigned int size)
{
  struct stat file;
  id_proc = "log_pol_"+id_proc+".dat";
  
  if (stat(id_proc.c_str(), &file) == -1)
    return true;
  
  if (file.st_size != size){
    std::remove(id_proc.c_str());
    return true;
  }
  
  return false;
}

void Sistema::flip_sigma(unsigned int idsigma){sigma[idsigma] *= -1;}
unsigned int Sistema::return_PNR(){return PNR;}
int Sistema::ret_sig(unsigned int i){return sigma[i];}
int8_t* Sistema::ret_sigarr(){ return sigma;}
