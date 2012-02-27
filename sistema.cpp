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
		 unsigned int dim,
		 bool polarizar){
  dimension = dim;
  L = lado;
  PNR = pow(lado,dimension);
  // Dimensionado de arreglos caraterísticos del sistema
  sigma.resize(PNR);
  mu_E.resize(PNR);

  J.resize(PNR);
  G.resize(PNR);
  vecinos = 2*dimension;
  for(unsigned int i = 0; i < J.size(); i++){
    J[i].resize(vecinos);
    G[i].resize(vecinos);
  }
  //Generar configuración espacial de PNR
  set_space_config();
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

double Sistema::Jex(gsl_rng* rng){
  std::vector< std::vector<double> > Jinter;
  Jinter.resize(PNR);
  /*Genera las matriz triangular superior de las
   * energías de intecambio segun una distribución
   * Gaussiana*/
  for(unsigned int i = 0; i<PNR; i++){
    Jinter[i].resize(PNR);
    for(unsigned int j = i+1; j<PNR; j++)
      Jinter[i][j] = gsl_ran_gaussian(rng,DeltaJ)+rho;
  }
  //Completa la parte inferior de la matriz de intercambio
  for(unsigned int i = 0; i<PNR; i++){
    for(unsigned int j = i+1; j<PNR; j++)
      Jinter[j][i] = Jinter[i][j];
  }
  // Elabora el arreglo de interacción de primeros vecinos
  for(unsigned int i=0; i<PNR; i++){
    for(unsigned int j=0; j<6; j++)
      J[i][j] = Jinter[i][G[i][j]];
  }
  //liberar mem
  for(unsigned int i=0; i<Jinter.size();i++)
    Jinter[i].clear();
  Jinter.clear();

  return stan_dev(J);
}

double Sistema::set_pol(gsl_rng* rng, bool polarizar){
  if (polarizar)
    sigma.assign(PNR,1);
  else{
    for(unsigned int i=0; i<PNR; i++)
      sigma[i] = (gsl_rng_uniform(rng)-0.5 > 0)? 1:-1;
  }

  return norm_pol();
}

double Sistema::set_mu(gsl_rng* rng, bool polarizar){
  for(unsigned int i=0; i<PNR; i++)
    mu_E[i]  = gsl_rng_uniform(rng);

//   //Tomar la suma de interacción de cada región con sus vecinos
//   for(unsigned int i=0; i<J.size(); i++){
//     for(unsigned int j=0; j<J[i].size(); j++)
//       mu_E[i]+=J[i][j];
//   }
//   //Encontrar el Maximo
//   double max=0;
//   for(unsigned int i=0;i<mu_E.size();i++)
//     if ( std::abs(mu_E[i]) > max) max=abs(mu_E[i]);
//   //normar
//   for(unsigned int i=0; i<mu_E.size();i++)
//     mu_E[i]/=max;
//   
//   array_print(mu_E,"polarizacion.dat");
  return set_pol(rng, polarizar);
}

void Sistema::init(gsl_rng* rng, double DJ, double p, bool polarizar, bool write){
  //Cambia las propidades del sistema
  DeltaJ = DJ;
  rho = p;
  // Genera las energías de intercambio de las PNR
  double Delta_J=Jex(rng);
  //Momentos dipolares y polarización
  double pol=set_mu(rng, polarizar);
  
  if (write){
    std::cout<<"Desvición Estandar Total= "<<Delta_J<<"\n";
    std::cout<<"Polarización inicial="<<pol<<"\n";
  }
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
  int s = sigma[idflip];
  for(unsigned int i = 0; i<vecinos; i++)
    dHamil += J[idflip][i]*s*sigma[G[idflip][i]];
  dHamil *=4;
  dHamil += 2*E*mu_E[idflip]*s;
  return dHamil;
}

double Sistema::norm_pol(){
  double P=0;
  for(unsigned int i=0; i<PNR; i++)
    P += mu_E[i]*sigma[i];

  return (double) P / PNR;
}

void Sistema::experimento(double T, double E, unsigned int tau, unsigned int Niter,
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
      for(unsigned int idflip = 0; idflip < PNR; idflip++){
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
}

void Gen_exp(std::vector< double >& Temps, std::vector< double >& Fields,
		      std::vector< double > tau, unsigned int numexps, double DJ,
		      double p, unsigned int L, unsigned int Equi_iter, unsigned int Exp_iter,
		      std::string Exp_ID, gsl_rng* rng)
{
  clock_t cl_start = clock();
  Sistema relaxor(L, rng);
  if (Exp_ID == "cool" || Exp_ID == "heat"){
    for(unsigned int t=0; t< tau.size() ; t++){
      for(unsigned int E=0; E<Fields.size(); E++){
	std::ostringstream id_proc;
	id_proc<<Exp_ID<<"_J"<<DJ<<"_p"<<p<<"_E"<<Fields[E]/DJ<<"_t"<<tau[t];
	for(unsigned int n=0; n<numexps; n++){
	  relaxor.init(rng,DJ,p,false);
	  for(unsigned int T=0; T<Temps.size(); T++){
	    relaxor.experimento(Temps[T],Fields[E],tau[t], Equi_iter,false,rng, id_proc.str());
	    relaxor.experimento(Temps[T],Fields[E],tau[t], Exp_iter,true,rng, id_proc.str());
	  }}}}}
  else {
    for(unsigned int t=0; t< tau.size() ; t++){
      for(unsigned int T=0; T<Temps.size(); T++){
	std::ostringstream id_proc;
	id_proc<<Exp_ID<<"_J"<<DJ<<"_p"<<p<<"_T"<<Temps[T]/DJ<<"_t"<<tau[t];
	for(unsigned int n=0; n<numexps; n++){
	  relaxor.init(rng,DJ,p,false);
	  for(unsigned int E=0; E<Fields.size(); E++){
	    relaxor.experimento(Temps[T],Fields[E],tau[t], Equi_iter,false,rng, id_proc.str());
	    relaxor.experimento(Temps[T],Fields[E],tau[t], Exp_iter,true,rng, id_proc.str());
	  }}}}}

  proces_data(Temps,Fields,tau,numexps, DJ,p,Exp_iter,Exp_ID);
  std::cout<<Exp_ID<<":"<<clock()-cl_start<<"\n";
}

void proces_data(std::vector< double >& Temps, std::vector< double >& Fields,
		 std::vector< double > tau, unsigned int numexps, double DJ,
		 double p, unsigned int Niter, std::string Exp_ID){
  //procesar los datos
  for (unsigned int t=0;t<tau.size();t++){
    for (unsigned int E=0;E<Fields.size();E++){
      for (unsigned int T=0; T<Temps.size();T++){
	std::ostringstream id_proc;
	std::vector<double> pol_stats, pol_int_avg;
	if (Exp_ID == "cool" || Exp_ID == "heat"){
	  id_proc<<Exp_ID<<"_J"<<DJ<<"_p"<<p<<"_E"<<Fields[E]/DJ<<"_t"<<tau[t];
	  pp_data(pol_stats,pol_int_avg,Temps.size(),numexps,tau[t],Niter,id_proc.str());
	  std::vector<double> intfield (1,Fields[E]);
	  eval_pol(pol_stats,numexps,DJ,Temps,id_proc.str(),true);
	  calc_sus(pol_int_avg,numexps,DJ,Temps,intfield,id_proc.str());
	  intfield.clear();
	} else {
	  id_proc<<Exp_ID<<"_J"<<DJ<<"_p"<<p<<"_T"<<Temps[T]/DJ<<"_t"<<tau[t];
	  pp_data(pol_stats,pol_int_avg,Fields.size(),numexps,tau[t],Niter,id_proc.str());
	  eval_pol(pol_stats,numexps,DJ,Fields,id_proc.str(),(Exp_ID=="hist_loop") ? false : true);
	  calc_sus(pol_int_avg,numexps,DJ,Fields,Fields,id_proc.str());
	}
	pol_int_avg.clear();
	pol_stats.clear();
  }}}
  
  //Graficar
  plot_sus(Exp_ID,DJ,p,Temps,Fields,tau);
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
  cos_wave = waves(Niter,tau,1.0,true);
  sin_wave = waves(Niter,tau,1.0,false);
  
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
void eval_pol(const std::vector<double>& pol_stats, unsigned int numexps, double unidad, const std::vector<double>& x_array, std::string id_proc, bool absolut) {

  //Polarización, o polarización absoluta y desviación estandar
  double * data_array = new double [numexps];
  double data_length = x_array.size();
  std::vector< std::vector<double> > pol_final;
  pol_final.resize(data_length);
  for(unsigned int x=0;x< data_length; x++){
    pol_final[x].resize(3);
    pol_final[x][0]=x_array[x]/unidad;
    
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
  for(unsigned int i=0;i<data_length;i++){
    pol_final[i].clear();
  }
  pol_final.clear();
}

void calc_sus(const std::vector<double>& pol_int_avg, unsigned int numexps, double unidad,
	      const std::vector<double>& x_array, const std::vector<double>& campo, std::string id_proc){

  /*Calcular susceptibildad más error*/
  double * data_arrayr = new double [numexps];
  double * data_arrayi = new double [numexps];
  bool fieldvec = (campo.size()>1) ? true : false;
  double data_length = x_array.size();
  std::vector< std::vector<double> > X_mat;
  X_mat.resize(data_length);
  for(unsigned int x=0;x<data_length;x++){
    X_mat[x].resize(5);
    X_mat[x][0]=x_array[x]/unidad;
    
    for(unsigned int n=0;n<numexps;n++){
      unsigned int ind = 2*(n*data_length+x);
      unsigned int field_ind=(fieldvec) ? x : 0;
      data_arrayr[n]=pol_int_avg[ind]/campo[field_ind];
      data_arrayi[n]=pol_int_avg[ind+1]/campo[field_ind];
    }
    X_mat[x][1]=gsl_stats_mean(data_arrayr,1,numexps);
    X_mat[x][2]=gsl_stats_sd_m(data_arrayr,1,numexps,X_mat[x][1]);
    X_mat[x][3]=gsl_stats_mean(data_arrayi,1,numexps);
    X_mat[x][4]=gsl_stats_sd_m(data_arrayi,1,numexps,X_mat[x][3]);
  }
  array_print(X_mat, "sus_"+id_proc+".dat");

  //liberar memoria
  for(unsigned int i=0; i<X_mat.size();i++)
    X_mat[i].clear();
  X_mat.clear();  
  delete[] data_arrayr;
  delete[] data_arrayi;
}

std::vector<double> step2vec(double unidad, double v_start, double v_end, double dv, std::vector<double> last){
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

std::vector<double> loop2vec(double unidad, double max, int divs){
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

std::vector<double> waves(unsigned int length, unsigned int tau, double amplitude, bool cossin){
  std::vector<double> wave;
  wave.resize(length);
  if (cossin) {
    for(unsigned int i=0; i<tau; i++)
      wave[i]=amplitude*cos(_2pi*i/tau);
  }
  else {
    for(unsigned int i=0; i<tau; i++)
      wave[i]=amplitude*sin(_2pi*i/tau);
  }
  
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
double stan_dev(const std::vector< std::vector<double> >& M){
  unsigned int celdas, columnas;
  columnas = M[1].size();
  celdas = M.size() * columnas;
  double * Aij = new double [celdas];
  for(unsigned int i = 0 ; i<M.size(); i++){
    for(unsigned int j = 0; j<columnas; j++)
      Aij[i*columnas + j] = M[i][j];
  }
  double sd = gsl_stats_sd (Aij, 1, celdas);
  delete[] Aij;
  return sd;
}