#include <cmath>
#include <ctime>
#include <fstream>
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

  // Inicializa al sistema, llenado de datos
  //Crear arreglos auxiliares
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
  //Generar configuración espacial de PNR
  set_space_config();
  
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

void Sistema::flip(unsigned int idflip, double T, double E, gsl_rng* rng){
  double dH = delta_E(idflip, E);
  if ( dH < 0)
    sigma[idflip] *= -1;
  else if ( exp(-dH/T) >= gsl_rng_uniform(rng) )
    sigma[idflip] *= -1;
}

int Sistema::experimento(double T, double E, unsigned int tau, unsigned int Niter,
			 bool grabar, gsl_rng* rng, std::string id_proc){
  std::vector<double> pol_mag;
  pol_mag.resize(Niter);
  
  std::vector<double> wave;
  wave.resize(tau);
  for(unsigned int i=0; i<tau; i++)
    wave[i]=E*std::cos( _2pi*i/tau );
    
  unsigned int periods = Niter/tau;
  unsigned int step = 0;
  for(unsigned int i = 0; i<periods; i++){
    for(unsigned int j = 0; j< tau; j++){
//       out(total_E(E), "energy_log.dat");
      for(unsigned int idflip = 0; idflip < sigma.size(); idflip++)
	flip(idflip, T, wave[j] , rng);
      
      if (grabar)
	pol_mag[step] = norm_pol();
      step++;
    }
  }
  if (grabar)
    array_print_bin(pol_mag,"polarizacion_"+id_proc+".dat");

  pol_mag.clear();
  wave.clear();
  
  return 1;
}

double stan_dev(const std::vector< std::vector<double> >& M){
  //Calcular la desviación standar del las energías de intercambio.
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

std::vector<double> temp_array(double unidad, double T_top, double dT, bool heat){
  std::vector<double> temp;
  temp.resize( (int) ceil(1/dT*T_top));
  for(unsigned int i = 0; i<temp.size() ;i++){
    if (heat)
      temp[i] = (i+1)*dT*unidad;
    else
      temp[i] = (T_top - i*dT)*unidad;
  }
  
  return temp;
}
std::vector<double> field_array(double unidad, double H_top, double dH){
  std::vector<double> field;
  field.resize(12);
  for(unsigned int i = 0 ; i< 5; i++)
    field[i] = (double) i*dH*unidad;
  for(unsigned int i = 0 ; i< 6; i++)
    field[i+5] = (double) (6+i*2)*dH*unidad;
  field[11] = 2*unidad;

  return field;
}

void calc_sus(unsigned int numexps, unsigned int tau, unsigned int Niter, unsigned int L, double unidad,
	     const std::vector<double>& Temperatura, const double campos, std::string id_proc){
  std::vector<double> pol_hist;
  pol_hist.resize(Niter);
  
  std::vector<double> cos_wave, sin_wave;
  cos_wave.resize(tau);
  sin_wave.resize(tau);
  for(unsigned int i=0; i<tau; i++){
    cos_wave[i]=std::cos( _2pi*i/tau );
    sin_wave[i]=std::sin( _2pi*i/tau );
  }
  std::string name ="polarizacion_"+id_proc+".dat";
  std::ifstream file (name.c_str());
  
  std::vector< std::vector<double> > X_mat;
  X_mat.resize(Temperatura.size());
  for(unsigned int T=0;T<Temperatura.size();T++){
    X_mat[T].resize(7);
    X_mat[T][0]=Temperatura[T]/unidad;
  }
  
  double * Freal = new double [numexps*Temperatura.size()];
  double * Fimag = new double [numexps*Temperatura.size()];
  unsigned int periods = Niter/tau;
    for(unsigned int n = 0; n < numexps; n++){    
      for(unsigned int T=0;T<Temperatura.size();T++){
	file.read((char * )&pol_hist[0],Niter*sizeof(double));
	
	unsigned int step = 0;
	double Int_cos=0, Int_sin=0;
	for(unsigned int i = 0; i<periods; i++){
	  for(unsigned int j = 0; j< tau; j++){
	    Int_cos+=pol_hist[step]*cos_wave[j];
	    Int_sin+=pol_hist[step]*sin_wave[j];
	    step++;
	  }
	}
	Freal[n*Temperatura.size()+T]=Int_cos;
	Fimag[n*Temperatura.size()+T]=Int_sin;
      }
    }
    
    for(unsigned int T=0;T<Temperatura.size();T++){
      X_mat[T][1]=X_mat[T][1]/Niter;
      X_mat[T][2]=X_mat[T][2]/Niter;
    }
  file.close();
  array_print(X_mat, "susceptibilidad_"+id_proc+".dat");
  /*//vaciar datos de ejecuciones anteriores
  file_wipe("Congelamiento.dat");
  file_wipe("Susceptibilidad.dat");
  file_wipe("Xmax_Tmax.dat");*/

//   //Procesar Dipolos congelados y Correspondiente Susceptibilidad
//   //declarar variables
//   std::vector< std::vector<double> > sigmas_time, S_frozen, Susceptibilidad;
// 
//   double Xmax, Tmax;
//   unsigned int L3 = L*L*L, temp_size = Temperatura.size(), field_size = campos.size(), index;
//   std::string file;
// 
//   //iniciar o llenar variables
//   import_data(sigmas_time, "sum_sigma_time_"+id_proc+".dat", field_size*numexps*temp_size , L3);
//   S_frozen.resize(temp_size);
//   Susceptibilidad.resize(temp_size);
//   for(unsigned int i = 0; i < S_frozen.size(); i++){
//     S_frozen[i].assign(slows.size(), 0);
//     Susceptibilidad[i].assign(slows.size(), 0);
//   }
//   
//   //Procesar para cada valor del campo
//   for(unsigned int E = 0 ; E < field_size; E++){
// 
//     //Promediador de los dipolos congelados en los experimentos
//     for(unsigned int n = 0; n < numexps; n++){
//       for(unsigned int T = 0; T < temp_size; T++){
// 	for(unsigned int k = 0; k < L3; k++){
// 	  index = (E*numexps+n)*temp_size+T;
// 	  sigmas_time[index][k] = std::abs(sigmas_time[index][k]/Niter);
// 	  //Sumar dipolos congelados en ponderación
// 	  for(unsigned int s = 0; s < slows.size(); s++){
// 	    if (sigmas_time[index][k] >= slows[s] )
// 	      S_frozen[T][s] += (double) 1/L3/numexps;
// 	  }//termina proporción congelada por celda
// 	}//termina todas las celdas
//       }//termina rango temperatura del experimento
//     }//termina ejecución del experimento
//     array_print(S_frozen, "Congelamiento_"+id_proc+".dat", true);
// 
//     //Encontrar la susceptibilidad de los experimentos promediados
//     for(unsigned int T = 0; T < temp_size; T++){
//       for(unsigned int s = 0; s < slows.size(); s++)
// 	Susceptibilidad[T][s] = unidad*(1 - S_frozen[T][s])/Temperatura[T];
//     }
//     array_print(Susceptibilidad, "Susceptibilidad_"+id_proc+".dat", true);
// 
//     //encontrar Xmax y Tmax
//     Tmax = 0;
//     Xmax = 0;
//     for(unsigned int T = 0; T< temp_size; T++){
//       if (Susceptibilidad[T][1] >= Xmax){
// 	Xmax = Susceptibilidad[T][1];
// 	Tmax = Temperatura[T];
//       }
//     }
//     out(campos[E]/unidad,"Xmax_Tmax"+id_proc+".dat", true, false);
//     out(Xmax,"Xmax_Tmax"+id_proc+".dat", true, false);
//     out(Tmax/unidad,"Xmax_Tmax"+id_proc+".dat");
// 
//     //Limpiar Matriz de dipolos Congelados
//     for(unsigned int f = 0; f < S_frozen.size(); f++){
//       S_frozen[f].assign(slows.size(), 0);
//     }
//   }//termina con los experimentos a un campo dado
// 
//   //liberar mem
//   for(unsigned int i=0;i<sigmas_time.size();i++)
//     sigmas_time[i].clear();
//   sigmas_time.clear();
//   
//   for(unsigned int i= 0;i<S_frozen.size();i++){
//     S_frozen[i].clear();
//     Susceptibilidad[i].clear();
//   }
//   S_frozen.clear();
//   Susceptibilidad.clear();
}
void Pol_proc(unsigned int numexps, unsigned int Niter, unsigned int L, double unidad,
	     const std::vector<double>& Temperatura, const std::vector<double>& campos, std::string id_proc){
  //Procesar Polarización del sistema
  //iniciar variables
  std::vector< std::vector<double> > sigmas_conf, Polarizacion;
  unsigned int L3 = L*L*L, temp_size = Temperatura.size(), field_size = campos.size(), index;

  import_data(sigmas_conf, "sum_sigma_conf_"+id_proc+".dat", field_size*numexps*temp_size , Niter);
  Polarizacion.resize(temp_size);
  for(unsigned int T = 0; T< temp_size; T++)
    Polarizacion[T].assign(field_size+1,0);
  for(unsigned int E = 0; E < field_size; E++){
    //Promediar la polarización
    for(unsigned int n = 0; n< numexps; n++){
      for(unsigned int T = 0; T< temp_size;T++){
	Polarizacion[T][0] = Temperatura[T]/unidad;
	for(unsigned int k = 0; k < Niter; k++){
	  index = (E*numexps+n)*temp_size+T;
	  Polarizacion[T][E+1] += (double) sigmas_conf[index][k]/L3/Niter/numexps;
	}}}
  }
  array_print(Polarizacion, "polarizacion_" + id_proc +".dat");
  
  //liberar mem
  for(unsigned int i=0;i<sigmas_conf.size();i++)
    sigmas_conf[i].clear();
  sigmas_conf.clear();
  for(unsigned int i=0;i<Polarizacion.size();i++)
    Polarizacion[i].clear();
  Polarizacion.clear();
}
void plot_data_sus(double unidad, const std::vector< double >& Temperatura, const std::vector< double >& campos, std::string id_proc){
  std::vector< std::vector<double> > Frozen, Susceptibilidad, Polarizacion;
  unsigned int temp_size = Temperatura.size(), field_size = campos.size();
  import_data(Frozen, "Congelamiento_"+id_proc+".dat", field_size*temp_size, 4);
  import_data(Susceptibilidad, "Susceptibilidad_"+id_proc+".dat", field_size*temp_size, 4);

  std::vector< std::vector<double> > plot_array;
  plot_array.resize(temp_size);
  //Congelados 90%
  for(unsigned int i = 0; i < temp_size; i++){
    plot_array[i].assign(field_size+1,0);
    plot_array[i][0] = Temperatura[i]/unidad;
    for(unsigned int j = 0; j < field_size; j++)
      plot_array[i][j+1] = Frozen[j*temp_size + i][1];
  }
  array_print(plot_array, "Dip_cong_90_" + id_proc + ".dat");
  
  //Suceptibilidad 90%
  for(unsigned int i=0; i<temp_size;i++){
    for(unsigned int j = 0; j < field_size; j++)
      plot_array[i][j+1] = Susceptibilidad[j*temp_size + i][1];
  }
  array_print(plot_array, "Sus_90_" + id_proc + ".dat");
  
  //Campo nulo congelados y Susceptibilidad
  for(unsigned int i = 0; i < temp_size; i++){
    plot_array[i].assign(5,0);
    plot_array[i][0] = Temperatura[i]/unidad;
    for(unsigned int j = 0; j < 4; j++)
      plot_array[i][j+1] = Frozen[i][j];
  }
  array_print(plot_array, "Dip_cong_Efix" + id_proc + ".dat");
  
  for(unsigned int i = 0; i < temp_size; i++){
    plot_array[i].assign(6,0);
    plot_array[i][0] = Temperatura[i]/unidad;
    plot_array[i][1] = unidad/Temperatura[i];
    for(unsigned int j = 0; j < 4; j++)
      plot_array[i][j+2] = Susceptibilidad[i][j];
  }
  array_print(plot_array, "Sus_Efix_" + id_proc + ".dat");
}