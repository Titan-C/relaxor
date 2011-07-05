#include <cmath>
#include "sistema.h"
#include "impresor.h"
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_randist.h>

/*Constructor:
Dimensiona y encera a los vectores del sistema. Llena sus datos iniciales */
Sistema::Sistema(unsigned int lado,
		 unsigned int Niter,
		 gsl_rng* rng,
		 double Delta_J,
		 unsigned int dim,
		 bool polarizar){
  dimension = dim;
  L = lado;
  // Dimensionado de arreglos caraterísticos del sistema
  sigma.resize(pow(lado,dimension));  
  sum_sigma_time.resize(sigma.size());
  sum_sigma_conf.resize(Niter);
  mu_E.resize(sigma.size());

  G.resize(sigma.size());
  J.resize(sigma.size());
  unsigned int vecinos = 2*dimension;
  for(unsigned int i = 0; i < J.size(); i++){
    G[i].resize(vecinos);
    J[i].resize(vecinos);
  }

  // Inicializa al sistema, llenado de datos
  //Crear arreglos auxiliares
  DeltaJ = init(rng, Delta_J, polarizar);  
}

/*Destructor:
libera la memoria asignada a los vectores del sistema*/
Sistema::~Sistema(){
  for(unsigned int i = 0; i < J.size(); i++){
    G[i].clear();
    J[i].clear();
  }
  J.clear();
  G.clear();
  mu_E.clear();
  sum_sigma_conf.clear();
  sum_sigma_time.clear();
  sigma.clear();
}

int Sistema::init_pol(gsl_rng* rng, bool polarizar){
  if (polarizar)
    sigma.assign(sigma.size(),1);
  else{
    for(unsigned int i=0; i<sigma.size(); i++)
      sigma[i] = (gsl_rng_uniform(rng)-0.5 > 0)? 1:-1;
  }
  //Calcular la polarización inicial
  int ini_pol=0;
  for(unsigned int i=0; i<sigma.size(); i++)
    ini_pol+=sigma[i];

  return ini_pol;
}
void Sistema::mul(gsl_rng* rng, bool print)
{
  for(unsigned int i=0; i<sigma.size(); i++)
    mu_E[i]  = gsl_rng_uniform(rng);
  if (print)
    array_print(mu_E, "mu_H.dat", false);
}

void Sistema::topografia(bool print){
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
  if (print){
    array_print(R, "posiciones.dat");
    array_print(G, "grafo_vecinos.dat");
  }
  //liberar mem
  for(unsigned int i=0; i<R.size();i++)
    R[i].clear();
  R.clear();
}

double Sistema::Jex(gsl_rng* rng, double Delta_J, bool print){
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
  double _std=stan_dev(Jinter);

  // Elabora el arreglo de interacción de primeros vecinos
  for(unsigned int i=0; i<J.size(); i++){
    for(unsigned int j=0; j<J[i].size(); j++)
      J[i][j] = Jinter[i][G[i][j]];
  }
  if (print){
    array_print(Jinter, "Matriz_Intercambio.dat");
    array_print(J, "J_vecinos.dat");
  }
  //liberar mem
  for(unsigned int i=0; i<Jinter.size();i++)
    Jinter[i].clear();
  Jinter.clear();

  return _std;
}

double Sistema::init(gsl_rng* rng, double Delta_J, bool polarizar){
  std::cout<<"Polarización inicial="<<init_pol(rng, polarizar)<<std::endl; //Estado inicial del Sistema y su topología
  mul(rng,false);//Momento dipolar eléctrico en eje principal
  topografia(false);

  // Genera las energías de intercambio de las PNR
  std::cout<<"Des stan Total= "<<Jex(rng,Delta_J,false)<<std::endl;
  double stan_dev_val = stan_dev(J);
  std::cout<<"Des stan Jveci= "<<stan_dev_val<<std::endl;

  return stan_dev_val;
}

double Sistema::total_E(double E){
  double H = 0;
  for(unsigned int i = 0; i < G.size(); i++){
    for(unsigned int j = 0; j < G[i].size(); j++)
      H -= J[i][j]*sigma[i]*sigma[G[i][j]];
    H -= E*mu_E[i]*sigma[i];
  }
  return H;
}

double Sistema::delta_E(unsigned int idflip, double E){
  double dH = 0;
  for(unsigned int i = 0; i<G[idflip].size(); i++)
    dH += J[idflip][i]*sigma[idflip]*sigma[G[idflip][i]];
  dH +=E*mu_E[idflip]*sigma[idflip];
  return 2*dH;
}

void Sistema::flip(unsigned int idflip, double T, double E, gsl_rng* rng){
  double dH = delta_E(idflip, E);
  if ( dH < 0)
    sigma[idflip] *= -1;
  else if ( exp(-dH/T) >= gsl_rng_uniform(rng) )
    sigma[idflip] *= -1;
}

int Sistema::experimento(double T, double E, unsigned int Niter,
			 bool grabar, gsl_rng* rng, std::string id_proc){
  for(unsigned int i = 0; i< Niter; i++){
    //out(total_E(E), "energy_log.dat");
    for(unsigned int idflip = 0; idflip < sigma.size(); idflip++)
      flip(idflip, T, E, rng);
    
    if (grabar) {
      for(unsigned int j = 0; j< sigma.size(); j++){
	sum_sigma_time[j]+=sigma[j];
	sum_sigma_conf[i]+=sigma[j];
      }
    }
  }
  
  if (grabar) {
    array_print(sum_sigma_time, "sum_sigma_time_"+id_proc+".dat");
    array_print(sum_sigma_conf, "sum_sigma_conf_"+id_proc+".dat");
    reset_sum_sigma();
  }
  return 1;
}

void Sistema::reset_sum_sigma(){
  sum_sigma_time.assign(sum_sigma_time.size(), 0);
  sum_sigma_conf.assign(sum_sigma_conf.size(), 0);
}

double stan_dev(const std::vector< std::vector<double> >& M){
  //Calcular la desviación standar del las energías de intercambio.
  unsigned int celdas, columnas;
  columnas = M[1].size();
  celdas = M.size() * columnas;
  double * Jij;
  Jij = new double [celdas];
  for(unsigned int i = 0 ; i<M.size(); i++){
    for(unsigned int j = 0; j<columnas; j++)
      Jij[i*columnas + j] = M[i][j];
  }
  return gsl_stats_sd (Jij, 1, celdas);
  delete[] Jij;
}

std::vector<double> temp_array(double unidad, double T_top, double dT, bool heat){
  std::vector<double> temp( (int) ceil(1/dT*T_top) , 0);
  for(unsigned int i = 0; i<temp.size() ;i++){
    if (heat)
      temp[i] = (i+1)*dT*unidad;
    else
      temp[i] = (T_top - i*dT)*unidad;
  }
  
  return temp;
}
std::vector<double> field_array(double unidad, double H_top, double dH){
  std::vector<double> field (12,0);
  for(unsigned int i = 0 ; i< 5; i++)
    field[i] = (double) i*dH*unidad;
  for(unsigned int i = 0 ; i< 6; i++)
    field[i+5] = (double) (6+i*2)*dH*unidad;
  field[11] = 2*unidad;

  return field;
}

void Sus_proc(unsigned int numexps, unsigned int Niter, unsigned int L, double unidad,
	     const std::vector<double>& Temperatura, const std::vector<double>& campos, std::string id_proc){
  /*//vaciar datos de ejecuciones anteriores
  file_wipe("Congelamiento.dat");
  file_wipe("Susceptibilidad.dat");
  file_wipe("Xmax_Tmax.dat");*/

  //Procesar Dipolos congelados y Correspondiente Susceptibilidad
  //declarar variables
  std::vector< std::vector<double> > sigmas_time, S_frozen, Susceptibilidad;
  std::vector<double> slows(4,1);
  slows[1]=0.9;
  slows[2]=0.8;
  slows[3]=0.6;

  double Xmax, Tmax;
  unsigned int L3 = L*L*L, temp_size = Temperatura.size(), field_size = campos.size(), index;
  std::string file;

  //iniciar o llenar variables
  import_data(sigmas_time, "sum_sigma_time_"+id_proc+".dat", field_size*numexps*temp_size , L3);
  S_frozen.resize(temp_size);
  Susceptibilidad.resize(temp_size);
  for(unsigned int i = 0; i < S_frozen.size(); i++){
    S_frozen[i].assign(slows.size(), 0);
    Susceptibilidad[i].assign(slows.size(), 0);
  }
  
  //Procesar para cada valor del campo
  for(unsigned int E = 0 ; E < field_size; E++){

    //Promediador de los dipolos congelados en los experimentos
    for(unsigned int n = 0; n < numexps; n++){
      for(unsigned int T = 0; T < temp_size; T++){
	for(unsigned int k = 0; k < L3; k++){
	  index = (E*numexps+n)*temp_size+T;
	  sigmas_time[index][k] = std::abs(sigmas_time[index][k]/Niter);
	  //Sumar dipolos congelados en ponderación
	  for(unsigned int s = 0; s < slows.size(); s++){
	    if (sigmas_time[index][k] >= slows[s] )
	      S_frozen[T][s] += (double) 1/L3/numexps;
	  }//termina proporción congelada por celda
	}//termina todas las celdas
      }//termina rango temperatura del experimento
    }//termina ejecución del experimento
    array_print(S_frozen, "Congelamiento_"+id_proc+".dat", true);

    //Encontrar la susceptibilidad de los experimentos promediados
    for(unsigned int T = 0; T < temp_size; T++){
      for(unsigned int s = 0; s < slows.size(); s++)
	Susceptibilidad[T][s] = unidad*(1 - S_frozen[T][s])/Temperatura[T];
    }
    array_print(Susceptibilidad, "Susceptibilidad_"+id_proc+".dat", true);

    //encontrar Xmax y Tmax
    Tmax = 0;
    Xmax = 0;
    for(unsigned int T = 0; T< temp_size; T++){
      if (Susceptibilidad[T][1] >= Xmax){
	Xmax = Susceptibilidad[T][1];
	Tmax = Temperatura[T];
      }
    }
    out(campos[E]/unidad,"Xmax_Tmax"+id_proc+".dat", true, false);
    out(Xmax,"Xmax_Tmax"+id_proc+".dat", true, false);
    out(Tmax/unidad,"Xmax_Tmax"+id_proc+".dat");

    //Limpiar Matriz de dipolos Congelados
    for(unsigned int f = 0; f < S_frozen.size(); f++){
      S_frozen[f].assign(slows.size(), 0);
    }
  }//termina con los experimentos a un campo dado

  //liberar mem
  for(unsigned int i=0;i<sigmas_time.size();i++)
    sigmas_time[i].clear();
  sigmas_time.clear();
  
  for(unsigned int i= 0;i<S_frozen.size();i++){
    S_frozen[i].clear();
    Susceptibilidad[i].clear();
  }
  S_frozen.clear();
  Susceptibilidad.clear();
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
void plot_data_sus(double unidad, const std::vector<double>& Temperatura, const std::vector<double>& campos, std::string id_proc){
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

// Aplica las condiciones de borde toroidales
void condborde ( std::vector <double>& R, int L){
  double bL=L/2.0;
  for(unsigned int i=0; i<R.size();i++){
    if (R[i]>bL)	R[i]-=L;
    else if (R[i]<-bL)	R[i]+=L;
  }
}
//Producto interno
double dot(const std::vector< double >& a, const std::vector< double >& b){
  double A=0;
  for(unsigned int i=0;i<a.size();i++)
    A+=a[i]*b[i];
  return A;
}