#include "sistema.h"
#include <cmath>
#include <ctime>
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
		 bool polarizar)
{
  dimension = dim;
  L = lado;
  P=1;
  // Dimensionado de arreglos caraterísticos del sistema
  sigma.resize(pow(lado,dimension));  
  sum_sigma_time.resize(sigma.size());
  sum_sigma_conf.resize(Niter);
  mu_H.resize(sigma.size());

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
Sistema::~Sistema()
{
  J.clear();
  G.clear();
  mu_H.clear();
  sum_sigma_conf.clear();
  sum_sigma_time.clear();
  sigma.clear();
}
// Genero los datos del sistema
double Sistema::init(gsl_rng* rng, double Delta_J, bool polarizar){
  //Creo variable auxiliares
  unsigned int ind_xy, L2=L*L;
  std::vector< std::vector<unsigned int> > R;
  R.resize(sigma.size());

  //Estado inicial del Sistema y su topología
  for(unsigned int i=0; i<sigma.size(); i++){

    // Polarización de sistema
    if (polarizar)
      sigma[i] = 1;
    else
      sigma[i] = (gsl_rng_uniform(rng)-0.5 > 0)? 1:-1;

    //Momento dipolar eléctrico en eje principal
    mu_H[i]  = gsl_rng_uniform(rng);

    // Coeficientes vector posición i-ésima PNR
    ind_xy = i % L2;
    R[i].resize(dimension);
    R[i][0] = ind_xy % L;
    R[i][1] = ind_xy / L;
    R[i][2] = i / L2;

    /*Encontrar índices de los primeros vecinos.
      solo existen 6: arriba y abajo(+z, -z), derecha e izquierda(+y, -y), adelante y atraz(+x, -x).
      También debo aplicar las condiciones de borde en este caso */
    G[i][0] = (R[i][2] == L-1 )	?i - (L-1)*L2	:i + L2;//arriba
    G[i][1] = (R[i][2] == 0 )	?i + (L-1)*L2	:i - L2;//abajo
    G[i][2] = (R[i][1] == L-1 )	?i - (L-1)*L	:i + L;//derecha
    G[i][3] = (R[i][1] == 0)	?i + (L-1)*L	:i - L;//izquierda
    G[i][4] = (R[i][0] == L-1 )	?i - L+1	:i + 1;//adelante
    G[i][5] = (R[i][0] == 0 )	?i + L-1	:i - 1;//atraz    
  }
  array_print(mu_H, "mu_H.dat", false);
  array_print(R, "posiciones.dat");
  R.clear();
  array_print(G, "grafo_vecinos.dat");

  // Calcular las energías de intercambio de las PNR
  std::vector< std::vector<double> > Jinter;
  Jinter.resize(sigma.size());
  for(unsigned int i = 0; i<Jinter.size(); i++){
    Jinter[i].resize(sigma.size());
    for(unsigned int j = i+1; j<Jinter[i].size(); j++){
      // Calcula la Energía de intercambio entre 2 PNR
      Jinter[i][j] = gsl_ran_gaussian(rng,Delta_J);
    }
  }

  //Completa la parte inferior de la matriz de intercambio
  for(unsigned int i = 0; i<Jinter.size(); i++){
    for(unsigned int j = i+1; j<Jinter.size(); j++)
      Jinter[j][i] = Jinter[i][j];
  }
  array_print(Jinter, "Matriz_Intercambio.dat");

  // Elabora el arreglo de interacción de primeros vecinos
  for(unsigned int i=0; i<J.size(); i++){
    for(unsigned int j=0; j<J[i].size(); j++)
      J[i][j] = Jinter[i][G[i][j]];
  }
  array_print(J, "J_vecinos.dat");

  double stan_dev_val = stan_dev(J);

  std::cout<<"Des stan Jveci= "<<stan_dev_val<<std::endl;
  std::cout<<"Des stan Total= "<<stan_dev(Jinter)<<std::endl;
  Jinter.clear();

  return stan_dev_val;
}
//Calcular la desviación standar del las energías de intercambio.
double stan_dev(const std::vector< std::vector<double> >& M)
{
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

//Calcula la energía total del sistema
double Sistema::total_E(double E)
{
  double H = 0;
  for(unsigned int i = 0; i < G.size(); i++){
    for(unsigned int j = 0; j < G[i].size(); j++)
      H -= J[i][j]*sigma[i]*sigma[G[i][j]];
    H -= E*mu_H[i]*sigma[i];
  }
  return H;
}
//Calcula la variación de energía del sistema debído a un cambio del spin dipolar
double Sistema::delta_E(unsigned int idflip, double E)
{
  double dH = 0;
  for(unsigned int i = 0; i<G[idflip].size(); i++)
    dH += J[idflip][i]*sigma[idflip]*sigma[G[idflip][i]];
  dH +=E*mu_H[idflip]*sigma[idflip];
  return 2*dH;
}
//realiza el cambio del spin dipolar en una ubicación dada
void Sistema::flip(unsigned int idflip, double T, double E, gsl_rng* rng)
{
  double dH = delta_E(idflip, E);
  if ( dH < 0) sigma[idflip] *= -1;
  else if ( exp(-dH/T) >= gsl_rng_uniform(rng) ) sigma[idflip] *= -1;
}

int Sistema::experimento(double T, double E, unsigned int Niter,
			 bool grabar, gsl_rng* rng)
{
  for(unsigned int i = 0; i< Niter; i++){
    out(total_E(E), "energy_log.dat");
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
    array_print(sum_sigma_time, "sum_sigma_time.dat");
    array_print(sum_sigma_conf, "sum_sigma_conf.dat");
    reset_sum_sigma();
  }  
  return 1;
}

void Sistema::reset_sum_sigma()
{
  sum_sigma_time.assign(sum_sigma_time.size(), 0);
  sum_sigma_conf.assign(sum_sigma_conf.size(), 0);
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
void temp_array(std::vector< double >& Temperatura, double unidad, double T, double dT)
{
  Temperatura.resize(int (T/dT));
  for(unsigned int i = 0; i<Temperatura.size() ;i++)
    Temperatura[i] = (T - i*dT)*unidad;
}
void field_array(std::vector< double >& campo, double unidad)
{
  campo.resize(23);
  for(unsigned int i = 0 ; i< 21; i++)
    campo[i] = (double) i/10*unidad;
  campo[21] = 5*unidad;
  campo[22] = 10*unidad;
}


void procesar(unsigned int numexps, unsigned int Niter, unsigned int L, double unidad, const std::vector<double>& Temperatura, const std::vector<double>& campos)
{
  std::vector< std::vector<double> > sigmas_time, S_frozen, Susceptibilidad;
  double Xmax, Tmax;
  unsigned int lentos = 4, L3 = L*L*L, temp_size = Temperatura.size(), field_size = campos.size(), index;

  import_data(sigmas_time, "sum_sigma_time.dat", numexps*temp_size*field_size , L3);
  file_wipe("Congelamiento.dat");
  file_wipe("Susceptibilidad.dat");
  file_wipe("Xmax_Tmax.dat");

  S_frozen.resize(temp_size);
  Susceptibilidad.resize(temp_size);
  for(unsigned int i = 0; i < S_frozen.size(); i++){
    S_frozen[i].assign(lentos+1, 0);
    S_frozen[i][lentos] = Temperatura[i]/unidad;
    Susceptibilidad[i].assign(S_frozen[i].begin(), S_frozen[i].end());
  }
  for(unsigned int E = 0 ; E < field_size; E++){
    //Promediador de los dipolos congelados en los experimentos
    for(unsigned int i = 0; i < numexps; i++){
      for(unsigned int j = 0; j < temp_size; j++){

	for(unsigned int k = 0; k < L3; k++){
	  index = (E*numexps+i)*temp_size+j;
	  sigmas_time[index][k] = std::abs(sigmas_time[index][k]/Niter);
	  //Sumar dipolos congelados en ponderación
	  for(unsigned int l = 0; l < lentos; l++){
	    if (sigmas_time[index][k] >= (1-0.1*l) )
	      S_frozen[j][l] += (double) 1/L3/numexps;
	  }
	}
      }
    }
    array_print(S_frozen, "Congelamiento.dat", true);

    //Encontrar la susceptibilidad
    for(unsigned int j = 0; j < temp_size; j++){
      for(unsigned int l = 0; l < lentos; l++)
	Susceptibilidad[j][l] = unidad*(1 - S_frozen[j][l])/Temperatura[j];
    }
    array_print(Susceptibilidad, "Susceptibilidad.dat", true);
    //encontrar sus max t tmax
    Tmax = 0;
    Xmax = 0;
    for(unsigned int j = 0; j< temp_size; j++){
      if (Susceptibilidad[j][1] >= Xmax){
	Xmax = Susceptibilidad[j][1];
	Tmax = Temperatura[j];
      }
    }
    out(campos[E]/unidad,"Xmax_Tmax.dat", true, false);
    out(Xmax,"Xmax_Tmax.dat", true, false);
    out(Tmax,"Xmax_Tmax.dat");

    //Limpiar Matriz de dipolos Congelados
    for(unsigned int f = 0; f < S_frozen.size(); f++){
      S_frozen[f].assign(lentos+1, 0);
      S_frozen[f][lentos] = Temperatura[f]/unidad;
    }
  }

  sigmas_time.clear();
  S_frozen.clear();
  Susceptibilidad.clear();


  std::vector< std::vector<double> > sigmas_conf, Polarizacion;
  import_data(sigmas_conf, "sum_sigma_conf.dat", numexps*temp_size*field_size , Niter);
  Polarizacion.resize(temp_size);
  for(unsigned int j = 0; j< temp_size; j++)
    Polarizacion[j].assign(field_size,0);

  for(unsigned int E = 0; E < field_size; E++){
    //Promediar la polarización
    for(unsigned int i = 0; i< numexps; i++){
      for(unsigned int j = 0; j< temp_size;j++){
	for(unsigned int k = 0; k < Niter; k++){
	  index = (E*numexps+i)*temp_size+j;
	  Polarizacion[j][E] += (double) sigmas_conf[index][k]/Niter/L3/numexps;
	}}}
  }
  array_print(Polarizacion, "polarizacion.dat");
  sigmas_conf.clear();
  Polarizacion.clear();
}

void graficos(double unidad, const std::vector<double>& Temperatura, const std::vector<double>& campos)
{
  std::vector< std::vector<double> > Frozen, Susceptibilidad, Polarizacion;
  unsigned int temp_size = Temperatura.size(), field_size = campos.size();
  import_data(Frozen, "Congelamiento.dat", temp_size*field_size, 5);
  import_data(Susceptibilidad, "Susceptibilidad.dat", temp_size*field_size, 5);
  out(2, "log");
  //fig 2 paper
  std::vector< std::vector<double> > plot_array;
  plot_array.resize(temp_size);
  for(unsigned int i = 0; i < temp_size; i++){
    plot_array[i].resize(3);
    plot_array[i][0] = Temperatura[i]/unidad;
    plot_array[i][1] = Frozen[i][0];
    plot_array[i][2] = Frozen[i][1];
  }
  array_print(plot_array, "Dip_cong.dat");
  out(3, "log");
  //fig 3 paper
  for(unsigned int i = 0; i < temp_size; i++){
    plot_array[i].assign(5,0);
    plot_array[i][0] = Temperatura[i]/unidad;
    plot_array[i][1] = unidad/Temperatura[i];
    plot_array[i][2] = Susceptibilidad[i][1];
    plot_array[i][3] = Susceptibilidad[i][2];
    plot_array[i][4] = Susceptibilidad[i][3];
  }
  array_print(plot_array, "Sus_E0.dat");
  out(4, "log");
  //fig 4 paper
  for(unsigned int i = 0; i < temp_size; i++){
    plot_array[i].assign(5,0);
    plot_array[i][0] = Temperatura[i]/unidad;
    plot_array[i][1] = Frozen[i][1];
    plot_array[i][2] = Frozen[8*temp_size +i][1];
    plot_array[i][3] = Frozen[20*temp_size +i][1];
    plot_array[i][4] = Frozen[22*temp_size +i][1];
  }
  array_print(plot_array, "Dip_cong_E.dat");
  out(5, "log");
  //fig 5 paper
  for(unsigned int i = 0; i < temp_size; i++){
    plot_array[i].assign(4,0);
    plot_array[i][0] = Temperatura[i]/unidad;
    plot_array[i][1] = Susceptibilidad[i][1];
    plot_array[i][2] = Susceptibilidad[8*temp_size +i][1];
    plot_array[i][3] = Susceptibilidad[20*temp_size +i][1];
  }
  array_print(plot_array, "Sus_E.dat");
  out(100, "log");
  //fig 6 ya esta
  //Polarización efriando
  import_data(Polarizacion, "polarizacion.dat", temp_size, field_size);
    for(unsigned int i = 0; i < temp_size; i++){
    plot_array[i].assign(5,0);
    plot_array[i][0] = Temperatura[i]/unidad;
    plot_array[i][1] = Polarizacion[i][0];
    plot_array[i][2] = Polarizacion[i][8];
    plot_array[i][3] = Polarizacion[i][20];
    plot_array[i][4] = Polarizacion[i][22];
  }
  array_print(plot_array, "Cool_pol.dat");
}
