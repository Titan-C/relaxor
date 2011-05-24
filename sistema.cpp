#include "sistema.h"
#include <cmath>
#include "impresor.h"
 #include <gsl/gsl_statistics.h>

/*Constructor:
Dimensiona y encera a los vectores del sistema. Llena sus datos iniciales */
Sistema::Sistema(unsigned int lado,
		 unsigned int Niter,
		 gsl_rng* rng,
		 double r_max,
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
  DeltaJ = init(rng, r_max, polarizar);
  
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

double Sistema::init(gsl_rng* rng, double r_max, bool polarizar){
  //Creo variable auxiliares
  unsigned int ind_xy, L2=L*L;
  double mu_mag, theta, mstheta, phi, PI = 4*atan(1);
  std::vector< std::vector<double> > mu, R;
  mu.resize(sigma.size());
  R.resize(sigma.size());
  
  // Creo las configuraciones de las Nanoregiones Polares (PNR)
  for(unsigned int i=0; i<mu.size(); i++){
   
    // Crear variables descriptivas de PNR
    mu_mag = r_max*gsl_rng_uniform(rng);
    mu_mag = 4*PI/3 * mu_mag*mu_mag*mu_mag * P;
    theta  = PI/2;
    theta *= gsl_rng_uniform(rng);
    mstheta = mu_mag * sin(theta);
    phi    = 2*PI*gsl_rng_uniform(rng);
    
    //Transformar en momentos dipolares
    mu[i].resize(dimension);
    
    mu[i][0] = mstheta * std::cos(phi);
    mu[i][1] = mstheta * std::sin(phi);
    mu[i][2] = std::cos(theta);
    
    // Almacenar Datos del sistema
    sigma[i] = (mu[i][2] > 0) ? 1 : -1 ;
    mu_H[i]  = std::abs(mu[i][2]);
    
    // Coeficientes vector posición i-ésima PNR para base 2*r_max
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
  array_print(mu, "mu_PNR.dat");
  array_print(R, "posiciones.dat");
  array_print(G, "grafo_vecinos.dat");
  // Calcular las enrgías de intercambio de las PNR
  double r, Jex;
  std::vector<double> delta_R;
  std::vector< std::vector<double> > Jinter;
  delta_R.resize(dimension);
  Jinter.resize(sigma.size());
  for(unsigned int i = 0; i<Jinter.size(); i++){
    Jinter[i].resize(sigma.size());
    for(unsigned int j = i+1; j<Jinter[i].size(); j++){
      /* Calcula vector diferencia entre PNR,
       * Aplica condiciones de borde,
       * Calcula el módulo */
      for(unsigned int k=0; k<dimension; k++)
	delta_R[k] = R[j][k]-R[i][k];
      condborde(delta_R,L);
      std::cout<<"("<<delta_R[0]<<", "<<delta_R[1]<<", "<<delta_R[2]<<") ";
      r=2*r_max*::sqrt(dot(delta_R,delta_R));
      std::cout<<r<<"|";
      // Calcula la Energía de intercambio entre 2 PNR
      Jex = 3*(4*r_max*r_max)*dot(mu[i],delta_R)*dot(mu[j],delta_R)/(r*r*r*r*r)-dot(mu[i],mu[j])/(r*r*r);
      Jinter[i][j] = Jex/2;      
    }
    std::cout<<std::endl;
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
 
  return stan_dev(Jinter);
  
  Jinter.clear();
  mu.clear();
  R.clear();
  delta_R.clear();
}
//Calcular la desviación standar del las energías de intercambio.
double Sistema::stan_dev(const std::vector< std::vector<double> >& M)
{
  unsigned int celdas, columnas;
  columnas = M[1].size();
  celdas = M.size() * columnas;
  double * Jij;
  Jij = new double [celdas];
  for(unsigned int i = 0 ; i<M.size(); i++){
    for(unsigned int j = 0; j<columnas; j++)
      Jij[i*columnas + j] = J[i][j];
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
  T *= DeltaJ;
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
  for(unsigned int i = 0; i < sum_sigma_time.size(); i++)
    sum_sigma_time[i] = 0;
  for(unsigned int i = 0; i < sum_sigma_conf.size(); i++)
    sum_sigma_conf[i] = 0;    
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
void temp_array(std::vector< double >& Temperatura, double T, double dT)
{
  Temperatura.resize(int (T/dT));
  for(unsigned int i = 0; i<Temperatura.size() ;i++)
    Temperatura[i] = T - i*dT;
}
void field_array(std::vector< double >& campo)
{
  int Mediciones = 1;
  campo.resize(Mediciones);
  campo[0]=0;
}


void procesar(unsigned int Niter, unsigned int L, const std::vector<double> & Temperatura)
{
  //Encontrar proporción de dipolos congelados
  std::vector< std::vector<double> > sigmas_time, S_frozen;
  unsigned int lentos = 4, L3 = L*L*L;

  import_data(sigmas_time, "sum_sigma_time.dat", Temperatura.size() , L3);

  S_frozen.resize(Temperatura.size());
  for(unsigned int i = 0; i < S_frozen.size(); i++){
    S_frozen[i].resize(lentos+1);
    S_frozen[i][lentos] = Temperatura[i];
    for(unsigned int j = 0; j < lentos; j++)
      S_frozen[i][j] = 0;
  }
  for(unsigned int i=0 ; i<sigmas_time.size(); i++){
    for(unsigned int j=0; j<sigmas_time[i].size(); j++){
      sigmas_time[i][j] = std::abs(sigmas_time[i][j]/Niter);
      for(unsigned int k=0;k<lentos;k++){
	if (sigmas_time[i][j] >= (1-0.1*k) )
	  S_frozen[i][k]+=(double) 1/L3;
      }
    }
  }
  array_print(S_frozen, "Congelamiento.dat");

  //Encontrar la susceptibilidad
  std::vector< std::vector<double> > Susceptibilidad;
  Susceptibilidad.resize(Temperatura.size());
  for(unsigned int i=0; i<Susceptibilidad.size(); i++){
    Susceptibilidad[i].resize(lentos+1);
    Susceptibilidad[i][lentos] = Temperatura[i];
    for(unsigned int j=0; j<lentos; j++)
      Susceptibilidad[i][j] = (1 - S_frozen[i][j])/Temperatura[i];
  }
  array_print(Susceptibilidad, "Susceptibilidad.dat");

  sigmas_time.clear();
  S_frozen.clear();
  Susceptibilidad.clear();
}