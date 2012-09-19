#ifndef SISTEMA_H
#define SISTEMA_H

#include <vector>
#include <gsl/gsl_rng.h>
#include <iostream>

class Material
{
private:
  gsl_rng * rng;	// GSL Random number Generator
  std::string ExpID;	// Experiment identification

  unsigned int PNR;	// Total PNRs
  double rho;		// Mean Ferroelectricity

  std::vector<int> sigma;	// Dipolar "Spin States"
  std::vector<double> mu_E;	// Projected Dipole magnitude on main axis

  // Topological configuration of PNRs
  std::vector< std::vector<unsigned int> > G;
  // Exchange energies between PNRs
  std::vector< std::vector<double> >  J;

public:
    // Constructor & destructor
    Material(unsigned int L,
	    bool polarizar = true);
    ~Material();

    /*Setup Topological configuration of
     * PNRs inside the material */
    void set_space_config(unsigned int L);
    /* Fills values for the interaction energy between
     * PNRs according to spatial configuration */
    void Jex();
    //Generate initial polarization
    void set_pol(bool polarize);
    //Generate dipolar moments
    void set_mu(bool polarize);
    //Initialize material
    void init(double p, std::string ID = 0, bool polarizar = true, bool write = false);

    // Evaluates system's total Energy
    double total_E (double E);
    // Evaluates Energy change after dipole flip
    double delta_E (unsigned int idflip, double E);
    // Evaluate Global polarization
    double norm_pol ();
    
    void update_log_sigma(std::vector< double >& log_sigma);
    // Attemps to flip every dipole within the array acoording to Boltzman Probability
    void MonteCarloStep(double T, double E_field);
    /*Evalúa el comportamiento del material a T[temperatura], E(t)[Campo alterno]
     * dados durante el tiempo otorgado(t) en [MCS/dipolo].
     * Graba los datos de ser necesario.*/
    void experimento(double T, double E, unsigned int tau, unsigned int Niter,
		    bool grabar);
    void flip_sigma(unsigned int idsigma);
    unsigned int return_PNR();
    int ret_sig(unsigned int i);
    std::vector<int> ret_sigarr();
};
//Funciones Para tratar los experimentos del sistema
void Gen_exp(unsigned int L, unsigned int numexps,std::vector<double> rho, std::vector<double>& Tdat,
	     std::vector<double>& Fields, std::vector<double> tau, std::string Exp_ID);
//Operaciones necesarias para tratar datos
void proces_data(std::vector< double >& Temps, double Field,
		 unsigned int tau, unsigned int numexps, unsigned int PNR,
		 double p, unsigned int Niter, std::string id_proc);
/*PreProcesa los datos almacenados de los experimentos realizados de acuerdo a sus
 * condicones y número de ejecuciones.*/
void pp_data(std::vector<double>& pol_stats, std::vector<double>& pol_int,
	     unsigned int data_length, unsigned int numexps, unsigned int tau,
	     unsigned int Niter, std::string id_proc);
/*Encuentra la polarización promedio de manera absoluta o respetando signo, incluye 
 * también la magnitud del error*/
void eval_pol(const std::vector< double >& pol_stats, unsigned int numexps,
	      const std::vector< double >& x_array,
	      std::string id_proc, bool absolut);
/*Calcula la susceptibilidad del sistema*/
void calc_sus(const std::vector<double>& pol_int_avg, unsigned int numexps,
	      const std::vector<double>& x_array,
	      const std::vector<double>& campo, std::string id_proc);
/*Calcula la polarización congelada*/
void eval_frozen(unsigned int array_size, unsigned int Niter, const std::vector<double>& Temps,
		 unsigned int numexps, std::string id_proc);

//Funciones adicionales
//Calcula la desviación estandar de los datos de toda la matriz
double stan_dev(std::vector< std::vector< double > > M);
//Genera el arreglo de temperaturas
std::vector<double> thermostat(unsigned int n, unsigned int i, double rho, double dT, double Tf);
//Genera un vector de datos double, entre dos números con cierto paso
std::vector<double> step2vec(double v_start, double v_end, double dv, std::vector< double > last, double unidad = 1);
//Genera un vector de datos para un lazo
std::vector<double> loop2vec(double max, int divs, double unidad=1);
//Genera un vector de datos double que contiene los números declarados en el string
std::vector<double> str2vec(std::string magnitudes, double unidad=1);
//Devuelve la cantidad de pasos necesarios durante la simulacion para cada frecuencia
unsigned int stepEstimator(unsigned int Niter, unsigned int tau, unsigned int min_periods);
//Realiza una integración por Simpson de la función f con un peso
double simpson_int(const double f_array[], const std::vector< double >& weight);
//Genera vectores de ondas cos- senoidales
std::vector<double> cosarray(unsigned int length, unsigned int tau, double amplitude,  double phase);
//Validador que justifica realizar una simulacion
bool needSimulation(std::string id_proc, unsigned int size);

#endif // SISTEMA_H
