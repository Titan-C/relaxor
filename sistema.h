#ifndef SISTEMA_H
#define SISTEMA_H

#include <vector>
#include <gsl/gsl_rng.h>
#include <iostream>

class Sistema
{
private:
    unsigned int dimension, L;			// Dimensionalidad del sistema, cantidad de PNR por lado
    double DeltaJ, rho;
    std::vector<int> sigma;			// Arreglo de spines dipolares
    std::vector<double> mu_E;			// Arreglo de proyeciones de momento al eje del campo

    std::vector< std::vector<double> > J;	// Energías de intercambio primeros vecinos
    std::vector< std::vector<unsigned int> > G;	// Configuración espacial de primeros vecinos

public:
    // Constructor y destructor
    Sistema(unsigned int lado,
	    gsl_rng * rng,
	    unsigned int dim = 3,
	    bool polarizar = true);
    ~Sistema();

    /*localiza a las PNR generando un vector espacial.
     * Luego encuentra los primeros vecinos*/
    void set_space_config();
    //Llena la matriz de energías de intercambio del sistema
    double Jex(gsl_rng* rng);
    //Genera la polarización inicial
    double set_pol(gsl_rng* rng, bool polarizar);
    //Genera los momentos dipolares
    double set_mu(gsl_rng* rng, bool polarizar);
    //Inicializa al sistema
    double init(gsl_rng* rng, double DJ = 1, double p=0,
		bool polarizar = true, bool write = false);

    //Calcula la energía total del sistema
    double total_E (double E);
    //Calcula la variación de energía del sistema debído a un cambio del spin dipolar
    double delta_E (unsigned int idflip, double E);
    //Calcula la polarización ponderada del sistema en un instante dado
    double norm_pol ();
    /*Evalúa el comportamiento del material a T[temperatura], E(t)[Campo alterno]
     * dados durante el tiempo otorgado(t) en [MCS/dipolo].
     * Graba los datos de ser necesario.*/
    void experimento(double T, double E, unsigned int tau, unsigned int Niter,
		    bool grabar, gsl_rng* rng, std::string id_proc);
    
    //Funciones Para tratar los experimentos del sistema
    void Gen_exp(std::vector<double>& temperaturas, std::vector<double>& campos,
		 std::vector<double> tau, unsigned int numexps, double DJ, double p,
		 unsigned int Equi_iter, unsigned int Exp_iter,std::string Exp_ID, gsl_rng* rng);
    //Variar temperatura del sistema, campos constantes. [Múltiples experimentos]
    void Var_Temp(std::vector<double>& temperaturas, std::vector<double>& campos,
		       std::vector<double> tau, unsigned int numexps, double DJ, double p,
		       unsigned int Equi_iter, unsigned int Exp_iter, gsl_rng* rng);
    //Variar Campo del sistema, temperaturas constantes. [Múltiples experimentos]
    void Var_Field(std::vector<double>& temperaturas, std::vector<double>& campos,
		       std::vector<double> tau, unsigned int numexps, double DJ, double p,
		       unsigned int Equi_iter, unsigned int Exp_iter, gsl_rng* rng);
    //Lazos de histéresis a temperaturas fijas
    void Hist_loop(std::vector<double>& temperaturas, double max_field, unsigned int numexps,
		   double DJ, double p,
		   unsigned int Equi_iter, unsigned int Exp_iter, gsl_rng* rng);
    
};

//Operaciones necesarias para tratar datos
/*PreProcesa los datos almacenados de los experimentos realizados de acuerdo a sus
 * condicones y número de ejecuciones.*/
void pp_data(std::vector<double>& pol_stats, std::vector<double>& pol_int,
	     unsigned int data_length, unsigned int numexps, unsigned int tau,
	     unsigned int Niter, std::string id_proc);
/*Encuentra la polarización promedio de manera absoluta o respetando signo, incluye 
 * también la magnitud del error*/
void eval_pol(const std::vector< double >& pol_stats, unsigned int numexps,
	      double unidad, const std::vector< double >& x_array,
	      std::string id_proc, bool absolut);
/*Calcula la susceptibilidad del sistema*/
void calc_sus(const std::vector<double>& pol_int_avg, unsigned int numexps,
	      double unidad, const std::vector<double>& x_array,
	      const std::vector<double>& campo, std::string id_proc);

//Funciones adicionales
//Calcula la desviación estandar de los datos de toda la matriz
double stan_dev(const std::vector< std::vector< double > >& M);
//Genera un vector de datos double, entre dos números con cierto paso
std::vector<double> step2vec(double unidad, double v_start, double v_end, double dv, std::vector< double > last);
//Genera un vector de datos para un lazo
std::vector<double> loop2vec(double unidad, double max, int divs);
//Genera un vector de datos double que contiene los números declarados en el string
std::vector<double> str2vec(double unidad, std::string magnitudes);
//Realiza una integración por Simpson de la función f con un peso
double simpson_int(const double f_array[], const std::vector< double >& weight);
//Genera vectores de ondas cos- senoidales
std::vector<double> waves(unsigned int length, unsigned int tau, double amplitude, bool cossin);

#endif // SISTEMA_H
