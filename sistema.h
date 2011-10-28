#ifndef SISTEMA_H
#define SISTEMA_H

#include <vector>
#include <gsl/gsl_rng.h>
#include <iostream>

class Sistema
{
private:
    unsigned int dimension, L;			// Dimensionalidad del sistema, cantidad de PNR por lado

    			// Arreglo de spines dipolares
    std::vector<double> mu_E;			// Arreglo de proyeciones de momento al eje del campo

    std::vector< std::vector<double> > J;	// Energías de intercambio primeros vecinos
    std::vector< std::vector<unsigned int> > G;	// Configuración espacial de primeros vecinos

public:
  std::vector<int> sigma;
    double DeltaJ;
    // Constructor y destructor
    Sistema(unsigned int lado,
	    gsl_rng * rng,
	    double Delta_J = 1,
	    unsigned int dim = 3,
	    bool polarizar = true);
    ~Sistema();


    double init(gsl_rng* rng, double Delta_J, bool polarizar);	//Inicializa al sistema
    double set_pol(gsl_rng* rng, bool polarizar);			//Genera la polarización inicial
    void set_mu(gsl_rng* rng);			//Genera los momentos dipolares
    double Jex(gsl_rng* rng, double Delta_J);	//Llena la matriz de energías de intercambio del sistema
    void set_space_config();			//localiza a las PNR generando un vector espacial. Luego encuentra los primeros vecinos


    double total_E (double E);					//Calcula la energía total del sistema
    double delta_E (unsigned int idflip, double E);		//Calcula la variación de energía del sistema debído a un cambio del spin dipolar
    double norm_pol ();						//Calcula la polarización ponderada del sistema en un instante dado

    int experimento(double T, double E, unsigned int tau, unsigned int Niter,
		    bool grabar, gsl_rng* rng, std::string id_proc);/* Realiza la simulación  del material a una temperatura dada. Para un campo alteno de amplitud E y periodo tau [MCS/dipolo] */
};
// Operaciones necesarias para tratar al sistema
double stan_dev(const std::vector< std::vector< double > >& M);	//Calsula la desviación estandar de los datos de toda la matriz

std::vector<double> step2vec(double unidad, double T_top, double dT, bool heat);
std::vector<double> str2vec(double unidad, std::string magnitudes);
void pp_data(std::vector<double>& pol_stats, std::vector<double>& pol_int, unsigned int data_length,
	      unsigned int numexps, unsigned int tau, unsigned int Niter, std::string id_proc);
void calc_sus(const std::vector<double>& pol_int_avg, unsigned int numexps, double unidad,
	      const std::vector<double>& x_array, const std::vector<double>& campo, std::string id_proc);
void eval_pol(const std::vector< double >& pol_stats, unsigned int numexps, double unidad, const std::vector< double >& x_array, std::string id_proc, bool absolut);
double simpson_int(const double f_array[], const std::vector< double >& weight);
std::vector<double> waves(unsigned int length, unsigned int tau, double amplitude, bool cossin);
#endif // SISTEMA_H
