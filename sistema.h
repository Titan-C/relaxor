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
    void flip (unsigned int idflip, double T, double E, gsl_rng* rng);//Realiza el cambio del spin dipolar en una ubicación dada

    int experimento(double T, double E, double tau, unsigned int Niter,
		    bool grabar, gsl_rng* rng, std::string id_proc);
};
// Operaciones necesarias para tratar al sistema
double stan_dev(const std::vector< std::vector< double > >& M);	//Calsula la desviación estandar de los datos de toda la matriz

std::vector<double> step2vec(double unidad, double T_top, double dT, bool heat);
std::vector<double> str2vec(double unidad, std::string magnitudes);
void calc_sus(unsigned int numexps, unsigned int tau, unsigned int Niter, unsigned int L, double unidad,
	      const std::vector< double >& Temperatura, const double campos, std::string id_proc);
void eval_pol(unsigned int Niter, unsigned int numexps, double unidad, const std::vector< double >& Temperatura, std::string id_proc);
void plot_data_sus(double unidad, const std::vector<double>& Temperatura, const std::vector<double>& campos, std::string id_proc);

#endif // SISTEMA_H
