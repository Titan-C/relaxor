#ifndef SISTEMA_H
#define SISTEMA_H

#include <vector>
#include <gsl/gsl_rng.h>
#include <iostream>

class Sistema
{
private:
    unsigned int dimension, L;

    std::vector<int> sigma;			// Arreglo de spines dipolares
    std::vector<int> sum_sigma_time;		// Suma de spines individuales en el tiempo
    std::vector<int> sum_sigma_conf;		// Suma de spines 
    std::vector<double> mu_E;			// Arreglo de proyeciones de momento al eje del campo

    std::vector< std::vector<unsigned int> > G;	// Topología de primeros vecinos
    std::vector< std::vector<double> > J;	// Energías de intercambio primeros vecinos

public:
    double DeltaJ;
    // Constructor y destructor
    Sistema(unsigned int lado,
	    unsigned int Niter,
	    gsl_rng * rng,
	    double Delta_J = 1.2,
	    unsigned int dim = 3,
	    bool polarizar = true);
    ~Sistema();


    double init(gsl_rng* rng, double Delta_J, bool polarizar);	//Inicializa al sistema
    int init_pol(gsl_rng* rng, bool polarizar);			//Genera la polarización inicial
    void mul(gsl_rng* rng, bool print);				//Genera los momentos dipolares
    double Jex(gsl_rng* rng, double Delta_J, bool print);	//Llena la matriz de energías de intercambio del sistema
    void topografia(bool print);				//localiza a las PNR generando un vector espacial. Luego encuentra los primeros vecinos



    double total_E (double E);					//Calcula la energía total del sistema
    double delta_E (unsigned int idflip, double E);		//Calcula la variación de energía del sistema debído a un cambio del spin dipolar
    void flip (unsigned int idflip, double T, double E, gsl_rng* rng);//Realiza el cambio del spin dipolar en una ubicación dada

    int experimento(double T, double E, unsigned int Niter,
		    bool grabar, gsl_rng* rng, std::string id_proc);
    void reset_sum_sigma();
};
// Operaciones necesarias para tratar al sistema
void condborde ( std::vector <double>& R, int L);
double dot(const std::vector<double>& a, const std::vector<double>& b);
double stan_dev(const std::vector< std::vector< double > >& M);

std::vector<double> temp_array(double unidad, double T_top, double dT, bool heat);
std::vector<double> field_array(double unidad, double H_top, double dH);
void Sus_proc(unsigned int numexps, unsigned int Niter, unsigned int L, double unidad,
	      const std::vector< double >& Temperatura, const std::vector< double >& campos, std::string id_proc);
void Pol_proc(unsigned int numexps, unsigned int Niter, unsigned int L, double unidad,
	      const std::vector< double >& Temperatura, const std::vector< double >& campos, std::string id_proc);
void plot_data_sus(double unidad, const std::vector<double>& Temperatura, const std::vector<double>& campos, std::string id_proc);

#endif // SISTEMA_H
