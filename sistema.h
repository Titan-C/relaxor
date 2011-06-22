#ifndef SISTEMA_H
#define SISTEMA_H

#include <vector>
#include <iostream>
#include <gsl/gsl_rng.h>

class Sistema
{
private:
    unsigned int dimension, L;
    double P;


    std::vector<int> sigma;			// Arreglo de spines dipolares
    std::vector<int> sum_sigma_time;		// Suma de spines individuales en el tiempo
    std::vector<int> sum_sigma_conf;		// Suma de spines 
    std::vector<double> mu_H;			// Arreglo de proyeciones de momento al eje del campo

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


    double init(gsl_rng* rng, double Delta_J, bool polarizar);
    // Operaciones de ejecución
    double total_E (double E);
    double delta_E (unsigned int idflip, double E);
    void flip (unsigned int idflip, double T, double E, gsl_rng* rng);
    // Operaciones de experimento
    int experimento(double T, double E, unsigned int Niter,
		    bool grabar, gsl_rng* rng);
    void reset_sum_sigma();
};
// Operaciones necesarias para tratar al sistema
void condborde ( std::vector <double>& R, int L);
double dot(const std::vector<double>& a, const std::vector<double>& b);
double stan_dev(const std::vector< std::vector< double > >& M);

void temp_array(std::vector< double >& Temperatura, double unidad, double T, double dT);
void field_array(std::vector< double >& campo, double unidad);
void procesar(unsigned int numexps, unsigned int Niter, unsigned int L, double unidad,
	      const std::vector< double >& Temperatura, const std::vector< double >& campos);
void graficos(double unidad, const std::vector< double >& Temperatura, const std::vector< double >& campos);
#endif // SISTEMA_H
