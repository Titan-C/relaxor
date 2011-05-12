#ifndef SISTEMA_H
#define SISTEMA_H

#include <vector>
#include <gsl/gsl_rng.h>

class Sistema
{
  private:
    unsigned int dimension, L;
    double Pol=1;
    
    
    std::vector<int> sigma;			// Arreglo de spines dipolares
    std::vector<int> sum_sigma;			// Suma de spines individuales
    std::vector<double> mu_H;			// Arreglo de proyeciones de momento al eje del campo
  
    std::vector< std::vector<int> > G;		// Topología de primeros vecinos
    std::vector< std::vector<double> > J;	// Energías de intercambio primeros vecinos

  public:
    // Constructor y destructor
    Sistema(unsigned int lado, gsl_rng * rng, unsigned double r_max = 2, unsigned int dim = 3, bool polarizar = true);
    ~Sistema();
    
    void init(gsl_rng* rng, unsigned double r_max, bool polarizar);
    
    
    
    
    
};
// Operaciones necesarias para tratar al sistema
void condborde ( vector <double>& R, int L);
double dot(const vector<double>& a, const vector<double>& b);

#endif // SISTEMA_H
