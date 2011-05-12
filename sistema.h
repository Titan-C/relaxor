#ifndef SISTEMA_H
#define SISTEMA_H

#include <vector>

class Sistema
{
  private:
    std::vector<int> sigma;			// Arreglo de spines dipolares
    std::vector<int> sum_sigma;			// Suma de spines individuales
  
    std::vector< std::vector<int> > G;		// Topología de primeros vecinos
    std::vector< std::vector<double> > J;	// Energías de intercambio primeros vecinos
    
  public:
    // Constructor y destructor
    Sistema(unsigned int lado, unsigned int dimension = 3, bool polarizar = true);
    ~Sistema();
    
    
    
};

#endif // SISTEMA_H
