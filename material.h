#ifndef MATERIAL_H
#define MATERIAL_H

#include "experiment.h"
#include "impresor.h"
#include <vector>
#include <ctime>
#include <cstdio>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

class Material
{
private:
  //Attributes

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

  //Methods

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

  // Evaluates system's total Energy
  double total_E (double E);
  // Evaluates Energy change after dipole flip
  double delta_E (unsigned int idflip, double E);
  // Evaluate Global polarization
  double norm_pol ();


  // Attemps to flip every dipole within the array acoording to Boltzman Probability
  void MonteCarloStep(double T, double E_field);
  /* Internal material operations */
  // Update material state log
  void update_log_sigma(std::vector< int >& log_sigma);

public:
  // Constructor & destructor
  Material(unsigned int L,
	   bool polarizar = true);
  ~Material();


  //Initialize material
  void init(double p, std::string ID, bool polarizar = true, bool write = false);

  /* Evaluates material behavior given T[temperature] and E(t)[external field]
   * during the given amount of time Niter. Records data if requested */
  void state(double T, std::vector< double >& field, bool record);

  // To allow test on material
  friend class tester;

};

#endif // MATERIAL_H