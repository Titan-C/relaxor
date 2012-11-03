#ifndef EXPERIMENT_H
#define EXPERIMENT_H


#include "material.h"

#define _2pi 8*atan(1)

//Perform Experiment according to certain parameters
void doExperiment(unsigned int repetitions,		//
		  unsigned int Equilibration_Iter,	//
		  std::vector<double>& Temperature_loop,//Temperature range
		  std::vector<double>& Electric_Field,	//External Electric Field
		  Material & relaxor,
		  bool polarize = false);

std::string ExpLabel(std::string Exp_ID,unsigned int L, double rho, double Field, double tau, unsigned int numexps,
		     std::vector<double>& Tdat, unsigned int Exp_iter, unsigned int Equi_iter);

void Gen_exp(unsigned int L, unsigned int numexps, std::vector<double> rho, std::vector<double>& Tdat,
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

//Genera un vector de datos double, entre dos números con cierto paso
void step2vec(double v_start, double v_end, double dv, std::vector< double >& array, double unidad = 1);
//Genera un vector de datos para un lazo
std::vector<double> loop2vec(double max, int divs, double unidad=1);
//Genera un vector de datos double que contiene los números declarados en el string
std::vector<double> str2vec(std::string magnitudes, double unidad=1);
//Devuelve la cantidad de pasos necesarios durante la simulacion para cada frecuencia
unsigned int stepEstimator(unsigned int Niter, unsigned int tau, unsigned int min_periods);
//Realiza una integración por Simpson de la función f con un peso
double simpson_int(const double f_array[], const std::vector< double >& weight);
//Genera vectores de ondas cos- senoidales
std::vector<double> wavearray(double amplitude, unsigned int tau, unsigned int length, double phase);
//Validador que justifica realizar una simulacion
bool needSimulation(std::string id_proc, unsigned int iter, unsigned int PNR, unsigned int temps, unsigned int numexps);

#endif // EXPERIMENT_H
