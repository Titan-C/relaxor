#ifndef IMPRESOR_H
#define IMPRESOR_H

#include <iostream>
#include <vector>

// Imprime datos de variables double
void out(double value, std::string ARCHIVO, bool app = true, bool br = true);
// Imprime datos de los arreglos vectoriales
void array_print_bin(int * V, unsigned int size, std::string ARCHIVO, bool app = true);
void array_print(const std::vector< int >& V, std::string ARCHIVO, bool app = true, bool br = true);
void array_print(const std::vector< double >& V, std::string ARCHIVO, bool app = true, bool br = true);
void array_print_bin(const std::vector< double >& V, std::string ARCHIVO, bool app = true);
// Imprime datos de los arreglos matricales
void array_print(const std::vector< std::vector<int> >& M,
		 std::string ARCHIVO, bool app=false, bool br = true);
void array_print(const std::vector< std::vector<unsigned int> >& M,
		 std::string ARCHIVO, bool app=false, bool br = true);
void array_print(const std::vector< std::vector<double> >& M,
		 std::string ARCHIVO, bool app=false, bool br = true);
void array_print_bin(const std::vector< double * >& V, std::string ARCHIVO, unsigned int cols);
void import_data(std::vector < std::vector< double > >& M,
		 std::string ARCHIVO,
                 unsigned int filas, unsigned int columnas);
#endif // IMPRESOR_H
