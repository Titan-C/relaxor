#ifndef IMPRESOR_H
#define IMPRESOR_H

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <sys/types.h>
#include <stdint.h>
#include <sys/stat.h>

// Imprime datos de variables double
void out(double value, std::string ARCHIVO, bool br = true);
void array_print(const std::vector< double >& V, std::string ARCHIVO, unsigned int colsize, double scale);
void array_print_bin(const std::vector< double * >& V, std::string ARCHIVO, unsigned int cols);
void import_data(std::vector < std::vector< double > >& M,
		 std::string ARCHIVO,
                 unsigned int filas, unsigned int columnas);


// Imprime arreglos vectoriales de todo tipo de variable en texto y binario
template <typename T>
void array_print(const std::vector<T> & V, std::string ARCHIVO) {
  std::ofstream file;
  file.open (ARCHIVO.c_str(), std::fstream::app);
  
  for(unsigned int i = 0; i<V.size(); i++)
    file<<V[i]<<"\t";
  file<<"\n";
  file.close();
}
template <typename T>
void array_print_bin(const std::vector< T >& V, std::string ARCHIVO) {
  std::ofstream file;
  file.open (ARCHIVO.c_str(), std::fstream::app);
  file.write((char * )&V[0],V.size()*sizeof(T));
  file.close();
}
// Imprime arreglos matriciales de todo tipo de variable en texto y binario
template <typename T>
void array_print(const std::vector< std::vector<T> >& M, std::string ARCHIVO) {
  std::ofstream file;
  file.open (ARCHIVO.c_str());

  for (unsigned int i=0;i<M.size();i++) {
    for (unsigned int j=0;j<M[i].size();j++) {
      file<<M[i][j]<<"\t";
    }
    file<<"\n";
  }
  file.close();
}

//Numpy Arrays
static const char LIBNPY_VERSION[] = "0.5";
static const char MAGIC[] = "\x93NUMPY";
static const unsigned int  MAJOR = 1;
static const unsigned int  MINOR = 0;
static const unsigned int  MAX_HDR_LEN = 256 * 256;
static const unsigned int  MAX_INT_LEN = 32;
static const unsigned int  PREAMBLE_LEN = 6 + 1 + 1 + 2;

#if __BYTE_ORDER == __LITTLE_ENDIAN
static const char ENDIAN_CHAR = '<';
#else
static const char ENDIAN_CHAR = '>';
#endif

int create_metadata(char preamble[PREAMBLE_LEN], char header[MAX_HDR_LEN],
                    char* descr, int fortran_order,
		    const std::vector<unsigned int> shape);

void npy_save(char* fname, char* descr, int fortran_order,
              const std::vector<unsigned int> shape, size_t sz, void* data);

void npy_save_double(char* fname, int fortran_order,
                     const std::vector<unsigned int> shape, double* data);


#endif // IMPRESOR_H
