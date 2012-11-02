#ifndef IMPRESOR_H
#define IMPRESOR_H

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <sys/types.h>
#include <typeinfo>
#include <stdint.h>
#include <sys/stat.h>

// Imprime datos de variables double
void out(double value, std::string ARCHIVO, bool br = true);
void array_print(const std::vector< double >& V, std::string ARCHIVO, unsigned int colsize, double scale);

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
template <typename T>
void array_print_bin(const std::vector< std::vector<T> >& M, std::string ARCHIVO){
  std::ofstream file;
  file.open (ARCHIVO.c_str(), std::fstream::app);

  for (unsigned int i=0;i<M.size();i++)
    file.write((char * )&M[i][0],M[i].size()*sizeof(T));

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

void create_metadata(std::string ARCHIVO,char* descr, int fortran_order,
		    const std::vector<unsigned int> shape);

template <typename T>
void npy_save_vector(std::string ARCHIVO, const std::vector<T>& data, int fortran_order=0)
{
    char descr[5];
    descr[0] = ENDIAN_CHAR;
    descr[1] = ( typeid(T)==typeid(double) ) ? 'f' : 'i';
    sprintf(descr+2, "%d", (int) sizeof(T));

    std::vector<unsigned int> shape (1,data.size());

    create_metadata(ARCHIVO,descr,fortran_order,shape);
    array_print_bin(data,ARCHIVO);
}

template <typename T>
void npy_save_matrix(std::string ARCHIVO, const std::vector< std::vector<T> >& data, int fortran_order=0)
{
    char descr[5];
    descr[0] = ENDIAN_CHAR;
    descr[1] = ( typeid(T)==typeid(double) ) ? 'f' : 'i';
    sprintf(descr+2, "%d", (int) sizeof(T));

    std::vector<unsigned int> shape (2,0);
    shape[0] = data.size();
    shape[1] = data[0].size();

    create_metadata(ARCHIVO,descr,fortran_order,shape);
    array_print_bin(data,ARCHIVO);
}
#endif // IMPRESOR_H
