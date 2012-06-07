#include "impresor.h"
#include <fstream>
#include <sstream>
#include <cstdlib>

void out(double value, std::string ARCHIVO, bool app, bool br) {
  std::fstream file;
    if (app)
    file.open (ARCHIVO.c_str(), std::fstream::out | std::fstream::app);
  else
    file.open (ARCHIVO.c_str(), std::fstream::out);

  file<<value<<"\t";
  if (br) file<<"\n";
  file.close();
}
// Imprime datos de los arreglos vectoriales
void array_print(const std::vector< int >& V, std::string ARCHIVO, bool app, bool br) {
  std::fstream file;
  if (app)
    file.open (ARCHIVO.c_str(), std::fstream::out | std::fstream::app);
  else
    file.open (ARCHIVO.c_str(), std::fstream::out);

  for(unsigned int i = 0; i<V.size(); i++)
    file<<V[i]<<"\n";

  if (br) file<<"\n";
  file.close();
}
void array_print_bin(int8_t* V, unsigned int size, std::string ARCHIVO, bool app) {
  std::ofstream file;
  if (app)
    file.open (ARCHIVO.c_str(), std::fstream::app);
  else
    file.open (ARCHIVO.c_str());
  
  file.write((char * )&V[0],size*sizeof(int8_t));
  file.close();
}
void array_print(double * V, unsigned int size, std::string ARCHIVO, bool app) {
  std::ofstream file;
  if (app)
    file.open (ARCHIVO.c_str(), std::fstream::app);
  else
    file.open (ARCHIVO.c_str());
  
  for(unsigned int i = 0; i<size; i++)
    file<<V[i]<<"\n";
  file.close();
}
void array_print(const std::vector< double >& V, std::string ARCHIVO, bool app, bool br) {
  std::fstream file;
  if (app)
    file.open (ARCHIVO.c_str(), std::fstream::out | std::fstream::app);
  else
    file.open (ARCHIVO.c_str(), std::fstream::out);

  for(unsigned int i = 0; i<V.size(); i++)
    file<<V[i]<<"\n";

  if (br) file<<"\n";
  file.close();
}
void array_print(const std::vector< double >& V, std::string ARCHIVO, unsigned int colsize, double scale){
  std::fstream file;
  file.open (ARCHIVO.c_str(), std::fstream::out);
  unsigned int rowsize = V.size()/colsize;

  for(unsigned int i = 0; i<rowsize; i++){
    for(unsigned int j = 0; j<colsize;j++)
      file<<V[i*colsize+j]/scale<<" ";
    file<<"\n";
  }
  file.close();
}
void array_print_bin(const std::vector< double >& V, std::string ARCHIVO, bool app) {
  std::ofstream file;
  if (app)
    file.open (ARCHIVO.c_str(), std::fstream::app);
  else
    file.open (ARCHIVO.c_str());
  
  file.write((char * )&V[0],V.size()*sizeof(double));
  file.close();
}
// Imprime datos de los arreglos matricales
void array_print(const std::vector< std::vector<int> >& M, std::string ARCHIVO, bool app, bool br) {
  std::fstream file;
  if (app)
    file.open (ARCHIVO.c_str(), std::fstream::out | std::fstream::app);
  else
    file.open (ARCHIVO.c_str(), std::fstream::out);

  for (unsigned int i=0;i<M.size();i++) {
    for (unsigned int j=0;j<M[i].size();j++) {
      file<<M[i][j]<<"\t";
    }
    if (br) file<<"\n";
  }
  file.close();
}


void array_print(unsigned int ** M, std::string ARCHIVO, unsigned int len) {
  std::fstream file;
  file.open (ARCHIVO.c_str(), std::fstream::out);

  for (unsigned int i=0;i<len;i++) {
    for (unsigned int j=0;j<6;j++) {
      file<<M[i][j]<<"\t";
    }
    file<<"\n";
  }
  file.close();
}

void array_print(double ** M, std::string ARCHIVO, unsigned int len) {
  std::fstream file;
  file.open (ARCHIVO.c_str(), std::fstream::out);

  for (unsigned int i=0;i<len;i++) {
    for (unsigned int j=0;j<6;j++) {
      file<<M[i][j]<<"\t";
    }
    file<<"\n";
  }
  file.close();
}

void array_print(const std::vector< std::vector< unsigned int > >& M, std::string ARCHIVO, bool app, bool br) {
  std::fstream file;
  if (app)
    file.open (ARCHIVO.c_str(), std::fstream::out | std::fstream::app);
  else
    file.open (ARCHIVO.c_str(), std::fstream::out);

  for (unsigned int i=0;i<M.size();i++) {
    for (unsigned int j=0;j<M[i].size();j++) {
      file<<M[i][j]<<"\t";
    }
    if (br) file<<"\n";
  }
  file.close();
}
void array_print(const std::vector< std::vector< double > >& M, std::string ARCHIVO, bool app, bool br) {
  std::fstream file;
  if (app)
    file.open (ARCHIVO.c_str(), std::fstream::out | std::fstream::app);
  else
    file.open (ARCHIVO.c_str(), std::fstream::out);

  for (unsigned int i=0;i<M.size();i++) {
    for (unsigned int j=0;j<M[i].size();j++) {
      file<<M[i][j]<<"\t";
    }
    if (br) file<<"\n";
  }
  file.close();
}
void array_print_bin(const std::vector< double * >& V, std::string ARCHIVO, unsigned int cols){
  std::ofstream file;
  file.open (ARCHIVO.c_str());
  
  for (unsigned int i=0;i<V.size();i++)
    file.write((char * )&V[i][0],cols*sizeof(double));
  
  file.close();
}
void import_data(std::vector< std::vector< double > >& M, std::string ARCHIVO, unsigned int filas, unsigned int columnas) {
  std::fstream file;
  file.open (ARCHIVO.c_str(), std::fstream::in);
  M.resize(filas);
  for(unsigned int i=0;i<filas;i++){
    M[i].resize(columnas);
    for(unsigned int j=0;j<columnas;j++)
      file>>M[i][j];
  }
  file.close();
}
