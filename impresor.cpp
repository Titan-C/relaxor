#include "impresor.h"

void out(double value, std::string ARCHIVO) {
  std::ofstream file;
  file.open (ARCHIVO.c_str(), std::fstream::app);
  file<<value<<"\n";
  file.close();
}
// Imprime datos de los arreglos vectoriales
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


// Imprime datos de los arreglos matricales
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
