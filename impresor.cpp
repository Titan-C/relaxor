#include "impresor.h"
#include <fstream>

void out(double value, std::string ARCHIVO, bool app, bool br) {
  std::fstream file;
    if (app)
    file.open (ARCHIVO.c_str(), std::fstream::out | std::fstream::app);
  else
    file.open (ARCHIVO.c_str(), std::fstream::out);

  file<<value;
  if (br); file<<std::endl;
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
    file<<V[i]<<"\t";

  if (br); file<<std::endl;
  file.close();
}
// Imprime datos de los arreglos matricales
void array_print(const std::vector< std::vector< int > >& M, std::string ARCHIVO, bool app, bool br) {
  std::fstream file;
  if (app)
    file.open (ARCHIVO.c_str(), std::fstream::out | std::fstream::app);
  else
    file.open (ARCHIVO.c_str(), std::fstream::out);
  
  for (unsigned int i=0;i<M.size();i++) {
    for (unsigned int j=0;j<M[i].size();j++) {
      file<<M[i][j]<<"\t";
    }
    if (br); file<<std::endl;
  }
  file.close();
}
// Imprime datos de los arreglos matricales
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
    if (br); file<<std::endl;
  }
  file.close();
}
//limpia los archivos
void file_wipe(std::string ARCHIVO){
  std::fstream file;
  file.open (ARCHIVO.c_str(), std::fstream::out);
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