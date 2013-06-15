#ifndef IMPRESOR_H
#define IMPRESOR_H

#include <vector>
#include <fstream>
#include <sstream>


// Imprime datos de variables double
void out(double value, std::string ARCHIVO, bool br = true);
void array_print(const std::vector< double >& V, std::string ARCHIVO, unsigned int colsize, double scale);

// Imprime arreglos vectoriales de todo tipo de variable en texto y binario
template <typename T>
void array_print(const std::vector<T> & V, std::string ARCHIVO) {
  std::ofstream file(ARCHIVO.c_str(), std::ofstream::app);
  
  for(unsigned int i = 0; i<V.size(); i++)
    file<<V[i]<<"\t";
  file<<"\n";
  file.close();
}
template <typename T>
void array_print_bin(const std::vector< T >& V, std::string ARCHIVO) {
  std::ofstream file(ARCHIVO.c_str(), std::ofstream::app);
  file.write((char * )&V[0],V.size()*sizeof(T));
  file.close();
}
// Imprime arreglos matriciales de todo tipo de variable en texto y binario
template <typename T>
void array_print(const std::vector< std::vector<T> >& M, std::string ARCHIVO) {
  std::ofstream file(ARCHIVO.c_str());

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
  std::ofstream file(ARCHIVO.c_str(), std::ofstream::app);

  for (unsigned int i=0;i<M.size();i++)
    file.write((char * )&M[i][0],M[i].size()*sizeof(T));

  file.close();
}
#endif // IMPRESOR_H
