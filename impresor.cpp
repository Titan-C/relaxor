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
//Cola de Gráficos
void plot_pol(std::string id_proc)
{
  std::ofstream file;
  file.open( "plots.p", std::ofstream::app);
  
  file<<"reset; \n set terminal png;\n";
  
  file<<"set output '"+id_proc+".png'\n";
  file<<"set xlabel \"Temperatura\'; \n set ylabel \"Polarizacion\";\n";
  file<<"set title \"Polarizacion media\";\n";

}
void plot_sus(std::string Exp_ID, double DJ, double rho,
	      const std::vector<double>& Temps, const std::vector<double>& Fields,
	      const std::vector<double>& tau){
  std::ofstream file;
  file.open( "plots.p");
  file<<"reset;\n set terminal png;\n";
  file<<"set ylabel\"$\\\\chi$\";\n";
  
  if (Exp_ID == "cool" || Exp_ID == "heat"){
    file<<"set xlabel \"$Temperatura [\\\\Delta J/k_B]$\";\n";
    /*Graficos de la susceptibilidad en función de la temperatura
     * para campo fijo y frecuencias disponibles*/
    for(unsigned int E=0; E<Fields.size() ; E++){
      std::ostringstream id_proc;
      id_proc<<"sus_"<<Exp_ID<<"_J"<<DJ<<"_p"<<rho<<"_E"<<Fields[E]/DJ;

      file<<"set output \'"<<id_proc.str()<<".png\'\n";
      file<<"set title \"Susceptibilidad ante Campos de amplitud constante $|E|="<<Fields[E]/DJ<<"[\\\\Delta J/\\\\overline{\\\\mu}]$ y distinta frecuencia\";\n";
      file<<"plot ";
      for(unsigned int t=0; t< tau.size(); t++){
	file<<"\""<<id_proc.str()<<"_t"<<tau[t]<<".dat\" w l title \"$\\\\tau= "<<tau[t];
	(t+1-tau.size()==0)? file<<"$\";\n\n" : file<<"$\", \\\n";
      }
    }
    /*Graficos de la susceptibilidad en función de la temperatura
     * para frecuencias fijas y campos disponibles*/
    for(unsigned int t=0;t<tau.size();t++){
      std::ostringstream id_proc;
      id_proc<<"sus_"<<Exp_ID<<"_J"<<DJ<<"_p"<<rho;
      file<<"set output \'"<<id_proc.str()<<"_t"<<tau[t]<<".png\'\n";
      file<<"set title \"Susceptibilidad ante Campos de frecuencia constante $\\\\tau="<<tau[t]<<"$ y amplitud variable\";\n";
      file<<"plot ";
      for(unsigned int E=0; E<Fields.size(); E++){
	file<<"\""<<id_proc.str()<<"_E"<<Fields[E]/DJ<<"_t"<<tau[t]<<".dat\" w l title \"$|E|= "<<Fields[E]/DJ;
	(E+1-Fields.size()==0)? file<<"$\";\n\n" : file<<"$\", \\\n";
      }
    }
  } else {
    file<<"set xlabel \"$Amplitud Campo [\\\\Delta J/\\\\overline{\\\\mu}]$\";\n";
    /*Graficos de la susceptibilidad en función de la Amplitud de campo
     * para temperatura fija y frecuencias disponibles*/
    for(unsigned int T=0; T<Temps.size() ; T++){
      std::ostringstream id_proc;
      id_proc<<"sus_"<<Exp_ID<<"_J"<<DJ<<"_p"<<rho<<"_T"<<Temps[T]/DJ;

      file<<"set output \'"<<id_proc.str()<<".png\'\n";
      file<<"set title \"Susceptibilidad a temperatura constante $T="<<Temps[T]/DJ<<"[\\\\Delta J/\\\\overline{\\\\mu}]$ y distinta frecuencia\";\n";
      file<<"plot ";
      for(unsigned int t=0; t< tau.size(); t++){
	file<<"\""<<id_proc.str()<<"_t"<<tau[t]<<".dat\" w l title \"$\\\\tau= "<<tau[t];
	(t+1-tau.size()==0)? file<<"$\";\n\n" : file<<"$\", \\\n";
      }
    }
    /*Graficos de la susceptibilidad en función de la Amplitud de campo
     * para frecuencias fijas y temperaturas disponibles*/
    for(unsigned int t=0;t<tau.size();t++){
      std::ostringstream id_proc;
      id_proc<<"sus_"<<Exp_ID<<"_J"<<DJ<<"_p"<<rho;
      file<<"set output \'"<<id_proc.str()<<"_t"<<tau[t]<<".png\'\n";
      file<<"set title \"Susceptibilidad ante Campos de frecuencia constante $\\\\tau="<<tau[t]<<"$ y distintas temperaturas\";\n";
      file<<"plot ";
      for(unsigned int T=0; T<Temps.size(); T++){
	file<<"\""<<id_proc.str()<<"_T"<<Temps[T]/DJ<<"_t"<<tau[t]<<".dat\" w l title \"$T= "<<Temps[T]/DJ;
	(T+1-Temps.size()==0)? file<<"$\";\n\n" : file<<"$\", \\\n";
      }
    }
  }
  file.close();
  std::system("gnuplot plots.p");
}