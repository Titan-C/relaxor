#include "sistema.h"
#include "impresor.h"
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>

using namespace std;

class tester
{
private:
  double errtol;
  Sistema relaxor;
public:
  tester(unsigned int L, double eps);
  void runAllTests();
  void test_deltaH();
  void rw_sigma();
  
};

tester::tester(unsigned int L, double eps):relaxor(L)
{
  errtol = eps;
  relaxor.init(0,false);
  runAllTests();  
}

void tester::runAllTests()
{
  cout<<"Ejecutando todas la pruebas"<<endl;
  rw_sigma();
  test_deltaH();
}

void tester::test_deltaH()
{
  clock_t cl_start = clock();
  cout<<"Calcular deltaH: ";
  for(unsigned int idflip = 0; idflip < relaxor.return_PNR(); idflip++){
    double H0 = relaxor.total_E(0);
    double dH = relaxor.delta_E(idflip,0);
    relaxor.flip_sigma(idflip);
    double H1 = relaxor.total_E(0);
    assert(abs(H1-H0-dH)<errtol);    
  }
 cout<<(double) (clock()-cl_start)/CLOCKS_PER_SEC<<"s\n";
}

void tester::rw_sigma()
{
  cout<<sizeof(int)<<""<<sizeof(double);
  clock_t cl_start = clock();
  cout<<"Almacenar arreglos de sigma tamaño "<<relaxor.return_PNR()<<": ";
  string savefile = "sigmasave.dat";
  int8_t * oldsigma = relaxor.ret_sigarr();
  unsigned int rep =2030;
  for(unsigned int i=0;i<rep;i++)
    array_print_bin(oldsigma,relaxor.return_PNR(),savefile,true);
  
  struct stat file;
  if (stat(savefile.c_str(), &file) == -1)
    cerr<<"no hay archivo";
  assert (file.st_size == rep*relaxor.return_PNR()*sizeof(int8_t));
  std::ifstream ifile(savefile.c_str());
  
  int8_t * newsigma = new int8_t[relaxor.return_PNR()];
  for(unsigned int i=0;i<rep;i++){
  ifile.read((char *)&newsigma[0],relaxor.return_PNR()*sizeof(int8_t));
  for(unsigned int s=0;s<relaxor.return_PNR();s++)
    assert(newsigma[s]==oldsigma[s]);
  }
  std::remove(savefile.c_str());
  cout<<(double) (clock()-cl_start)/CLOCKS_PER_SEC<<"s\n";
}

int main(int argc, char **argv) {
  //Parametros de entrada
  unsigned int L=16;
  double eps	=5e-12;
  
  if (argc==3){
    L=atoi(argv[1]);
    eps	=atof(argv[2]);
  }
  tester test(L,eps);
  
  return 0;
}