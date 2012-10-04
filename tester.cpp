#include "tester.h"

using namespace std;

tester::tester(unsigned int L, double eps):relaxor(L)
{
  errtol = eps;
  relaxor.init(0,"test",false);
  runAllTests();  
}

void tester::runAllTests()
{
  cout<<"Ejecutando todas la pruebas"<<endl;
  sizes();
  rw_sigma();
  test_deltaH();
}

void tester::sizes(){
  clock_t cl_start = clock();
  cout<<"Verificar longitud de arreglos: ";
  unsigned int PNR = relaxor.PNR;
  assert(relaxor.sigma.size() == PNR);
  assert(relaxor.mu_E.size() == PNR);
  assert(relaxor.G.size() == PNR);
  assert(relaxor.J.size() == PNR);
 cout<<(double) (clock()-cl_start)/CLOCKS_PER_SEC<<"s\n";
}

void tester::test_deltaH()
{
  clock_t cl_start = clock();
  cout<<"Calcular deltaH: ";
  for(unsigned int idflip = 0; idflip < relaxor.PNR; idflip++){
    double H0 = relaxor.total_E(0);
    double dH = relaxor.delta_E(idflip,0);
    relaxor.sigma[idflip] *= -1;
    double H1 = relaxor.total_E(0);
    assert(abs(H1-H0-dH)<errtol);    
  }
 cout<<(double) (clock()-cl_start)/CLOCKS_PER_SEC<<"s\n";
}

void tester::rw_sigma()
{
  clock_t cl_start = clock();
  
  cout<<"Almacenar arreglos de sigma tamaño "<<relaxor.PNR<<": ";
  string savefile = "sigmasave.dat";
  std::vector<int> oldsigma = relaxor.sigma;
  array_print_bin(oldsigma,savefile);
  
  struct stat file;
  if (stat(savefile.c_str(), &file) == -1)
    cerr<<"no hay archivo";
  assert (file.st_size == relaxor.PNR*sizeof(int));
  std::ifstream ifile(savefile.c_str());
  
  int * newsigma = new int[relaxor.PNR];
  ifile.read((char *)&newsigma[0],relaxor.PNR*sizeof(int));
  for(unsigned int s=0;s<relaxor.PNR;s++)
    assert(newsigma[s]==oldsigma[s]);

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
