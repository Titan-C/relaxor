#include "sistema.h"
#include <assert.h>
#include <cstdlib>
#include <cmath>
#include <ctime>

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
  void print_sigma();
  
};

tester::tester(unsigned int L, double eps):relaxor(L)
{
  errtol = eps;
  relaxor.init();
  runAllTests();  
}

void tester::runAllTests()
{
  cout<<"Ejecutando todas la pruebas"<<endl;
  test_deltaH();  
}

void tester::test_deltaH()
{
  clock_t cl_start = clock();
  cout<<"Test: deltaH: ";
  for(unsigned int idflip = 0; idflip < relaxor.return_PNR(); idflip++){
    double H0 = relaxor.total_E(0);
    double dH = relaxor.delta_E(idflip,0);
    relaxor.flip_sigma(idflip);
    double H1 = relaxor.total_E(0);
    assert(abs(H1-H0-dH)<errtol);    
  }
 cout<<(double) (clock()-cl_start)/CLOCKS_PER_SEC<<"s\n";
}

void tester::print_sigma()
{
  cout<<relaxor.return_PNR()<<" ";
  for(unsigned int i=0; i<relaxor.return_PNR();i++)
    cout<<relaxor.ret_sig(i)<<" ";
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