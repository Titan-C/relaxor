#include "tester.h"

using namespace std;

tester::tester(unsigned int L, double eps):relaxor(L,0,"test")
{
  errtol = eps;
  relaxor.set_interaction_dipole_config(false);
  runAllTests();  
}

void tester::runAllTests()
{
  cout<<"Ejecutando todas la pruebas"<<endl;
  material();
  rw_sigma();
  test_deltaH();
}

void tester::material(){
  clock_t cl_start = clock();
  cout<<"Test array sizes: "<<endl;
  unsigned int PNR = relaxor.PNR;
  assert(relaxor.sigma.size() == PNR);
  assert(relaxor.mu_E.size() == PNR);
  assert(relaxor.G.size() == PNR);
  assert(relaxor.J.size() == PNR);
  cout<<"Test interaction energy: "<<endl;
  double mean, sd;
  mean_sd_stats(relaxor.J,mean,sd);
  assert(abs(mean - relaxor.rho)<0.02);
  assert(abs(sd -1)<0.02);
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
  assert (file.st_size == int (relaxor.PNR*sizeof(int)));
  std::ifstream ifile(savefile.c_str());
  
  int * newsigma = new int[relaxor.PNR];
  ifile.read((char *)&newsigma[0],relaxor.PNR*sizeof(int));
  for(unsigned int s=0;s<relaxor.PNR;s++)
    assert(newsigma[s]==oldsigma[s]);

  std::remove(savefile.c_str());
  cout<<(double) (clock()-cl_start)/CLOCKS_PER_SEC<<"s\n";
}

//Calcular la desviación estandar de una matriz
void mean_sd_stats(const vector< std::vector< double > >& M, double& mean, double& sd)
{
  unsigned int entries;
  entries = M.size() * M[0].size();
  double * Aij = new double [entries];
  for(unsigned int i = 0 ; i<M.size(); i++){
    for(unsigned int j = 0; j<M[0].size(); j++)
      Aij[i*M[0].size() + j] = M[i][j];
  }
  mean = gsl_stats_mean(Aij,1,entries);
  sd = gsl_stats_sd_m (Aij, 1, entries,mean);
  delete[] Aij;
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
