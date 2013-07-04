#include <cstdlib>
#include "experiment.h"

using namespace std;

int main(int argc, char **argv) {
  
  //Entry parameters
  unsigned int L=	atoi(argv[1]);
  unsigned int numexps=	atoi(argv[2]);
  vector<double> rho = 	str2vec(argv[3]);
  vector<double> Temp=	str2vec(argv[4]);
  vector<double> campos=str2vec(argv[5]);
  vector<double> tau =	str2vec(argv[6]);
  
  Gen_exp(L,numexps,rho,Temp,campos,tau,"Liufix");
  
  return 0;
}

