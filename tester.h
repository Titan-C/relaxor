#ifndef TESTER_H
#define TESTER_H

#include "material.h"
#include <cassert>

class tester
{
private:
  double errtol;
  Material relaxor;
public:
  tester(unsigned int L, double eps);
  void runAllTests();
  void test_deltaH();
  void rw_sigma();
  void material();
  
};
void mean_sd_stats(const std::vector< std::vector< double > >& M, double& mean, double& sd);
#endif // TESTER_H