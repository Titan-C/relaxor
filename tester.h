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
  void sizes();
  
};

#endif // TESTER_H