#include <iostream>
#include <cmath>
#include <armadillo>
#include <cstdio>
#include <random>
#include <ctime>
#include <string>
#include <typeinfo>

#include "MainClass.hpp"

using namespace std;
using namespace arma;

//used to solve problem c

int main(int argc, char const *argv[])
{
  //from command line
  int N = atoi(argv[1]);
  double exp_MC_cyc = atof(argv[2]);

  double b = 100000; //temperatures we are interested in
  double cLL = 0.1; //concentration times L*L of rods.

  MainClass test(N, "test", pow(10, exp_MC_cyc-2)); //initialize random start

  //run for b = 1/T
  test.initialize_random();
  test.equilibrate(b, cLL, pow(10, exp_MC_cyc-2));
  test.Run(b, cLL, pow(10, exp_MC_cyc));
  test.reset();


  return 0;
}
