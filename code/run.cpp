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

  double b = 1e2; //temperatures we are interested in
  double c = 4.0; //concentration of rods.
  double L = 1.0; //length of rods.

  MainClass test(N, "test", pow(10, exp_MC_cyc-2)); //initialize random start

  test.initialize_uniform();
  //run for b = 1/T
  for (double A = 1; A <= 10; A += 1)
  {
    b = pow(10, A);
    test.equilibrate(b, c, L, pow(10, exp_MC_cyc-2));
    test.Run(b, c, L, pow(10, exp_MC_cyc));
    test.reset();
  }

  return 0;
}
