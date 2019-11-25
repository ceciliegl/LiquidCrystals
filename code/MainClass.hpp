#ifndef MAINCLASS
#define MAINCLASS

#include <iostream>
#include <cmath>
#include <armadillo>
#include <cstdio>
#include <string>
#include <iomanip>
#include <fstream>

using namespace std;
using namespace arma;

class MainClass
{
public:
  double E;
  int N;
  double b;
  double cLL;
  double dang;
  //save as double to avoid overflow

  //for writing to file
  string filename;
  int data_lines;

  //vectors
  Col<double> abssingamma;
  vec expectationvalues;
  vec finalvalues;
  Col<double> f;
  Col<int> index_vector;

  //for timing
  double start;
  double stop;
  double time_used;

  //random number generators
  std::mt19937 generator;
  std::uniform_real_distribution<double> zero_to_one_distribution;
  std::uniform_int_distribution<int> zero_to_L_distribution;

  MainClass();
  MainClass(int size, string save_name, int amount_of_data);

  //the function names are self explanetory
  void initialize_random();
  void initialize_uniform();

  //calcualting initial energy and magnetization
  double calc_energy();
  double delta_energy(int chosen_i, int leftright, double s);

  //run metropolis, two different ways, with and without counting accepts
  //metroplis is called in the Run function
  void Metropolis();
  void Run(double b0, double cLL0, int nr_cycles);

  //time and equilibrate data, as well as reseting expectation values
  void equilibrate(double b0, double cLL0, int nr_cycles);
  double normalization();
  void reset();

  //write data to file
  void write_occupations(ofstream& OutputFile);
};

#endif
