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
  double c;
  double L;
  double D;
  double cLLD;
  double dx;
  //save as double to avoid overflow

  //for writing to file
  string filename;
  int data_lines;

  //vectors
  Mat<double> abssingamma;
  Col<double> x;
  vec expectationvalues;
  vec finalvalues;
  Col<double> f;
  Col<int> index_vector;

  //exchange
  double J1, J2, J3;
  double Qlen;

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
  void set_exchange(double J10, double J20, double J30);

  //calcualting initial energy and change in energy
  double calc_energy();
  double spin_wave();
  double J(double qx, double qy);
  double delta_energy(int chosen_i, int leftright, double s);
  double delta_energy_sw(int chosen_i, int leftright, double s);

  //run metropolis, two different ways, with and without counting accepts
  //metroplis is called in the Run function
  void Metropolis();
  void Run(double b0, double c0, double L0, double D0, int nr_cycles);

  //time and equilibrate data, as well as reseting expectation values
  void equilibrate(double b0, double c0, double L0, double D0, int nr_cycles);
  double normalization();
  void reset();

  //write data to file
  void write_occupations(ofstream& OutputFile);
};

#endif
