#include "MainClass.hpp"

const double PI = 3.1415926535897932384626433832795028841971693993751;

MainClass::MainClass()
{}

MainClass::MainClass(int size, string save_name, int amount_of_data)
{
  //initalize everything which is decleared in .hpp
  N = size;
  E = 0;
  b = 100;
  dang = 2*PI/float(N);

  f = Col<double>(N);
  abssingamma = Col<double>(N);

  expectationvalues = vec(2, fill::zeros);
  finalvalues = vec(2, fill::zeros);
  filename = save_name;
  data_lines = amount_of_data;

  std::mt19937 generator (std::clock());
  std::uniform_real_distribution<double> zero_to_one_distribution(0.0, 1.0);
  std::uniform_int_distribution<int> zero_to_L_distribution(0, N-1);

  //create the index vector for periodic boundary conditions
  index_vector = Col<int>(N+2);
  for (int i = 1; i < N+1; i++)
  {
    index_vector(i) = i-1;
  }
  index_vector(0) = N-1;
  index_vector(N+1) = 0;


  for (int i = 0; i < N; i++)
  {
    abssingamma(i) = abs(sin(dang*i));
  }
}

//gives a random normalized initialization
void MainClass::initialize_random()
{
  double s = 0;
  for (int i = 0; i < N; i++)
  {
    f(i) = zero_to_one_distribution(generator);
    s += f(i);
  }
  for (int i = 0; i < N; i++)
  {
    f(i) /= s*dang;
  }
}


//gives an uniform normalized initialization
void MainClass::initialize_uniform()
{
  double val = 1.0/double(2*PI);
  for (int i = 0; i < N; i++)
  {
    f(i) = val;
  }
}


//calculates total energy of system
double MainClass::calc_energy()
{
  double energy = 0;

  //calculate contributions from the integrals.

  double int1 = 0;
  for (int i = 0; i < N; i++)
  {
    int1 += f(i)*log(2*PI*f(i))*dang; //*dang?
  }

  double int2 = 0;
  for (int i= 0; i < N; i++)
  {
    for (int j = i+1; j < N; j++)
    {
      int2 += 2*f(i)*f(j)*abssingamma(abs(i-j)); //dang = delta phi
    }
  }
  int2 *= 0.5*cLL*dang*dang; //*dang*dang? I så fall må du endre normaliseringen av f.

  energy = int1 + int2;

  return double(energy);
}

double MainClass::delta_energy(int chosen_i, int leftright, double s)
{
  //can I calculate delta_E as a function of s?
  int a, a1, a2;
  double K;
  double dE;

  dE = 0;
  int indices [2] = {chosen_i, chosen_i + leftright};

  //single integral
  for (int k = 0; k < 2; k++)
  {
    a = (k == 0) ? 1 : -1;
    dE += (f(index_vector(indices[k])) + a*s)*log(2*PI*(f(index_vector(indices[k])) + a*s))*dang;
    dE -= f(index_vector(indices[k]))*log(2*PI*f(index_vector(indices[k])))*dang;
  }

  //double integral
  for (int k = 0; k < 2; k++)
  {
    for (int l = 1; l < N+1; l++)
    {
      if (l != k)
      {
        a1 = (k == 0) ? 1 : -1;
        if (l == indices[0] || l == indices[1])
        {
          a2 = (l == indices[0]) ? 1 : -1;
          K = 0.5;
        }
        else
        {
          a2 = 0;
          K = 1;
        }
        dE += K*cLL*dang*dang*((f(index_vector(indices[k])) + a1*s)*(f(index_vector(l)) + a2*s) - f(index_vector(indices[k]))*f(index_vector(l)))*abssingamma(abs(index_vector(indices[k])-index_vector(l)));
      }
    }
  }
  return dE;
}

//Can calculate the change in energy much simpler.


//the selection of a random f, and using the metropolis requirement for checking if we flip or not
void MainClass::Metropolis()
{
  double delta_E;
  int chosen_i;
  int leftright;

  double r;
  double s;

  for (int i = 0; i < N; i++)
  {
    //choose random site
    chosen_i = zero_to_L_distribution(generator)%N + 1; //we plus with one since we will use them on "index_vector", which is shifted by +1
    leftright = (zero_to_one_distribution(generator) > 0.5) ? 1 : -1;

    //beregn den nye okkupasjonen
    r = zero_to_one_distribution(generator);
    s = (f(index_vector(chosen_i)) + f(index_vector(chosen_i + leftright)))*r - f(index_vector(chosen_i));

    delta_E = delta_energy(chosen_i, leftright, s);
    //cout << double(delta_E) << endl;

    //metropolis check if we want to accept step or not
    if (delta_E < 0 || zero_to_one_distribution(generator) <= exp(-b*delta_E))
    {
      //cout << double(delta_E) << endl;
      E += double(delta_E);     //update energy (now it is not "per site")
      f(index_vector(chosen_i)) += s;
      f(index_vector(chosen_i + leftright)) -= s;   //update configuration
    }
  }
}

void MainClass::Run(double b0, double cLL0, int nr_cycles)
{
  b = b0;
  cLL = cLL0;

  //start by calculating the energy and magnetization of configuration
  E = calc_energy();

  //open data files
  ofstream outfile_occupations("../data.nosync/" +  filename + "_occupations.txt", std::ios_base::app);
  if (!outfile_occupations.is_open())
     cout<<"Could not open file" << endl;

  ofstream outfile_final("../data.nosync/" + filename + "_final.txt", std::ios_base::app);
  if (!outfile_final.is_open())
     cout<<"Could not open file" << endl;


  ofstream outfile_expect("../data.nosync/" + filename + "_expect.txt", std::ios_base::app);
  if (!outfile_expect.is_open())
    cout<<"Could not open file" << endl;

  ofstream outfile_norm("../data.nosync/" + filename + "_normalization.txt", std::ios_base::app);
  if (!outfile_norm.is_open())
    cout<<"Could not open file" << endl;

  //write f to file before we start, so we know initial configuration
  //write_occupations(outfile_occupations);

  int denominator = nr_cycles/data_lines;

  for (int cycle = 1; cycle <= nr_cycles; cycle++)
  {
    //call metropolis
    Metropolis();

    //update expectation values after metropolis
    expectationvalues(0) = E;
    expectationvalues(1) = E*E;

    //cout << E << endl;

    //write data to file
    if ( (cycle) % denominator == 0)
    {
      outfile_expect << b << " " << expectationvalues(0) << " " << expectationvalues(1) << endl;

      finalvalues += expectationvalues/nr_cycles;

      outfile_norm << normalization() << endl;

      //write_occupations(outfile_occupations);
    }

  }
  outfile_final << b << " " << finalvalues(0) << " " << finalvalues(1) << "\n";
  write_occupations(outfile_occupations);
}


//equilibrate by simply running a given number of MC cycles
void MainClass::equilibrate(double b0, double cLL0, int nr_cycles)
{
  b = b0;
  cLL = cLL0;
  for (int i = 0; i < nr_cycles/10; i++) //divide by 10
    {
      Metropolis();
    }
}

double MainClass::normalization()
{
  double s = 0;
  for (int i = 0; i < N; i++)
  {
    s += f(i);
  }
  s *= dang;
  return s;
}

void MainClass::reset()
{
  //Function used to reset the expectationvalues
  expectationvalues(0) = 0;
  expectationvalues(1) = 0;
  finalvalues(0) = 0;
  finalvalues(1) = 0;
}


void MainClass::write_occupations(ofstream& OutputFile)
{
  //Function used to write the current configuration to file.
  for (int i = 0; i < N; i++)
  {
    OutputFile << f(i) << " ";
  }
  OutputFile << " \n";
}
