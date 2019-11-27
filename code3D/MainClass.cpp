#include "MainClass.hpp"



const double PI = 3.1415926535897932384626433832795028841971693993751;
const double SQRTTHREEOVERTWO = 0.8660254038;



MainClass::MainClass()
{}



MainClass::MainClass(int size, string save_name, int amount_of_data)
{
  //initalize everything which is decleared in .hpp
  N = size;
  E = 0;
  b = 100;
  dang = PI/float(N);

  f = Col<double>(N);
  abssingamma = Col<double>(N);
  cosdang = Col<double>(N);
  sindang = Col<double>(N);

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
    cosdang(i) = cos(dang*i);
    sindang(i) = sin(dang*i);
  }
}



//gives a random normalized initialization
void MainClass::initialize_random()
{
  double s = 0;
  for (int i = 0; i < N; i++)
  {
    f(i) = zero_to_one_distribution(generator);
    s += 2*PI*f(i)*sindang(i);
  }
  for (int i = 0; i < N; i++)
  {
    f(i) /= s*dang;
  }
}



//gives an uniform normalized initialization
void MainClass::initialize_uniform()
{
  double val = 1.0/double(4*PI);
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
    int1 += f(i)*log(4*PI*f(i))*sindang(i)*dang; //*dang?
  }

  double int2 = 0;
  for (int i= 0; i < N; i++)
  {
    for (int j = i+1; j < N; j++)
    {
      int2 += f(i)*f(j)*abssingamma(abs(i-j))*sindang(i)*sindang(j); //dang = delta phi
    }
  }
  int2 *= cLL*dang*dang; //0.5 forsvinner fordi vi ganger med 2.

  energy = 2*PI*int1 + 2*PI*2*PI*int2;

  return double(energy);
}



double MainClass::delta_energy(int chosen_i, int leftright, double s)
{
  //can I calculate delta_E as a function of s?
  int a, a1, a2;
  double K;
  double dE;
  double sk, sl;

  dE = 0;
  int indices [2] = {chosen_i, chosen_i + leftright};

  //single integral
  for (int k = 0; k < 2; k++)
  {
    sk = s/sindang(index_vector(indices[k]));
    a = (k == 0) ? 1 : -1;
    dE += 2*PI*(f(index_vector(indices[k])) + a*sk)*log(4*PI*(f(index_vector(indices[k])) + a*sk))*sindang(index_vector(indices[k]))*dang;
    dE -= 2*PI*f(index_vector(indices[k]))*log(4*PI*f(index_vector(indices[k])))*sindang(index_vector(indices[k]))*dang;
  }

  //cout << a*s << endl;

  //double integral
  for (int k = 0; k < 2; k++)
  {
    for (int l = 1; l < N+1; l++)
    {
      sk = s/sindang(index_vector(indices[k]));
      sl = s/sindang(index_vector(l));
      a1 = (k == 0) ? 1 : -1;
      if (index_vector(l) == index_vector(indices[0]) || index_vector(l) == index_vector(indices[1]))
      {
        a2 = (index_vector(l) == index_vector(indices[0])) ? 1 : -1;
        K = 0.5;
      }
      else
      {
        a2 = 0;
        K = 1;
      }
      dE += 2*PI*2*PI*K*cLL*dang*dang*((f(index_vector(indices[k])) + a1*sk)*(f(index_vector(l)) + a2*sl) - f(index_vector(indices[k]))*f(index_vector(l)))*abssingamma(abs(index_vector(indices[k])-index_vector(l)))*sindang(index_vector(indices[k]))*sindang(index_vector(l));
    }
  }
  return dE;
}



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

    //calculate new f
    r = zero_to_one_distribution(generator);
    s = (f(index_vector(chosen_i))*sindang(index_vector(chosen_i)) + f(index_vector(chosen_i + leftright))*sindang(index_vector(chosen_i + leftright)))*r - f(index_vector(chosen_i))*sindang(index_vector(chosen_i));

    //s = 0.01;

    delta_E = delta_energy(chosen_i, leftright, s);
    //cout << double(delta_E) << endl;

    //metropolis check if we want to accept step or not
    if (delta_E < 0 || zero_to_one_distribution(generator) <= exp(-b*delta_E))
    {
      //cout << double(delta_E) << endl;
      E += double(delta_E);     //update energy (now it is not "per site")
      //cout << E << endl;
      f(index_vector(chosen_i)) += s/sindang(index_vector(chosen_i));
      f(index_vector(chosen_i + leftright)) -= s/sindang(index_vector(chosen_i + leftright));   //update f
    }
  }
}

void MainClass::Run(double b0, double c0, double L0, int nr_cycles)
{
  b = b0;
  c = c0;
  L = L0;
  cLL = c*L*L;

  //start by calculating the energy and magnetization of configuration
  E = calc_energy();

  //cout << E << endl;

  //open data files
  ofstream outfile_occupations("../data3D.nosync/" +  filename + "_occupations.txt", std::ios_base::app);
  if (!outfile_occupations.is_open())
     cout<<"Could not open file" << endl;

  ofstream outfile_final("../data3D.nosync/" + filename + "_final.txt", std::ios_base::app);
  if (!outfile_final.is_open())
     cout<<"Could not open file" << endl;


  ofstream outfile_expect("../data3D.nosync/" + filename + "_expect.txt", std::ios_base::app);
  if (!outfile_expect.is_open())
    cout<<"Could not open file" << endl;

  ofstream outfile_norm("../data3D.nosync/" + filename + "_normalization.txt", std::ios_base::app);
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

    finalvalues += expectationvalues/nr_cycles;

    //write data to file

    if ( (cycle) % denominator == 0)
    {
      outfile_expect << c << " " << L << " " << b << " " << expectationvalues(0) << " " << expectationvalues(1) << endl;

      outfile_norm << normalization() << endl;

      //write_occupations(outfile_occupations);
    }


  }
  outfile_final << c << " " << L << " " << b << " " << finalvalues(0) << " " << finalvalues(1) << "\n";
  write_occupations(outfile_occupations);
}


//equilibrate by simply running a given number of MC cycles
void MainClass::equilibrate(double b0, double c0, double L0, int nr_cycles)
{
  b = b0;
  c = c0;
  L = L0;
  cLL = c*L*L;
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
    s += f(i)*sindang(i);
  }
  s *= 2*PI*dang;
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
