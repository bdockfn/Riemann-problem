#include "alg_functions.h"
#include "solver_algorithm.h"
#include "reader.h"
using namespace std;

int main() {
   string filename = "input/input1.txt";
   string outputname = "output/output1_1000_c1.05.txt";
   int N = 1000;
   double Courant = 1.05;
   double gamma_const = 5.0 / 3.0;
   double time_break = 0.10;
   double a = -1, b = 1;

   // reading input data
   double  rho_L, v_L, p_L, rho_R, v_R, p_R;
   vector<double> rho(N), p(N), v(N), x(N);
   read_input(filename, rho_L, v_L, p_L, rho_R, v_R, p_R);

   // algorithm
   scheme(a, b, N, rho_L, v_L, p_L, rho_R, v_R, p_R, p, v, rho, x, gamma_const, Courant, time_break);

   // writing results to output file
   write_output(outputname, N, x, p, v, rho);
}


