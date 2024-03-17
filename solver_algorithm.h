#include "alg_functions.h"

void scheme(double a, double b, int N, double rho_L, double v_L, double p_L, 
            double rho_R, double v_R, double p_R, vector<double> &p, 
            vector<double> &v, vector<double> &rho, vector<double> &x, 
            double gamma_const, double Courant, double time_break);