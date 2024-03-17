#pragma once
using namespace std;
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <vector>

double determine_v_star_function(double p, double rho_L, double v_L, double p_L, double rho_R, double v_R, double p_R,
            double c_L, double c_R, double v_A, double v_B, double v_C, string configuration, double gamma_const);

double determine_v_star_function_derivative(double p, double rho_L, double v_L, double p_L, double rho_R, double v_R, double p_R,
            double c_L, double c_R, double v_A, double v_B, double v_C, string configuration, double gamma_const);

double Newthon_method_for_Riemann_problem(double rho_L, double v_L, double p_L, double rho_R, double v_R, double p_R,
            double c_L, double c_R, double v_A, double v_B, double v_C, string configuration, double gamma_const);

void Riemann_problem_solver(string filename, int n, double moment_of_time, string outputname);