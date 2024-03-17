using namespace std;
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>


void initial_conditions_grid(double a, double b, int N, double  rho_L, double  v_L,
                             double  p_L, double  rho_R, double  v_R, double  p_R,
                             vector<double> &p, vector<double> &v, vector<double> &rho, 
                             vector<double> &x);

double adiabatic_sound_speed(double p, double rho, double gamma_const);

double stability_condition_vmax(vector<double> rho, vector<double> v, vector<double> p, 
                                double gamma_const);

void vectors_u_F(vector<double> &rho, vector<double> &v, vector<double> &p, 
                 vector<vector<double>> &u,vector<vector<double>> &F, double gamma_const);
                
void vect_u_j(int i,vector<double> &u_j, vector<double> &rho, vector<double> &v, 
              vector<double> &p, double gamma_const);

void vect_F_j(int i,vector<double> &F_j, vector<double> &rho, vector<double> &v, 
              vector<double> &p, double gamma_const);

void vectors_FStar(vector<vector<double>> &FStar, vector<double> &rho, vector<double> &v, 
                   vector<double> &p, double gamma_const, double vmax);

void new_time_layer(vector<vector<double>> &u, vector<vector<double>> FStar, double vmax, 
                    double Courant);
                    
void vectors_rho_v_p(vector<vector<double>> u, vector<double> &rho, vector<double> &v, 
                     vector<double> &p, double gamma_const);
