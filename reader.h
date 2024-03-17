#pragma once
using namespace std;
#include <fstream>
#include <iostream>
#include <string>
#include <vector>



// reading input data from filename
void read_input(string filename, double  &rho_L, double  &v_L,
                double  &p_L, double  &rho_R, double  &v_R, double  &p_R);

void write_output(string &outputname, int N, vector<double> &x, vector<double> &p, 
            vector<double> &v, vector<double> &rho);

