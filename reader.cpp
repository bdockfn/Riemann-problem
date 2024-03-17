#include "reader.h"


void read_input(string filename, double  &rho_L, double  &v_L,
                double  &p_L, double  &rho_R, double  &v_R, double  &p_R){
    ifstream fin(filename);
    fin >> rho_L >> v_L >> p_L >> rho_R >> v_R >> p_R;
    fin.close();
}

void write_output(string &outputname, int N, vector<double> &x, vector<double> &p, vector<double> &v, vector<double> &rho)
{
   ofstream fout(outputname);
   fout << "data: x, p, v, rho"
        << "\n";
   for (int i = 0; i < N; i++){
      fout << x[i] << " " << p[i] << " " << v[i] << " " << rho[i] << "\n";
   }
   fout.close();
}