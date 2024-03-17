#include "alg_functions.h"


void initial_conditions_grid(double a, double b, int N, double  rho_L, double  v_L,
                double  p_L, double  rho_R, double  v_R, double  p_R,
                vector<double> &p, vector<double> &v, vector<double> &rho, vector<double> &x){
    x[0] = a;
    for (int i = 1; i < N; i++){
    x[i] = x[0] + (b - a) / (N - 1)*i;
    }
    for (int i = 0; i < N; i++){
        if (x[i] >= 0){
           rho[i] = rho_R;
           v[i] = v_R;
            p[i] = p_R;
        }else{
           rho[i] = rho_L;
           v[i] = v_L;
           p[i] = p_L;
        }
   }
}

double adiabatic_sound_speed(double p, double rho, double gamma_const){
    return sqrt(gamma_const*p/rho);
}

double stability_condition_vmax(vector<double> rho, vector<double> v, vector<double> p, double gamma_const){
    int N = p.size();
    vector<double> v_csj(N);
    for (int i = 0; i < N; i++){
        v_csj[i] = abs(v[i]) + adiabatic_sound_speed(p[i], rho[i], gamma_const);
    }
    double v_max = *max_element(v_csj.begin(), v_csj.end());
    return v_max;
}

void vectors_u_F(vector<double> &rho, vector<double> &v, vector<double> &p, vector<vector<double>> &u,vector<vector<double>> &F, double gamma_const){
    int N = p.size();
    for (int i = 0; i < N; i++){
        double eps = p[i] / rho[i] / (gamma_const - 1.0);
        u[i][0] = rho[i];
        u[i][1] = rho[i] * v[i];
        u[i][2] = rho[i] * (eps + pow(v[i], 2) / 2.0);

        F[i][0] = rho[i] * v[i];
        F[i][1] = rho[i] * pow(v[i], 2) + p[i];
        F[i][2] = rho[i] * v[i] * (eps + pow(v[i], 2) / 2.0 + p[i] / rho[i]);
    }
}

void vect_u_j(int i, vector<double> &u_j, vector<double> &rho, vector<double> &v, vector<double> &p, double gamma_const){
    int N = p.size();
    double eps = p[i] / rho[i] / (gamma_const - 1.0);
    u_j[0] = rho[i];
    u_j[1] = rho[i] * v[i];
    u_j[2] = rho[i] * (eps + pow(v[i], 2) / 2.0);
}

void vect_F_j(int i, vector<double> &F_j, vector<double> &rho, vector<double> &v, vector<double> &p, double gamma_const){
    int N = p.size();
    double eps = p[i] / rho[i] / (gamma_const - 1.0);
    F_j[0] = rho[i] * v[i];
    F_j[1] = rho[i] * pow(v[i], 2) + p[i];
    F_j[2] = rho[i] * v[i] * (eps + pow(v[i], 2) / 2.0 + p[i] / rho[i]);
}

void vectors_FStar(vector<vector<double>> &FStar, vector<double> &rho, vector<double> &v, vector<double> &p, double gamma_const, double vmax){
    vector<double> F_left(3), F_right(3), u_left(3), u_right(3);
    int N = rho.size();
    vect_F_j(0, FStar[0], rho, v, p, gamma_const);
    vect_F_j(N-1, FStar[N], rho, v, p, gamma_const);

    for (int i = 1; i < N; i++){
        vect_F_j(i-1, F_left, rho, v, p, gamma_const);
        vect_F_j(i, F_right, rho, v, p, gamma_const);
        vect_u_j(i-1, u_left, rho, v, p, gamma_const);
        vect_u_j(i, u_right, rho, v, p, gamma_const);
            for (int l=0; l<3; l++){
               FStar[i][l] = (F_left[l]+F_right[l])/2.0 - abs(vmax)/2.0*(u_right[l]-u_left[l]);
            }
    }   
}

void new_time_layer(vector<vector<double>> &u, vector<vector<double>> FStar, double vmax, double Courant){
    int N = u.size();
    for (int i = 0; i < N-1; i++){
        for (int l=0; l<3; l++){
            u[i][l] = u[i][l] - Courant/vmax*( FStar[i+1][l]-FStar[i][l]);
        }
    }
        for (int l=0; l<3; l++){
        u[N-1][l] = u[N-1][l] - Courant/vmax* (FStar[N][l]-FStar[N-1][l]);
        }
}

void vectors_rho_v_p(vector<vector<double>> u, vector<double> &rho, vector<double> &v, vector<double> &p, double gamma_const){
    int N = p.size();
    for (int i = 0; i < N; i++){
        rho[i] = u[i][0];
        v[i] = u[i][1] / rho[i];
        p[i] = (u[i][2] / rho[i] - pow(v[i], 2) / 2.0) * rho[i] * (gamma_const - 1);
    }
}
