#include "alg_functions.h"
#include "solver_algorithm.h"


void scheme(double a, double b, int N, double rho_L, double v_L, double p_L, double rho_R, 
            double v_R, double p_R, vector<double> &p, vector<double> &v, vector<double> &rho, 
            vector<double> &x, double gamma_const, double Courant, double time_break)
{
   // задать сетку
   initial_conditions_grid(a, b, N, rho_L, v_L, p_L, rho_R, v_R, p_R,
                           p, v, rho, x);
   double x_step = (b - a) / (N - 1);

   // условие устойчивости
   double vmax;
   vmax = stability_condition_vmax(rho, v, p, gamma_const);
   double t_step;
   t_step = Courant * x_step / vmax;
   double t = 0;

   vector<double> inside_vect(3);
   vector<vector<double>> u(N, inside_vect), F(N, inside_vect);

   vector<double> F_left(3), F_right(3), u_left(3), u_right(3);
   vector<vector<double>> FStar(N + 1, inside_vect);

   while (t < time_break){
      // алг. для одного вр. слоя
      vectors_u_F(rho, v, p, u, F, gamma_const);
      vectors_FStar(FStar, rho, v, p, gamma_const, vmax);
      new_time_layer(u, FStar, vmax, Courant);

      // распаковываем rho, v, p для нового слоя
      vectors_rho_v_p(u, rho, v, p, gamma_const);
      vmax = stability_condition_vmax(rho, v, p, gamma_const);

      t_step = Courant * x_step / vmax;
      t = t + t_step;
   }
}