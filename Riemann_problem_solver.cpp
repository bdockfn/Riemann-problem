#include "Riemann_problem_solver.h"

double determine_v_star_function(double p, double rho_L, double v_L, double p_L, double rho_R, double v_R, double p_R,
            double c_L, double c_R, double v_A, double v_B, double v_C, string configuration, double gamma_const) {
   double v_star_L, v_star_R;

   if (configuration == "left rarefaction, right shock") {
      v_star_L = v_L - 2.0 *c_L /(gamma_const - 1.0) *(pow((p /p_L), (gamma_const - 1.0) /(2.0 *gamma_const)) - 1.0); //left rarefaction wave
      v_star_R = v_R + 1.0 /(rho_R *c_R) *(p - p_R) /sqrt((gamma_const + 1.0) /(2.0 *gamma_const) *p /p_R + 
               (gamma_const - 1.0) /(2.0 *gamma_const));                                                              //right shock wave
   } else if (configuration == "both shock") {
      v_star_L = v_L - 1.0 /(rho_L *c_L) *(p - p_L) /sqrt((gamma_const + 1) /(2 *gamma_const) *p /p_L + 
               (gamma_const - 1.0) /(2.0 *gamma_const));                                                              //left shock wave
      v_star_R = v_R + 1.0 /(rho_R *c_R) *(p - p_R) /sqrt((gamma_const + 1.0) /(2.0 *gamma_const) *p /p_R + 
               (gamma_const - 1.0) /(2.0 *gamma_const));                                                              //right shock wave
   } else if (configuration == "both rarefaction") {
      v_star_L = v_L - 2.0 *c_L /(gamma_const - 1.0) *(pow(p /p_L, (gamma_const - 1.0) /(2.0 *gamma_const)) - 1.0);   //left rarefaction wave
      v_star_R = v_R + 2.0 *c_R /(gamma_const - 1.0) *(pow(p /p_R, (gamma_const - 1.0) /(2.0 *gamma_const)) - 1.0);   //right rarefaction wave
   } else if (configuration == "left shock, right rarefaction") {
      v_star_L = v_L - 1.0 /(rho_L *c_L) *(p - p_L) /sqrt((gamma_const + 1) /(2 *gamma_const) *p /p_L + 
               (gamma_const - 1.0) /(2.0 *gamma_const));                                                              //left shock wave
      v_star_R = v_R + 2.0 *c_R /(gamma_const - 1.0) *(pow(p /p_R, (gamma_const - 1.0) /(2.0 *gamma_const)) - 1.0);   //right rarefaction wave
   }

   return v_star_L - v_star_R;
}

double determine_v_star_function_derivative(double p, double rho_L, double v_L, double p_L, double rho_R, double v_R, double p_R,
            double c_L, double c_R, double v_A, double v_B, double v_C, string configuration, double gamma_const) {
   double v_star_L_derivative, v_star_R_derivative;

   if (configuration == "left rarefaction, right shock") {
      v_star_L_derivative = -2.0 *c_L /(2.0 *gamma_const) *
                           (pow(p /p_L, -1.0 *(gamma_const + 1.0) /(2.0 *gamma_const)));                                       //left rarefaction wave
      v_star_R_derivative = 1.0 /(rho_R *c_R) *1 /sqrt((gamma_const + 1) /(2.0 *gamma_const) *p /p_R + 
                           (gamma_const - 1.0) /(2.0 *gamma_const)) +
                           1.0 /(rho_R *c_R) *(p - p_R) *1.0 /(2.0 *sqrt((gamma_const + 1.0) /(2.0 *gamma_const) *p /p_R +
                           (gamma_const - 1.0) /(2.0 *gamma_const))) *((gamma_const + 1.0) /(2.0 *gamma_const) /p_R);          //right shock wave
   } else if (configuration == "both shock") {
      v_star_L_derivative = -1.0 /(rho_L *c_L) *1.0 /sqrt((gamma_const + 1.0) /(2.0 *gamma_const) *p /p_L + 
                           (gamma_const - 1.0) /(2.0 *gamma_const)) - 1.0 /(rho_L *c_L) *(p - p_L) *1.0 /(2.0 *sqrt((gamma_const + 1.0) 
                           /(2.0 *gamma_const) *p /p_L + (gamma_const - 1.0) /(2.0 *gamma_const))) 
                           *((gamma_const + 1.0) /(2.0 *gamma_const) /p_L);                                                    //left shock wave
      v_star_R_derivative = 1.0 /(rho_R *c_R) *1.0 /sqrt((gamma_const + 1) /(2.0 *gamma_const) *p /p_R +
                           (gamma_const - 1.0) /(2.0 *gamma_const)) + 1.0 /(rho_R *c_R) *(p - p_R) *1.0 /(2.0 *sqrt((gamma_const + 1.0) 
                           /(2.0 *gamma_const) *p /p_R + (gamma_const - 1.0) /(2.0 *gamma_const))) *((gamma_const + 1.0) 
                           /(2.0 *gamma_const) /p_R);                                                                          //right shock wave
   } else if (configuration == "both rarefaction") {
      v_star_L_derivative = -2.0 *c_L /(2.0 *gamma_const) *(pow(p /p_L, -1.0 *(gamma_const + 1.0) /(2.0 *gamma_const)));       //left rarefaction wave
      v_star_R_derivative = 2.0 *c_R /(2.0 *gamma_const) *(pow(p /p_R, -1.0 *(gamma_const + 1.0) /(2.0 *gamma_const)));        //right rarefaction wave
   }
   else if (configuration == "left shock, right rarefaction") {
      v_star_L_derivative = -1.0 /(rho_L *c_L) *1.0 /sqrt((gamma_const + 1.0) /(2.0 *gamma_const) *p /p_L + (gamma_const - 1.0) 
      /(2.0 *gamma_const)) - 1.0 /(rho_L *c_L) *(p - p_L) *1.0 /(2.0 *sqrt((gamma_const + 1.0) /(2.0 *gamma_const) *p /p_L +
      (gamma_const - 1.0) /(2.0 *gamma_const))) *((gamma_const + 1.0) /(2.0 *gamma_const) /p_L);                               //left shock wave
      v_star_R_derivative = 2.0 *c_R /(2.0 *gamma_const) *(pow(p /p_R, -1.0 *(gamma_const + 1.0) /(2.0 *gamma_const)));        //right rarefaction wave
   }

   return v_star_L_derivative - v_star_R_derivative;
}

double Newthon_method_for_Riemann_problem(double rho_L, double v_L, double p_L, double rho_R, double v_R, double p_R,
            double c_L, double c_R, double v_A, double v_B, double v_C, string configuration, double gamma_const){
   double p_0, p_n, eps;
   int i;
   double func, func_derivative;

   eps = 1e-8;
   p_0 = 1.0 /2.0 *(p_L + p_R);
   p_n = p_0 + 1.0;

      i = 0;
   while (abs(p_n - p_0) > eps | i == 0) {
      p_0 = p_n;
      func = determine_v_star_function(p_0, rho_L, v_L, p_L, rho_R, v_R, p_R, c_L, c_R, v_A, v_B, v_C, configuration, gamma_const);
      func_derivative = determine_v_star_function_derivative(p_0, rho_L, v_L, p_L, rho_R, v_R, p_R, c_L, c_R, v_A, v_B, v_C, configuration, gamma_const);
      p_n = p_0 - func /func_derivative;
      i = i + 1;
   }

   return p_n;
}

void Riemann_problem_solver(string filename, int n, double moment_of_time, string outputname){

    // reading input data from filename
    std::ifstream fin(filename);
    double  rho_L, v_L, p_L, rho_R, v_R, p_R;
    fin >> rho_L >> v_L >> p_L >> rho_R >> v_R >> p_R;
    fin.close();

    // determine config
    double gamma_const = 5.0 /3;
    string configuration;
    double  c_L = sqrt(gamma_const *p_L /rho_L);
    double  c_R = sqrt(gamma_const *p_R /rho_R);
    double v_C = -2.0 /(gamma_const - 1.0) *(c_L + c_R);
    double v_A = 2.0 *c_L /(gamma_const - 1.0) *(pow((p_R /p_L), (gamma_const - 1.0)/(2.0 *gamma_const)) - 1.0);
    double v_B = 1 /(rho_R *c_R) *(p_L - p_R) /
                sqrt((gamma_const + 1.0) /(2.0 *gamma_const) *p_L /p_R + (gamma_const - 1.0) /(2.0 *gamma_const));
    if (v_L - v_R < v_C) {
        configuration = "vacuum?";
    } else if ((v_L - v_R) < v_B&  (v_L - v_R) > v_A) {
       configuration = "left rarefaction, right shock";
    } else if ((v_L - v_R) < v_A&  (v_L - v_R) > v_B) {
       configuration = "left shock, right rarefaction";
    } else if ((v_L - v_R) > v_B) {
       configuration = "both shock";
    } else if ((v_L - v_R) < v_A) {
       configuration = "both rarefaction";
    }
    std::cout << configuration << "\n";

    // declare variables for data vectors
    std::vector<double> x, p, v, rho;
    x.resize(n);
    p.resize(n);
    v.resize(n);
    rho.resize(n);

    // calculating a, b, vector x
    double a, b;


    if (configuration == "left rarefaction, right shock" | configuration == "left shock, right rarefaction") {
       b = 0.5;
       a = -0.5;
    } else {
       b = 0.3;
       a = -0.3;
    }
    x[0] = a;

    int i;
    for (i = 0; i < x.size(); i++) {
       x[i + 1] = x[i] + (b - a) /(x.size() - 1);
    }

    // calculating star region parameters
    double v_star, p_star, rho_star_L, rho_star_R;

    p_star = Newthon_method_for_Riemann_problem(rho_L, v_L, p_L, rho_R, v_R, p_R, c_L, c_R, v_A, v_B, v_C, configuration, gamma_const);

    if (configuration == "left rarefaction, right shock") {
       v_star = v_L - 2.0 *c_L /(gamma_const - 1.0) *(pow(p_star /p_L, (gamma_const - 1.0) /(2.0 *gamma_const)) - 1.0);
       rho_star_L = rho_L *pow(p_star /p_L, 1.0 /gamma_const);
       rho_star_R = rho_R *((gamma_const + 1.0) *p_star /p_R + (gamma_const - 1.0)) /
                   ((gamma_const + 1.0) + (gamma_const - 1.0) *p_star /p_R);
    }
    if (configuration == "both shock") {
       v_star = v_L - 1.0 /(rho_L *c_L) *(p_star - p_L) /
                         sqrt((gamma_const + 1.0) /(2.0 *gamma_const) *p_star /p_L + (gamma_const - 1.0) /(2.0 *gamma_const));
       rho_star_R = rho_R *((gamma_const + 1.0) *p_star /p_R + gamma_const - 1.0) /((gamma_const + 1.0) + (gamma_const - 1.0) *p_star /p_R);
       rho_star_L = rho_L *((gamma_const + 1.0) *p_star /p_L + gamma_const - 1.0) /((gamma_const + 1.0) + (gamma_const - 1.0) *p_star /p_L);
    }
    if (configuration == "both rarefaction") {
       v_star = v_L - 2.0 *c_L /(gamma_const - 1.0) *(pow(p_star /p_L, (gamma_const - 1.0) /(2.0 *gamma_const)) - 1.0);
       rho_star_L = rho_L *pow(p_star /p_L, 1.0 /gamma_const);
       rho_star_R = rho_R *pow(p_star /p_R, 1.0 /gamma_const);
    }
    if (configuration == "left shock, right rarefaction") {
       rho_star_R = rho_R *pow(p_star /p_R, 1.0 /gamma_const);
       rho_star_L = rho_L *((gamma_const + 1.0) *p_star /p_L + gamma_const - 1.0) /
                   ((gamma_const + 1.0) + (gamma_const - 1.0) *p_star /p_L);
       v_star = v_R + 2.0 *c_R /(gamma_const - 1.0) *(pow(p_star /p_R, (gamma_const - 1.0) /(2.0 *gamma_const)) - 1.0);
    }
    
    double x_break;
    x_break = v_star *moment_of_time;

    // calculating Head and Tail parameters
 
    
    double v_RWHead, v_RWTail, v_RW, rho_RW, D_R, D_L, p_RW;
    double v_RWHead_R, v_RWHead_L;
    double v_RWTail_R, v_RWTail_L;

    //Newthon_method_for_Riemann_problem(data_sample, extra_data, configuration);

   if (configuration == "left rarefaction, right shock") {
      //std::cout << "configuration == 'left rarefaction, right shock'" << "\n";

      v_RWHead = v_L - c_L;                                                                                                           //left rarefaction wave
      v_RWTail = v_star - c_L *pow((p_star /p_L), (gamma_const - 1) /(2 *gamma_const));

      D_R = v_R + c_R *pow((gamma_const + 1) /(2 *gamma_const) *p_star /p_R + (gamma_const - 1) /(2 *gamma_const), 1.0 /2.0);         //right shock wave
      for (i = 0; i < x.size(); i++) {
         if (x[i] >= D_R *moment_of_time) {                                                                                           //behind right shock wave
            p[i] = p_R;
            v[i] = v_R;
            rho[i] = rho_R;
         } if (x[i] <= v_RWHead *moment_of_time) {                                                                          //inside left rarefaction wave
            p[i] = p_L;
            v[i] = v_L;
            rho[i] = rho_L;
         } if ((x[i] >= v_RWHead *moment_of_time)&  (x[i] <= v_RWTail *moment_of_time)){                                   //inside  left rarefaction wave
            v_RW = v_L *(gamma_const - 1) /(gamma_const + 1) + 2 /(gamma_const + 1) *(x[i] /moment_of_time + c_L);

            rho_RW = rho_L *pow((2 /(gamma_const + 1) + (gamma_const - 1) /(c_L *(gamma_const + 1)) *(v_L - x[i] /moment_of_time)), 2.0 /(gamma_const - 1));

            p_RW = p_L *pow((2 /(gamma_const + 1) + (gamma_const - 1) /(c_L *(gamma_const + 1)) *(v_L - x[i] /moment_of_time)), 2.0 *gamma_const /(gamma_const - 1));

            p[i] = p_RW;
            v[i] = v_RW;
            rho[i] = rho_RW;
         } if ((x[i] <= D_R *moment_of_time)&  (x[i] >= v_RWTail *moment_of_time)) {       //between  left rarefaction wave and right shock wave, star region
            p[i] = p_star;
            v[i] = v_star;
            if (x[i] < x_break) {
               rho[i] = rho_star_L;
            } else {
               rho[i] = rho_star_R;
            }
         }
      }
   } if (configuration == "left shock, right rarefaction") {
      //std::cout << "configuration == 'left shock, right rarefaction'" << "\n";
      D_L = v_L - c_L *pow((gamma_const + 1.0) /(2.0 *gamma_const) *p_star /p_L + (gamma_const - 1.0) /(2.0 *gamma_const), 1.0 /2.0);
      v_RWHead_R = v_R + c_R;
      v_RWTail_R = v_star + c_R *pow((p_star /p_R), (gamma_const - 1.0) /(2.0 *gamma_const));
      
      for (i = 0; i < x.size(); i++) {
         if (x[i] <= D_L *moment_of_time) {                                                                    //behind the left shock wave
            p[i] = p_L;
            v[i] = v_L;
            rho[i] = rho_L;
         } else if ((x[i] > D_L *moment_of_time)&  (x[i] < v_RWTail_R *moment_of_time)) {
            p[i] = p_star;
            v[i] = v_star;
            if (x[i] < x_break) {
               rho[i] = rho_star_L;
            } else {
               rho[i] = rho_star_R;
            }
         } else if ((x[i] <= v_RWHead_R *moment_of_time)&  (x[i] >= v_RWTail_R *moment_of_time)) {               //into right rarefaction wave
            v_RW = v_R *(gamma_const - 1.0) /(gamma_const + 1.0) + 2.0 /(gamma_const + 1.0) *(x[i] /moment_of_time - c_R);

            rho_RW = rho_R *pow((2.0 /(gamma_const + 1.0) - (gamma_const - 1.0) /(c_R *(gamma_const + 1.0)) *(v_R - x[i] /moment_of_time)), 2.0 /(gamma_const - 1.0));

            p_RW = p_R *pow((2.0 /(gamma_const + 1.0) - (gamma_const - 1) /(c_R *(gamma_const + 1.0)) *(v_R - x[i] /moment_of_time)), 2.0 *gamma_const /(gamma_const - 1.0));
            p[i] = p_RW;
            v[i] = v_RW;
            rho[i] = rho_RW;
         } else if (x[i] > v_RWHead_R *moment_of_time) {
            p[i] = p_R;
            rho[i] = rho_R;
            v[i] = v_R;
         }
      }
   }  if (configuration == "both shock") {
      //std::cout << "configuration == 'both shock'" << "\n";
 
      D_R = v_R + c_R *pow((gamma_const + 1.0) /(2.0 *gamma_const) *p_star /p_R + (gamma_const - 1.0) /(2.0 *gamma_const), 1.0 /2.0);
      D_L = v_L - c_L *pow((gamma_const + 1.0) /(2.0 *gamma_const) *p_star /p_L + (gamma_const - 1.0) /(2.0 *gamma_const), 1.0 /2.0);
      for (i = 0; i < x.size(); i++) {
         if (x[i] >= D_R *moment_of_time) {                                                             //behind the right shock wave
            p[i] = p_R;
            v[i] = v_R;
            rho[i] = rho_R;
         } else if (x[i] <= D_L *moment_of_time) {                                                      //behind the left shock wave
            p[i] = p_L;
            v[i] = v_L;
            rho[i] = rho_L;
         } else if ((x[i] <= D_R *moment_of_time)&  (x[i] >= D_L *moment_of_time)) {                    //between left and right shock wave, star region
            p[i] = p_star;
            v[i] = v_star;
            if (x[i] < x_break) {
               rho[i] = rho_star_L;
            } else {
               rho[i] = rho_star_R;
            }
         }
      }
   } if (configuration == "both rarefaction") {
      //std::cout << "configuration == 'both rarefaction'" << "\n";
 
      v_RWHead_R = v_R + c_R;                                                                           //right rarefaction wave
      v_RWTail_R = v_star + c_R *pow((p_star /p_R), (gamma_const - 1.0) /(2.0 *gamma_const));
 
      v_RWHead_L = v_L - c_L;                                                                           //left rarefaction wave
      v_RWTail_L = v_star - c_L *pow((p_star /p_L), (gamma_const - 1.0) /(2.0 *gamma_const));
      for (i = 0; i < x.size(); i++) {
         if (x[i] < v_RWHead_L *moment_of_time) {
            p[i] = p_L;
            v[i] = v_L;
            rho[i] = rho_L;
         } if (x[i] > v_RWHead_R *moment_of_time) {
            p[i] = p_R;
            v[i] = v_R;
            rho[i] = rho_R;
         } if ((x[i] <= v_RWHead_R *moment_of_time)&  (x[i] >= v_RWTail_R *moment_of_time)) {                              //into right rarefaction wave
            v_RW = v_R *(gamma_const - 1.0) /(gamma_const + 1.0) + 2.0 /(gamma_const + 1.0) *(x[i] /moment_of_time - c_R);
            rho_RW = rho_R *pow((2.0 /(gamma_const + 1.0) - (gamma_const - 1.0) /(c_R *(gamma_const + 1.0)) *(v_R - x[i] /moment_of_time)), 2.0 /(gamma_const - 1.0));
            p_RW = p_R *pow((2.0 /(gamma_const + 1.0) - (gamma_const - 1) /(c_R *(gamma_const + 1.0)) *(v_R - x[i] /moment_of_time)), 2.0 *gamma_const /(gamma_const - 1.0));
            p[i] = p_RW;
            v[i] = v_RW;
            rho[i] = rho_RW;
         } if ((x[i] >= v_RWHead_L *moment_of_time)&  (x[i] <= v_RWTail_L *moment_of_time)) {                              //into left rarefaction wave
            v_RW = v_L *(gamma_const - 1.0) /(gamma_const + 1.0) + 2.0 /(gamma_const + 1.0) *(x[i] /moment_of_time + c_L);
            rho_RW = rho_L *pow((2.0 /(gamma_const + 1) + (gamma_const - 1.0) /(c_L *(gamma_const + 1.0)) *(v_L - x[i] /moment_of_time)), 2.0 /(gamma_const - 1.0));
            p_RW = p_L *pow((2.0 /(gamma_const + 1.0) + (gamma_const - 1.0) /(c_L *(gamma_const + 1.0)) *(v_L - x[i] /moment_of_time)), 2.0 *gamma_const /(gamma_const - 1.0));
            p[i] = p_RW;
            v[i] = v_RW;
            rho[i] = rho_RW;
         } if ((x[i] >= v_RWTail_L *moment_of_time)&  (x[i] <= v_RWTail_R *moment_of_time)) {
            p[i] = p_star;
            v[i] = v_star;
            if (x[i] < x_break) {
               rho[i] = rho_star_L;
            } else {
               rho[i] = rho_star_R;
            }
         }
      }
   }

   // writing results to output file
   std::ofstream fout(outputname);
   fout << "data: x, p, v, rho"<< "\n";
   fout << configuration<< "\n";
   for (int i = 0; i < n; i++) {
      fout << x[i] << " " << p[i] << " " << v[i] << " " << rho[i] << "\n";
   }
   fout.close();  
   std::cout << "output: x, p, v, rho"<< "\n"<< "\n";
}