#include <math.h>
#include <stdio.h>

// Ce code correspond au problème P3 dans le papier de Xu 2022 :
// Analytical solutions for heat conduction problems with 3 kinds of periodic
// boundary conditions and their application
// Applied Mathematics and Computation

// Ce code calcule le profil de température analytique d'un slab avec
// une condition de flux en sinusoïde temporelle en haut (Neumann, Eq 3)
// et température constante en bas (Dirichlet).


void compute_constants(double T,double alpha, double q0,double A,double ua,double x0, double k,
                       double* p_a, double* p_b, double* p_c, double* p_d, double* p_phi)
{
  double K=sqrt(M_PI/(alpha*T))*x0;
  double P=(1-exp(2*K))*sin(K) - (1+exp(2*K))*cos(K); // Eq 33
  double Q=(1-exp(2*K))*sin(K) + (1+exp(2*K))*cos(K); // Eq 33
  double phi = atan(-Q/P); // Eq 36
  double b = A*sqrt(alpha*T/M_PI)/(k*(Q*sin(phi)-P*cos(phi))); // Eq 34
  double a = -b*exp(2*K);  // Eq 24
  double c = -q0/k; // Eq 35
  double d = ua+q0*x0/k; // Eq 36
  printf("%lf,%lf\n",c,d);
  *p_a = a; *p_b = b; *p_c = c; *p_d = d; *p_phi = phi;
}

double temperature(double x, double t, double T, double alpha,
                   double q0, double A, double ua, double x0,
                   double k)
{
  double a,b,c,d,phi;
  double x_adim=x/x0;
  double t_adim=t/T;
  double K=sqrt(M_PI/(alpha*T))*x0;
  compute_constants(T, alpha, q0, A, ua, x0, k, &a,&b,&c,&d,&phi);

  return a*exp(-K*x/x0)*cos(2*M_PI*t/T - K*x/x0 + K + phi)
        +b*exp( K*x/x0)*cos(2*M_PI*t/T + K*x/x0 - K + phi)
        + c*x + d; // Eq 16
}

int main(void)
{
  // Table 1 - physical parameters
  double k  = 1.72; // W /m /°C
  double rho = 1870; // kg /m3
  double C = 1076; // J /kg /°C
  double alpha =  k / (rho * C); // thermal diffusivity  
  // Table 4 - boundary conditions
  double q0 = 3.5e4; // J /d /m2  -> in days
  double q0_si = q0 / 3600 / 24; // J /s /m2 -> in seconds
  double A = 1.6e6;  // J /d /m2 -> in days
  double A_si = A / 3600 / 24;   // J /s /m2 -> in seconds
  double T = 3.1536e7; // s
  double x0 = 4; // m
  double ua = 1; // °C
  // Table 5 - probe position
  double x_obs= 0;
  double t_obs= 5*T;

  // Check results of Table 5 can be reproduced
  double a,b,c,d,phi;
  double u;
  printf("Result a,b,c,d,phi: %f\t%f\t%f\t%f\t%f\n", a,b,c,d,phi);
  printf("Expect a,b,c,d,phi: %f\t%f\t%f\t%f\t%f\n", -23.709,1.545,-0.236,1.942, -0.552);
  u = temperature(1, t_obs, T, alpha, q0_si, A_si, ua, x0, k);
  printf("Result u: %f\n", u);
  printf("Expect u: %f\n", 1.0);
}
