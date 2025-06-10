#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double eqa(double u_old, double u_x_prec_old, double u_x_next_old, double alpha, double pas_x, double pas_t){
    return u_old + alpha*((u_x_prec_old + u_x_next_old - 2*u_old)/(pas_x*pas_x))*pas_t;
}

double eqb(double u_x_prec, double tk, double q0, double A, double T, double pas_t, double pas_x, double k){
    return u_x_prec + (q0 + A*sin(M_PI*2/T * tk*pas_t))*pas_x/k;
}

int main(){
	double rho = 1870; // kg /m3
	double C = 1076; // J /kg /°C
	double ua = 1; // °C
    double pas_t, pas_x, q0, A, alpha, T, k, L;
    L = 4;
    int n = 200;
    k = 1.72;
    alpha = k / (rho * C);
    A =  1.6e6;
    double A_si = 1.6e6/3600/24;
    q0 = 3.5e4;
    double q0_si = 3.5e4/3600/24;

    T = 3.1536e7;
    pas_x = L/(n+1);
    pas_t = 0.1*pow(pas_x,2)/alpha;
    //printf("%lf,%lf,%lf,%lf,%lf,%lf\n",k,alpha,q0,T,pas_x,pas_t);

    double u_old[n+2];
    double u[n+2];
    
    double K=sqrt(M_PI/(alpha*T))*L;
    double P=(1-exp(2*K))*sin(K) - (1+exp(2*K))*cos(K); // Eq 33
    double Q=(1-exp(2*K))*sin(K) + (1+exp(2*K))*cos(K); // Eq 33
    double phi = atan(-Q/P); // Eq 36
    double b = A_si*sqrt(alpha*T/M_PI)/(k*(Q*sin(phi)-P*cos(phi))); // Eq 34
    double a = -b*exp(2*K);  // Eq 24
    double c = -q0_si/k; // Eq 35
    double d = ua+q0_si*L/k; // Eq 36
    printf("%lf,%lf,%lf,%lf,%lf\n",a,b,c,d,phi);

    //init
    for(int i=0; i<n+2; i++){
        u[i] = ua;
        u_old[i] = ua;
    }
    u_old[n+1] = ua + q0_si*pas_x/k;

    int Nt = 5.5 * T;
    for(int i=1; i<Nt+1; i++){
        for(int j=1; j<n+1; j++){
            u[j] = eqa(u_old[j],u_old[j-1],u_old[j+1],alpha,pas_x,pas_t);
                    //printf("%lf", u[j]);

        }
        u[n+1] = eqb(u[n],i,q0_si,A_si,T,pas_t,pas_x,k);
        for(int j=0; j<n+2; j++){
            u_old[j] = u[j];
        }
    }
    for(int ind=0; ind<n+2; ind++){
        printf("%lf\n",u[ind]);
    }
}