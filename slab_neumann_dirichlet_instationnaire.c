#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <rsys/cstr.h>	
#include <star/ssp.h>
#include <star/swf.h>

/*               (2,2,2)
 *      __________
 *     /         /|
 *    /    u1   / |
 *   /_________/  |   Temperature u1 sur le plan opposé (ox,oy)
 *   |         |  |   
 *   |     x   |  |      z
 *   |         |  /      |  y
 *   |         | /       | /
 *   |_________|/        |/____ x
 *(0,0,0)
 */
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
  *p_a = a; *p_b = b; *p_c = c; *p_d = d; *p_phi = phi;
}

double temperature(double * x, double t, double T, double alpha,
                   double q0, double A, double ua, double x0,
                   double k)
{
  double a,b,c,d,phi;
  double z = x[3];
  double z_adim=z/x0;
  double t_adim=t/T;
  double K=sqrt(M_PI/(alpha*T))*x0;
  compute_constants(T, alpha, q0, A, ua, x0, k, &a,&b,&c,&d,&phi);

  return a*exp(-K*z/x0)*cos(2*M_PI*t/T - K*z/x0 + K + phi)
        +b*exp( K*z/x0)*cos(2*M_PI*t/T + K*z/x0 - K + phi)
        + c*z + d; // Eq 16
}

struct swf_tabulation* Propagateur;

double fif(double *x_cur, double t_cur, double x0, FILE* fptr, double k, double q0, double T, double alpha, double eps);

// Température sur le plan inferieur
double u1(double * x){
	return 1;
}

// Flux sur le plan superieur
double q1(double t_cur, double k, double q0, double T){
	return -1/k*(q0 + sin(2*M_PI/T*t_cur));
}

double sum_G(double t, double tau, double R, double alpha){
	double seuil = 10e-10;
	double s = 0;
	int k = 1;
	double terme = 0;
	while(terme > seuil){
		terme = pow(-1,k)*pow(k,2)*exp(-alpha*pow(k*M_PI/R,2)*(t-tau));
		s += terme;
		k++;
	}
	printf("%lf\n",t);
	return s*M_PI/2*pow(R,4);
}

double wos(double *x_cur, double t_cur, double x0, FILE *fptr, 
	double k, double q0, double T, double alpha, double eps){

	//printf("%lf,%lf,%lf\n",x_cur[0],x_cur[1],x_cur[2]);

	//Echantillonnage
	double u = (double) rand() / (double) RAND_MAX;
	double phi = u*2.0*M_PI;

	double v = (double) rand() / (double) RAND_MAX;
	double theta = acos(-2.0*v+1.0);

	double w = (double) rand() / (double) RAND_MAX;
	double tirage_propagateur = swf_tabulation_inverse(Propagateur, SWF_QUADRATIC, w);
	
	// Calcul du rayon de la plus grande sphère inscrite dans le cube
    double radius = fmin(x_cur[2],x0-x_cur[2]);
    if(radius>x0/2){
		printf("Erreur rayon");
		return 1;
	}

	double tau = radius*radius*tirage_propagateur/alpha;
	//printf("%lf\n",tau);
	// Mise à jour de la position et du temps
    x_cur[0] = x_cur[0]+radius*sin(theta)*cos(phi);
  	x_cur[1] = x_cur[1]+radius*sin(theta)*sin(phi);
	x_cur[2] = x_cur[2]+radius*cos(theta);

	t_cur = t_cur - tau;
	
	#ifdef CHEMINS
		fprintf(fptr,"%lf,%lf,%lf",x_cur[0],x_cur[1],x_cur[2]);
	#endif
	if(t_cur < 0){
		printf("Erreur temps\n");
	}
	//Test position
	//Eps-shell Dirichlet
	if(x_cur[2] < eps){
		#ifdef CHEMINS
			fprintf(fptr,"e");
		#endif
		return -alpha*u1(x_cur);
	}	
	//Eps-shell Neumann
	else if(x0-x_cur[2] < eps){
		#ifdef CHEMINS
			fprintf(fptr,";");
		#endif
		return fif(x_cur, t_cur, x0, fptr, k, q0, T, alpha, eps);
	}
	else{
		#ifdef CHEMINS
			fprintf(fptr,";");
		#endif
		return wos(x_cur, t_cur, x0, fptr, k, q0, T, alpha, eps);
	}
}


double fif(double *x_cur, double t_cur, double x0, FILE* fptr, 
	double k, double q0, double T, double alpha, double eps){
	
	//Echantillonnage
	double u = (double) rand() / (double) RAND_MAX;
	double r = x0*(1-sqrt(u));
	
	double u2 = (double) rand() / (double) RAND_MAX;
	double theta_q = 2*M_PI*u2;

	//double x_disc = x_cur[0]+cos(theta_q)*r; 
	double weight = q1(t_cur, k, q0, T);
	//printf("%lf\n", t_cur);

	double v = (double) rand() / (double) RAND_MAX;
	double theta = acos(v);

	double w = (double) rand() / (double) RAND_MAX;
	double phi = 2.0*M_PI*w;

	double s = (double) rand() / (double) RAND_MAX;
	double tirage_propagateur = swf_tabulation_inverse(Propagateur, SWF_QUADRATIC, w);

	//Calcul rayon plus grande sphère
   	double radius = x0;
	
	double tau = radius*radius*tirage_propagateur/alpha;

	// Calcul poids
	weight *= sum_G(tau,t_cur,radius,alpha);
	weight *= 2*M_PI*radius*alpha*t_cur;

	//Mise à jour de la position et du temps
	t_cur = t_cur - tau;
	
    x_cur[0] = x_cur[0]+radius*sin(theta)*cos(phi);
	x_cur[1] = x_cur[1]+radius*sin(theta)*sin(phi);
  	x_cur[2] = x_cur[2]-radius*cos(theta);

	
	#ifdef CHEMINS
		fprintf(fptr,"%lf,%lf,%lf",x_cur[0],x_cur[1],x_cur[2]);
	#endif
	if(t_cur < 0){
		printf("Erreur temps\n");
	}
	//Test position
	//Eps-shell Dirichlet
	if(x_cur[2] < eps){
		#ifdef CHEMINS
			fprintf(fptr,"e");
		#endif
		return alpha*u1(x_cur);
	}	
	//Eps-shell Neumann
	else if(x0-x_cur[2] < eps){
		#ifdef CHEMINS
			fprintf(fptr,";");
		#endif
		return fif(x_cur, t_cur, x0, fptr, k, q0, T, alpha, eps) + weight;
	}
	else{
		#ifdef CHEMINS
			fprintf(fptr,";");
		#endif
		return wos(x_cur, t_cur, x0, fptr, k, q0, T, alpha, eps) + weight;
	}
}


int main() {
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
	double x0 = 4; // m hauteur du slab
	double ua = 1; // °C
	// Table 5 - probe position
	double x_obs[3] = {1.0,1.0,x0};
	double t_obs = 5*T;

	double a,b,c,d,phi;
	double u;
  	compute_constants(T, alpha, q0_si, A_si, ua, x0, k, &a,&b,&c,&d,&phi);

	
	// Paramètres numériques Monte Carlo
	double eps = 0.00001;
	int N = 100000; // nombre d'échantillons

	swf_H3d_tabulate(&SWF_H3D_TABULATE_ARGS_DEFAULT, &Propagateur);

	// Variables observation et MC
	double t_cur = t_obs;
	double sum = 0;
	double sum_2 = 0;

	FILE* fptr = fopen("echantillons_pos_slab.txt","w");

	//Pour chaque réalisation de Monte-Carlo
	for (int i = 0; i < N; i++) {
		double weight = 0;
		#ifdef CHEMINS
			fprintf(fptr,"%lf,%lf,%lf;",x_obs[0],x_obs[1],x_obs[2]);
		#endif
		double x_cur[3] = {x_obs[0],x_obs[1],x_obs[2]};
		double t_cur = t_obs;
		
		//Test position
		//Eps-shell Dirichlet
		if(x_cur[2] < eps){
			weight = u1(x_cur);
		}	
		//Eps-shell Neumann
		else if(x0-x_cur[2] < eps){
			weight = fif(x_cur, t_cur, x0, fptr, k, q0, T, alpha, eps);
		}
		else{
			weight = wos(x_cur, t_cur, x0, fptr, k, q0, T, alpha, eps);
		}
		printf("%lf\n",weight);
		sum += weight;
		sum_2 += pow(weight, 2);
		
	}

	fclose(fptr);

	double average = sum/(double)N;
	double variance=sum_2/(double)N-pow(sum/(double)N,2);
	double std_dev=sqrt(variance/(double)N);

	printf("Estimation : %f +- %f K\n", average, std_dev);
	u = temperature(x_obs, t_obs, T, alpha, q0_si, A_si, ua, x0, k);
	printf("Solution analytique dans le slab : %lf \n",u);

	return 0;
}


