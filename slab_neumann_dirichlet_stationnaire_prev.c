#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <limits.h>
#include <stdlib.h>

/*               (2,2,2)
 *      __________
 *     /         /|
 *    /    q1   / |
 *   /_________/  |   Temperature u1 sur le plan opposé (ox,oy)
 *   |         |  |   
 *   |     x   |  |      z
 *   |         |  /      |  y
 *   |         | /       | /
 *   |_________|/        |/____ x
 *(0,0,0)
 */

double fif(double *x_cur,double L,double eps);

double u1(double x, double z){
	return sin(x)*exp(-z);
}

double q1(double x, double z){
	return -sin(x)*exp(-z);
}

double wos(double *x_cur, double L,double eps){

	//Echantillonnage
	double u = (double) rand() / (double) RAND_MAX;
	double phi = u*2.0*M_PI;

	double v = (double) rand() / (double) RAND_MAX;
	double theta = acos(-2.0*v+1.0);
	
 	// Calcul du rayon de la plus grande sphère inscrite dans le cube
      	double radius = fmin(x_cur[2],L-x_cur[2]);
      	//printf("%lf,%lf,%lf\n",x_cur[0],x_cur[1],x_cur[2]);
	// Mise à jour de la position 
      	x_cur[0] = x_cur[0]+radius*sin(theta)*cos(phi);
      	x_cur[1] = x_cur[1]+radius*sin(theta)*sin(phi);
      	x_cur[2] = x_cur[2]+radius*cos(theta);
	
	//Test position
	//Eps-shell Dirichlet
	if(x_cur[2] < eps){
		x_cur[2] = 0;
		return u1(x_cur[0],0);
	}	
	//Eps-shell Neumann
	else if(L-x_cur[2] < eps){
		x_cur[2] = L;
		return fif(x_cur,L,eps);
	}
	else{
		return wos(x_cur,L,eps);
	}
}


double fif(double *x_cur, double L, double eps){
	//Echantillonnage
	double u = (double) rand() / (double) RAND_MAX;
	double r = L*(1-sqrt(u));
	
	double u2 = (double) rand() / (double) RAND_MAX;
	double theta_q = 2*M_PI*u2;

	double v = (double) rand() / (double) RAND_MAX;
	double theta = acos(v);

	double w = (double) rand() / (double) RAND_MAX;
	double phi = 2.0*M_PI*w;

	//Calcul nouveau rayon	
      	double radius = L;

	//Mise à jour de la position	
      	x_cur[0] = x_cur[0]+radius*sin(theta)*cos(phi);
      	x_cur[1] = x_cur[1]+radius*sin(theta)*sin(phi);
      	x_cur[2] = x_cur[2]-radius*cos(theta);

	//Test position
	//Eps-shell Dirichlet
	if(x_cur[2] < eps){
		x_cur[2] = 0;
		return u1(x_cur[0],0);
	}	
	//Eps-shell Neumann
	else if(L-x_cur[2] < eps){
		x_cur[2] = L;
		return fif(x_cur,L,eps);
	}
	else{
		double x_disc = x_cur[0]+cos(theta_q)*r; 
		return wos(x_cur,L,eps) + q1(x_disc,L);
	}


}

int main() {
  // Paramètres physiques du système
  double L = 2; //Hauteur du slab

  // Observation
  double x_obs[3] = {1,0.5,L}; // point initial au centre du cube
  // Paramètres numériques
  double eps = 0.00001;
  int N = 100000; // nombre d'échantillons

  double sum = 0;
  double sum_2 = 0;

  

  // Pour chaque réalisation de Monte-Carlo
  for (int i = 0; i < N; i++) {
    double x_cur[3] = {x_obs[0],x_obs[1],x_obs[2]};
    double weight = wos(x_cur,L,eps);

    sum += weight;
    sum_2 += pow(weight, 2);
       
  }

  double average = sum/(double)N;
  double variance=sum_2/(double)N-pow(sum/(double)N,2);
  double std_dev=sqrt(variance/(double)N);

  printf("Estimation : %f +- %f K\n", average, std_dev);
  printf("Solution analytique au centre du cube : %f \n",u1(x_obs[0],x_obs[2]));
  return 0;
}


