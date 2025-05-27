#include <math.h>
#include <stdbool.h>
#include <stdio.h>
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

double fif(double *x_cur,double L,double eps, double tp_init, FILE* fptr);

double u1(double x, double z){
	return sin(x)*exp(-z);
}

double q1(double x, double z){
	return -sin(x)*exp(-z);
}

double wos(double *x_cur, double L, double eps, double tp_init, FILE *fptr){

	//Echantillonnage
	double u = (double) rand() / (double) RAND_MAX;
	double phi = u*2.0*M_PI;

	double v = (double) rand() / (double) RAND_MAX;
	double theta = acos(-2.0*v+1.0);
	
 	// Calcul du rayon de la plus grande sphère inscrite dans le cube
      	double radius = fmin(x_cur[2],L-x_cur[2]);
      	if(radius>1){
		printf("Erreur rayon");
		return 1;
	}
	// Mise à jour de la position 
      	x_cur[0] = x_cur[0]+radius*sin(theta)*cos(phi);
      	x_cur[1] = x_cur[1]+radius*sin(theta)*sin(phi);
      	x_cur[2] = x_cur[2]+radius*cos(theta);
	
	#ifndef CHEMINS
		fprintf(fptr,"%lf,%lf,%lf",x_cur[0],x_cur[1],x_cur[2]);
	#endif

	//Test position
	//Eps-shell Dirichlet
	if(x_cur[2] < eps){
		#ifndef CHEMINS
			fprintf(fptr,"e");
		#endif
		return u1(x_cur[0],0);
	}	
	//Eps-shell Neumann
	else if(L-x_cur[2] < eps){
		#ifndef CHEMINS
			fprintf(fptr,";");
		#endif
		return fif(x_cur,L,eps,tp_init,fptr);
	}
	else{
		#ifndef CHEMINS
			fprintf(fptr,";");
		#endif
		return wos(x_cur,L,eps,tp_init,fptr);
	}
}


double fif(double *x_cur, double L, double eps, double tp_init, FILE* fptr){
	
	//Echantillonnage
	double u = (double) rand() / (double) RAND_MAX;
	double r = L*(1-sqrt(u));
	
	double u2 = (double) rand() / (double) RAND_MAX;
	double theta_q = 2*M_PI*u2;

	double x_disc = x_cur[0]+cos(theta_q)*r; 
	double weight = q1(x_disc,L)*L/2;

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
	
	#ifndef CHEMINS
		fprintf(fptr,"%lf,%lf,%lf",x_cur[0],x_cur[1],x_cur[2]);
	#endif
	
	//Test position
	//Eps-shell Dirichlet
	if(x_cur[2] < eps){
		#ifndef CHEMINS
			fprintf(fptr,"e");
		#endif
		return u1(x_cur[0],0);
	}	
	//Eps-shell Neumann
	else if(L-x_cur[2] < eps){
		#ifndef CHEMINS
			fprintf(fptr,";");
		#endif
		return fif(x_cur,L,eps,tp_init,fptr)+weight;
	}
	else{
		#ifndef CHEMINS
			fprintf(fptr,";");
		#endif
		return wos(x_cur,L,eps,tp_init,fptr)+weight;
	}
}

int main() {
  // Paramètres physiques du système
  double L = 2; //Hauteur du slab

  // Observation
  double x_obs[3] = {1,0.5,1}; // point initial au centre du cube
  // Paramètres numériques
  double eps = 0.00001;
  int N = 100000; // nombre d'échantillons
  double tp_init = 100000000;  

  double sum = 0;
  double sum_2 = 0;

  FILE* fptr = fopen("echantillons_pos_slab.txt","w");

  // Pour chaque réalisation de Monte-Carlo
  for (int i = 0; i < N; i++) {
	double weight = 0;
	#ifndef CHEMINS
		fprintf(fptr,"%lf,%lf,%lf;",x_obs[0],x_obs[1],x_obs[2]);
	#endif
    	double x_cur[3] = {x_obs[0],x_obs[1],x_obs[2]};
    
	//Test position
	//Eps-shell Dirichlet
	if(x_cur[2] < eps){
		weight = u1(x_cur[0],0);
	}	
	//Eps-shell Neumann
	else if(L-x_cur[2] < eps){
		weight = fif(x_cur,L,eps,tp_init,fptr);
	}
	else{
		weight = wos(x_cur,L,eps,tp_init,fptr);
	}
    sum += weight;
    sum_2 += pow(weight, 2);
       
  }
  fclose(fptr);
  double average = sum/(double)N;
  double variance=sum_2/(double)N-pow(sum/(double)N,2);
  double std_dev=sqrt(variance/(double)N);

  printf("Estimation : %f +- %f K\n", average, std_dev);
  printf("Solution analytique au centre du cube : %f \n",u1(x_obs[0],x_obs[2]));
  return 0;
}


