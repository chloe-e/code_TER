#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>

int main() {

  
  double L = 2; // Largeur du carré
  double t1 = 300; //température sur le bord haut
  double t2 = 100; //température sur les 3 côtés restants
  double x_obs[2] = {1,1}; //point initial
  double eps = 0.001;
  // Variables du MC
  int N = 1000000; // nombre d'échantillons
  double sum = 0;
  double sum_2 = 0;

  for (int i = 0; i < N; i++) {
   double x_cur[2] = {x_obs[0],x_obs[1]};
   double weight = 0.;
   bool stop = false;
    while(!stop){

	    // echantillonner r
	    // double u = (double) rand() / (double) RAND_MAX;
	    // double theta = acos(u);
	    double v = (double) rand() / (double) RAND_MAX;
	    double phi = 2*M_PI*v;
	    // Générer la nouvelle position
	    double r = fmin(fmin(x_cur[0],x_cur[1]),fmin(L-x_cur[0],L-x_cur[1]));
	    x_cur[0] = x_cur[0]+r*cos(phi);
	    x_cur[1] = x_cur[1]+r*sin(phi);

	    if(x_cur[0]<eps || x_cur[0]>L-eps || x_cur[1]<eps){
		    weight = t2;
		    stop = true;
	    }
	    if(x_cur[1]>L-eps){
		    weight = t1;
		    stop = true;
	    }
        }
    	sum += weight;
        sum_2 += pow(weight, 2);
   }  // for i=1,N
  double average = sum/(double)N;
  double variance=sum_2/(double)N-pow(sum/(double)N,2);
  double std_dev=sqrt(variance/(double)N);
  printf("Estimation : %e +- %e K\n", average, std_dev);
  printf("Solution analytique : %e \n",(1.*t1+3.*t2)/4.);
  return 0;
}
