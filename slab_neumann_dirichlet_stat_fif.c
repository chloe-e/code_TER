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


double u1(double x, double z){
	return sin(x)*exp(-z);
}

double q1(double x_disc[3]){
	return -sin(x_disc[0])*exp(-x_disc[2]);
}


double fif(double *x_cur, double L, double eps, FILE* fptr){
    //printf("%lf,%lf,%lf\n",x_cur[0],x_cur[1],x_cur[2]);	
    double tmp[3] = {0,0,0};

    //Echantillonage angles nouvelle pos
    double v = (double) rand() / (double) RAND_MAX;
    double theta = acos(-2.0*v+1.0);
    double w = (double) rand() / (double) RAND_MAX;
    double phi = 2.0*M_PI*w;

    //Calcul nouveau rayon	
    double radius = x_cur[2];
    if(radius > L){
        printf("Erreur rayon");
        return 1;
    }

    //Mise à jour de la position	
    tmp[0] = x_cur[0]+radius*sin(theta)*cos(phi);
    tmp[1] = x_cur[1]+radius*sin(theta)*sin(phi);
    tmp[2] = x_cur[2]-radius*cos(theta);
	
    #ifdef CHEMINS
	    fprintf(fptr,"%lf,%lf,%lf",x_cur[0],x_cur[1],x_cur[2]);
    #endif
	
    //TEST POSITION
    //Pas de disque
    if(x_cur[2]<=L/2){

        x_cur[0] = tmp[0];
        x_cur[1] = tmp[1];
        x_cur[2] = tmp[2];
	//printf("pas de disque\n");

        //Eps-shell Dirichlet
        if(x_cur[2] < eps){
            #ifdef CHEMINS
                fprintf(fptr,"e");
            #endif

            return u1(x_cur[0],0);
        }
        else{
            #ifdef CHEMINS
                fprintf(fptr,";");
            #endif

            return fif(x_cur,L,eps,fptr);
        }
    }
    //Disque
    else if(x_cur[2]>L/2){
        //printf("disque\n");
        //Echantillonnage sur le disque
        double u = (double) rand() / (double) RAND_MAX;
        double r = L*(1-sqrt(u));
        
        double s = (double) rand() / (double) RAND_MAX;
        double theta_q = 2*M_PI*s;

        double x_disc[3] = {x_cur[0]+cos(theta_q)*r,x_cur[1]+sin(theta_q)*r,L}; 
        double weight = q1(x_disc)*radius/2;

        //TEST EN DEHORS
        if(tmp[2]>L){
            //printf("dehors\n");
            tmp[0] = x_cur[0] + (tmp[0]-x_cur[0])*(L-x_cur[2])/(tmp[2]-x_cur[2]);
            tmp[1] = x_cur[1] + (tmp[1]-x_cur[1])*(L-x_cur[2])/(tmp[2]-x_cur[2]);
            tmp[2] = L;
        }

        x_cur[0] = tmp[0];
        x_cur[1] = tmp[1];
        x_cur[2] = tmp[2];
        
        //Eps-shell Dirichlet
        if(x_cur[2] < eps){
            #ifdef CHEMINS
                fprintf(fptr,"e");
            #endif
            //printf("%lf\n",x_cur[2]);
            return u1(x_cur[0],0);
        }	
        else{
            #ifdef CHEMINS
                fprintf(fptr,";");
            #endif
            return fif(x_cur,L,eps,fptr) + weight;
        }
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
	     
  double sum = 0;
  double sum_2 = 0;

  FILE* fptr = fopen("echantillons_pos_slab_fif.txt","w");

  // Pour chaque réalisation de Monte-Carlo
  for (int i = 0; i < N; i++) {
	double weight = 0;
	#ifdef CHEMINS
		fprintf(fptr,"%lf,%lf,%lf;",x_obs[0],x_obs[1],x_obs[2]);
	#endif
	double x_cur[3] = {x_obs[0],x_obs[1],x_obs[2]};
    
	//Test position
	//Eps-shell Dirichlet
	if(x_cur[2] < eps){
		weight = u1(x_cur[0],0);
	}	
	else{
		weight = fif(x_cur,L,eps,fptr);
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


