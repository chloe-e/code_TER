#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>

/*               (2,2,2)
 *      __________
 *     /         /|
 *    /    T1   / |
 *   /_________/  |   Temperature T2 sur les autres faces
 *   |         |  |   
 *   |     x   |  |      z
 *   |         |  /      |  y
 *   |         | /       | /
 *   |_________|/        |/____ x
 *(0,0,0)
 */

int main() {
  // Paramètres physiques du système
  double L = 2; // Largeur du carré
  double temperature_1 = 100; //température sur le bord haut
  double temperature_2 = 300; //température des autres bords
  // Observation
  double x_obs[3] = {1,1,1}; // point initial au centre du cube
  // Paramètres numériques
  double eps = 0.00001;
  int N = 10; // nombre d'échantillons
  
  double sum = 0;
  double sum_2 = 0;
  FILE *fptr = fopen("echantillons_pos_cube.txt","w");  
  

  // Pour chaque réalisation de Monte-Carlo
  for (int i = 0; i < N; i++) {
    double x_cur[3] = {x_obs[0],x_obs[1],x_obs[2]};
    double weight = 0.;
    bool stop = false;
    fprintf(fptr,"%lf,%lf,%lf;",x_obs[0],x_obs[1],x_obs[2]);
  
    while(!stop)
    {
      // Echantillonnages servant à calculer la prochaine position
      double u = (double) rand() / (double) RAND_MAX;
      double theta = acos(2.0*u-1.0);
      double v = (double) rand() / (double) RAND_MAX;
      double phi = 2.0*M_PI*v;

      // Calcul du rayon de la plus grande sphère inscrite dans le cube
      double radius = fmin(fmin(fmin(x_cur[0],x_cur[1]),x_cur[2]),
                           fmin(fmin(L-x_cur[0],L-x_cur[1]),L-x_cur[2]));

      // Mise à jour de la position 
      x_cur[0] = x_cur[0]+radius*sin(theta)*cos(phi);
      x_cur[1] = x_cur[1]+radius*sin(theta)*sin(phi);
      x_cur[2] = x_cur[2]+radius*cos(theta);

      // Eventuel arrêt si la paroi est atteinté
      if(    x_cur[0]<eps || x_cur[0]>L-eps
          || x_cur[1]<eps || x_cur[1]>L-eps
          || x_cur[2]<eps )
      {
        weight = temperature_2;
        stop = true;
      }
      if(x_cur[2]>L-eps) {
        weight = temperature_1;
        stop = true;
      }

      fprintf(fptr,"%lf,%lf,%lf",x_cur[0],x_cur[1],x_cur[2]);
      fprintf(fptr,stop ? "e" : ";");
    }
    sum += weight;
    sum_2 += pow(weight, 2);
  }  // end for i=1,N

  fclose(fptr);
  
  double average = sum/(double)N;
  double variance=sum_2/(double)N-pow(sum/(double)N,2);
  double std_dev=sqrt(variance/(double)N);
  printf("Estimation : %f +- %f K\n", average, std_dev);
  printf("Solution analytique au centre du cube : %f \n",(1.*temperature_1+5.*temperature_2)/6.);
  return 0;
}
