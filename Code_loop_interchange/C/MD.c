/*
 *  Simple molecular dynamics code.
 *  2022
 */
#include <stdio.h>
#include <math.h>
#include "coord.h"
#include <immintrin.h>


void vis_forces(int N,double *f, double *vis, double *vel);
void add_norms(int N,double *r, double *delta);
double forces(double W, double delta, double r);
void wind_forces(int N,double *f, double *vis, double vel);

void evolve(int count,double dt){
  int step;
  int i,j,k,l;
  int have_collided;
  double size;
/*
 * Loop over timesteps.
 */
  for(step = 1;step<=count;step++){
    printf("timestep %d\n",step);
    printf("collisions %d\n",collisions);

    /* set the viscosity term in the force calculation */
    for(j=0;j<Ndim;j++){


      vis_forces(Nbody,f[j],vis,velo[j]);
    }
    /* add the wind term in the force calculation */
    for(j=0;j<Ndim;j++){
      

      wind_forces(Nbody,f[j],vis,wind[j]);
    }
    
    /* calculate distance from central mass */
    for(k=0;k<Nbody;k++){
      r[k] = 0.0;
    }
    
    for(i=0;i<Ndim;i++){

      for(k=0;k<Nbody;k++){
        r[k] += (pos[i][k] * pos[i][k]);
      }
      // add_norms(Nbody,r,pos[i]);

    }
    for(k=0;k<Nbody;k++){
      r[k] = sqrt(r[k]);
    }
    /* calculate central force */
    for(l=0;l<Ndim;l++){
      for(i=0;i<Nbody;i++){
      

        f[l][i] = f[l][i] - forces(G*mass[i]*M_central,pos[l][i],r[i]);
        
      }
	  }
    
    /* calculate pairwise separation of the particles */

    for(l=0;l<Ndim;l++){
      k = 0;
      for(i=0;i<Nbody;i++){
        for(j=i+1;j<Nbody;j++){
        // k = (2 * Nbody - 1 - i) * i/2 + j - i - 1;
        delta_pos[l][k] = pos[l][i] - pos[l][j];
        k = k + 1;
        }
      }
    }

    /* calculate norm of separation vector */
    for(k=0;k<Npair;k++){
      delta_r[k] = 0.0;
    }
    
    for(i=0;i<Ndim;i++){

      // add_norms(Npair,delta_r,delta_pos[i]);
      for(k=0;k<Npair;k++){
        delta_r[k] += (delta_pos[i][k] * delta_pos[i][k]);
      }

    }
    for(k=0;k<Npair;k++){
      delta_r[k] = sqrt(delta_r[k]);
    }

    /*
    * add pairwise forces.
    */

   
    k = 0;
    for(i=0;i<Nbody;i++){
      for(j=i+1;j<Nbody;j++){
        
        size = radius[i] + radius[j];
        have_collided=0;

        for(l=0;l<Ndim;l++){
          // k = (2 * Nbody - 1 - i) * i/2 + j - i - 1;
          /*  flip force if close in */
          if( delta_r[k] >= size ){
            f[l][i] = f[l][i] - forces(G*mass[i]*mass[j],delta_pos[l][k],delta_r[k]);
            f[l][j] = f[l][j] + forces(G*mass[i]*mass[j],delta_pos[l][k],delta_r[k]);
            
          }else{
            f[l][i] = f[l][i] + forces(G*mass[i]*mass[j],delta_pos[l][k],delta_r[k]);
            f[l][j] = f[l][j] - forces(G*mass[i]*mass[j],delta_pos[l][k],delta_r[k]);
            
		        have_collided=1;
          }
        }
        if( have_collided == 1 ){
          collisions++;
        }
        k = k + 1;
      }
    }

    /* update positions */
    for(j=0;j<Ndim;j++){
      for(i=0;i<Nbody;i++){
      
        pos[j][i] = pos[j][i] + dt * velo[j][i];
      }
    }

    /* update velocities */
    for(j=0;j<Ndim;j++){
      for(i=0;i<Nbody;i++){
      
        velo[j][i] = velo[j][i] + dt * (f[j][i]/mass[i]);
      }
    }
  }

}




