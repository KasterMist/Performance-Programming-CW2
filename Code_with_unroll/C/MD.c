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
    // for(j=0;j<Ndim;j++){
      // vis_forces(Nbody,f[j],vis,velo[j]);
      vis_forces(Nbody,f[0],vis,velo[0]);
      vis_forces(Nbody,f[1],vis,velo[1]);
      vis_forces(Nbody,f[2],vis,velo[2]);
    // }
    /* add the wind term in the force calculation */
    // for(j=0;j<Ndim;j++){
      
      // wind_forces(Nbody,f[j],vis,wind[j]); 
      wind_forces(Nbody,f[0],vis,wind[0]); 
      wind_forces(Nbody,f[1],vis,wind[1]); 
      wind_forces(Nbody,f[2],vis,wind[2]); 
    // }
    /* calculate distance from central mass */
    for(k=0;k<Nbody;k++){
      r[k] = 0.0;
    }
    
    // for(i=0;i<Ndim;i++){

      // for(k=0;k<Nbody;k++){
      //   r[k] += (pos[i][k] * pos[i][k]);
      // }

      for(k=0;k<Nbody;k++){
        r[k] += (pos[0][k] * pos[0][k]);
      }
      for(k=0;k<Nbody;k++){
        r[k] += (pos[1][k] * pos[1][k]);
      }
      for(k=0;k<Nbody;k++){
        r[k] += (pos[2][k] * pos[2][k]);
      }
      // add_norms(Nbody,r,pos[i]);

    // }
    for(k=0;k<Nbody;k++){
      r[k] = sqrt(r[k]);
    }
    /* calculate central force */
    // for(l=0;l<Ndim;l++){
      for(i=0;i<Nbody;i++){
        // f[l][i] = f[l][i] - forces(G*mass[i]*M_central,pos[l][i],r[i]);
        f[0][i] = f[0][i] - forces(G*mass[i]*M_central,pos[0][i],r[i]);
        f[1][i] = f[1][i] - forces(G*mass[i]*M_central,pos[1][i],r[i]);
        f[2][i] = f[2][i] - forces(G*mass[i]*M_central,pos[2][i],r[i]);
      }
	  // }
    
    /* calculate pairwise separation of the particles */
    k = 0;
    
    for(i=0;i<Nbody;i++){
      for(j=i+1;j<Nbody;j++){
        // for(l=0;l<Ndim;l++){
        // delta_pos[l][k] = pos[l][i] - pos[l][j];
        delta_pos[0][k] = pos[0][i] - pos[0][j];
        delta_pos[1][k] = pos[1][i] - pos[1][j];
        delta_pos[2][k] = pos[2][i] - pos[2][j];
        // }
      k = k + 1;
      }
    }

    /* calculate norm of separation vector */
    for(k=0;k<Npair;k++){
      delta_r[k] = 0.0;
    }
    
    // for(i=0;i<Ndim;i++){

      // add_norms(Npair,delta_r,delta_pos[i]);
      // for(k=0;k<Npair;k++){
      //   delta_r[k] += (delta_pos[i][k] * delta_pos[i][k]);

      // }
      for(k=0;k<Npair;k++){
        delta_r[k] += (delta_pos[0][k] * delta_pos[0][k]);
      }
      for(k=0;k<Npair;k++){
        delta_r[k] += (delta_pos[1][k] * delta_pos[1][k]);
      }
      for(k=0;k<Npair;k++){
        delta_r[k] += (delta_pos[2][k] * delta_pos[2][k]);
      }

    // }
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


        if(delta_r[k] >= size){
          // for(l=0;l<Ndim;l++){
            // f[l][i] = f[l][i] - forces(G*mass[i]*mass[j],delta_pos[l][k],delta_r[k]);
            // f[l][j] = f[l][j] + forces(G*mass[i]*mass[j],delta_pos[l][k],delta_r[k]);
          // }
          f[0][i] = f[0][i] - forces(G*mass[i]*mass[j],delta_pos[0][k],delta_r[k]);
          f[0][j] = f[0][j] + forces(G*mass[i]*mass[j],delta_pos[0][k],delta_r[k]);

          f[1][i] = f[1][i] - forces(G*mass[i]*mass[j],delta_pos[1][k],delta_r[k]);
          f[1][j] = f[1][j] + forces(G*mass[i]*mass[j],delta_pos[1][k],delta_r[k]);

          f[2][i] = f[2][i] - forces(G*mass[i]*mass[j],delta_pos[2][k],delta_r[k]);
          f[2][j] = f[2][j] + forces(G*mass[i]*mass[j],delta_pos[2][k],delta_r[k]);
        }else{
          // for(l=0;l<Ndim;l++){
            // f[l][i] = f[l][i] + forces(G*mass[i]*mass[j],delta_pos[l][k],delta_r[k]);
            // f[l][j] = f[l][j] - forces(G*mass[i]*mass[j],delta_pos[l][k],delta_r[k]);
		        // have_collided=1; 

            f[0][i] = f[0][i] + forces(G*mass[i]*mass[j],delta_pos[0][k],delta_r[k]);
            f[0][j] = f[0][j] - forces(G*mass[i]*mass[j],delta_pos[0][k],delta_r[k]);
		        have_collided=1; 

            f[1][i] = f[1][i] + forces(G*mass[i]*mass[j],delta_pos[1][k],delta_r[k]);
            f[1][j] = f[1][j] - forces(G*mass[i]*mass[j],delta_pos[1][k],delta_r[k]);
		        have_collided=1; 

            f[2][i] = f[2][i] + forces(G*mass[i]*mass[j],delta_pos[2][k],delta_r[k]);
            f[2][j] = f[2][j] - forces(G*mass[i]*mass[j],delta_pos[2][k],delta_r[k]);
		        have_collided=1; 
          // }
        }

        // for(l=0;l<Ndim;l++){
        //   /*  flip force if close in */
        //   if( delta_r[k] >= size ){
        //     f[l][i] = f[l][i] - forces(G*mass[i]*mass[j],delta_pos[l][k],delta_r[k]);
        //     f[l][j] = f[l][j] + forces(G*mass[i]*mass[j],delta_pos[l][k],delta_r[k]);
            
        //   }else{
        //     f[l][i] = f[l][i] + forces(G*mass[i]*mass[j],delta_pos[l][k],delta_r[k]);
        //     f[l][j] = f[l][j] - forces(G*mass[i]*mass[j],delta_pos[l][k],delta_r[k]);
            
		    //     have_collided=1;
        //   }
        // }


        if( have_collided == 1 ){
          collisions++;
        }
        k = k + 1;
      }
    }

    /* update positions */
    // for(j=0;j<Ndim;j++){
      // for(i=0;i<Nbody;i++){
      //   pos[j][i] = pos[j][i] + dt * velo[j][i];
      // }
      for(i=0;i<Nbody;i++){
        pos[0][i] = pos[0][i] + dt * velo[0][i];
      }
      for(i=0;i<Nbody;i++){
        pos[1][i] = pos[1][i] + dt * velo[1][i];
      }
      for(i=0;i<Nbody;i++){
        pos[2][i] = pos[2][i] + dt * velo[2][i];
      }
    // }

    /* update velocities */
    // for(j=0;j<Ndim;j++){
      // for(i=0;i<Nbody;i++){
      //   velo[j][i] = velo[j][i] + dt * (f[j][i]/mass[i]);
      // }
      for(i=0;i<Nbody;i++){
        velo[0][i] = velo[0][i] + dt * (f[0][i]/mass[i]);
      }
      for(i=0;i<Nbody;i++){
        velo[1][i] = velo[1][i] + dt * (f[1][i]/mass[i]);
      }
      for(i=0;i<Nbody;i++){
        velo[2][i] = velo[2][i] + dt * (f[2][i]/mass[i]);
      }
    // }
  }

}




