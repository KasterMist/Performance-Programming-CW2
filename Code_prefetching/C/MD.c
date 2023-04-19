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

  __assume_aligned(r, 16);
  __assume_aligned(delta_r, 16);
  __assume_aligned(mass, 16);
  __assume_aligned(radius, 16);
  __assume_aligned(vis, 16);

  __assume_aligned(f[0], 16);
  __assume_aligned(pos[0], 16);
  __assume_aligned(velo[0], 16);
  __assume_aligned(delta_pos[0], 16);

  for(step = 1;step<=count;step++){
    printf("timestep %d\n",step);
    printf("collisions %d\n",collisions);

    /* set the viscosity term in the force calculation */
    /* add the wind term in the force calculation */
    for(j=0;j<Ndim;j++){
      // #pragma vector aligned
      vis_forces(Nbody,f[j],vis,velo[j]);
      wind_forces(Nbody,f[j],vis,wind[j]);
    }


    /* calculate distance from central mass */
    // #pragma vector aligned
    for(k=0;k<Nbody;k++){
      _mm_prefetch(&pos[0][k], _MM_HINT_T0);
      _mm_prefetch(&pos[1][k], _MM_HINT_T0);
      _mm_prefetch(&pos[2][k], _MM_HINT_T0);
      _mm_prefetch(&r[k], _MM_HINT_T0);

      r[k] = (pos[0][k] * pos[0][k]);
      r[k] += (pos[1][k] * pos[1][k]);
      r[k] += (pos[2][k] * pos[2][k]);
      r[k] = sqrt(r[k]);
    }


    /* calculate central force */
    // #pragma vector aligned
    for(l=0;l<Ndim;l++){
      for(i=0;i<Nbody;i++){
        _mm_prefetch(&f[l][i], _MM_HINT_T0);
        
        f[l][i] = f[l][i] - forces(G*mass[i]*M_central,pos[l][i],r[i]);
        
      }
	  }
    
    /* calculate pairwise separation of the particles */
    
    
    // #pragma vector aligned
    for(l=0;l<Ndim;l++){
      k = 0;
      for(i=0;i<Nbody;i++){
        for(j=i+1;j<Nbody;j++){
          _mm_prefetch(&delta_pos[l][k], _MM_HINT_T0);
          _mm_prefetch(&pos[l][i], _MM_HINT_T0);
          _mm_prefetch(&pos[l][j], _MM_HINT_T0);
          delta_pos[l][k] = pos[l][i] - pos[l][j];
          k = k + 1;
        }
      
      }
    }

    /* calculate norm of separation vector */
    // #pragma vector aligned
    for(k=0;k<Npair;k++){
      _mm_prefetch(&delta_r[k], _MM_HINT_T0);
      _mm_prefetch(&delta_pos[0][k], _MM_HINT_T0);
      _mm_prefetch(&delta_pos[1][k], _MM_HINT_T0);
      _mm_prefetch(&delta_pos[2][k], _MM_HINT_T0);

      delta_r[k] = (delta_pos[0][k] * delta_pos[0][k]);
      delta_r[k] += (delta_pos[1][k] * delta_pos[1][k]);
      delta_r[k] += (delta_pos[2][k] * delta_pos[2][k]);
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
          for(l=0;l<Ndim;l++){
            f[l][i] = f[l][i] - forces(G*mass[i]*mass[j],delta_pos[l][k],delta_r[k]);
            f[l][j] = f[l][j] + forces(G*mass[i]*mass[j],delta_pos[l][k],delta_r[k]);
          }
        }else{
          for(l=0;l<Ndim;l++){
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
    /* update velocities */
    // #pragma vector aligned
    for(j=0;j<Ndim;j++){
      for(i=0;i<Nbody;i++){
        _mm_prefetch(&pos[j][i], _MM_HINT_T0);
        _mm_prefetch(&velo[j][i], _MM_HINT_T0);
        _mm_prefetch(&f[j][i], _MM_HINT_T0);
        _mm_prefetch(&mass[i], _MM_HINT_T0);
        _mm_prefetch(&dt, _MM_HINT_T0);

        pos[j][i] = pos[j][i] + dt * velo[j][i];
        velo[j][i] = velo[j][i] + dt * (f[j][i]/mass[i]);
      }
    }


    
  }

}




