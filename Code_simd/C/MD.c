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

// __m256d approximate_sqrt(__m256d);
// void simd_assign(double*, double*, int);

__m256d approximate_sqrt(__m256d x) {
    __m256d xhalf = _mm256_mul_pd(x, _mm256_set1_pd(0.5));
    __m256d invsqrt = _mm256_rsqrt_pd(x);
    __m256d three = _mm256_set1_pd(3.0);
    __m256d invsqrt_half = _mm256_mul_pd(invsqrt, _mm256_set1_pd(0.5));
    __m256d y = _mm256_mul_pd(invsqrt_half, _mm256_mul_pd(invsqrt_half, x));
    y = _mm256_mul_pd(three, _mm256_sub_pd(y, _mm256_set1_pd(1.0)));
    return _mm256_mul_pd(xhalf, _mm256_mul_pd(invsqrt, y));
}

void simd_assign(double* a, double* b, int size) {
    for (int i = 0; i < size; i += 4) {
      // __m256d a_vec = _mm256_load_pd(&a[i]); 
      __m256d b1_vec = _mm256_load_pd(&b[i]);
      __m256d b2_vec = _mm256_load_pd(&b[i+size]);
      __m256d b3_vec = _mm256_load_pd(&b[i+2*size]);
      __m256d b1_mul = _mm256_mul_pd(b1_vec, b1_vec);
      __m256d b2_mul = _mm256_mul_pd(b2_vec, b2_vec);
      __m256d b3_mul = _mm256_mul_pd(b3_vec, b3_vec);
      
      __m256d result_vec_1 = _mm256_add_pd(b1_mul, b2_mul); 
      __m256d result_vec = _mm256_add_pd(b3_vec, result_vec_1);
      // __m256d result_sqrt_vec = _mm256_sqrt_pd(result_vec);
      _mm256_store_pd(&a[i], result_vec); 
    }
}

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
    simd_assign(r, pos[0], Nbody);
    for(k=0;k<Nbody;k++){

      // r[k] = (pos[0][k] * pos[0][k]);
      // r[k] += (pos[1][k] * pos[1][k]);
      // r[k] += (pos[2][k] * pos[2][k]);
      r[k] = sqrt(r[k]);
    }


    // for(k=0;k<Nbody;k++){
    //   r[k] = 0.0;
    // }
    
    
    // for(i=0;i<Ndim;i++){
    //   // simd_assign(r, pos[i], Nbody);
    //   for(k=0;k<Nbody;k++){
    //     r[k] += (pos[i][k] * pos[i][k]);
    //   }
    // //   // add_norms(Nbody,r,pos[i]);

    // }
    // for(k=0;k<Nbody;k++){
    //   r[k] = sqrt(r[k]);
    // }


    /* calculate central force */
    // #pragma vector aligned
    for(l=0;l<Ndim;l++){
      for(i=0;i<Nbody;i++){
      
        f[l][i] = f[l][i] - forces(G*mass[i]*M_central,pos[l][i],r[i]);
        
      }
	  }
    
    /* calculate pairwise separation of the particles */
    
    
    // #pragma vector aligned
    for(l=0;l<Ndim;l++){
      k = 0;
      for(i=0;i<Nbody;i++){
        for(j=i+1;j<Nbody;j++){
          delta_pos[l][k] = pos[l][i] - pos[l][j];
          k = k + 1;
        }
      
      }
    }

    /* calculate norm of separation vector */

    // #pragma vector aligned
    // for(k=0;k<Npair;k++){
    //   delta_r[k] = (delta_pos[0][k] * delta_pos[0][k]);
    //   delta_r[k] += (delta_pos[1][k] * delta_pos[1][k]);
    //   delta_r[k] += (delta_pos[2][k] * delta_pos[2][k]);
    //   delta_r[k] = sqrt(delta_r[k]);
    // }

    for(k=0;k<Npair;k++){
      delta_r[k] = 0.0;
    }
    
    for(i=0;i<Ndim;i++){
      // simd_assign(delta_r, delta_pos[i], Npair);
      add_norms(Npair,delta_r,delta_pos[i]);
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
    // #pragma vector aligned
    // for(j=0;j<Ndim;j++){
    //   for(i=0;i<Nbody;i++){
    //     pos[j][i] = pos[j][i] + dt * velo[j][i];
    //     velo[j][i] = velo[j][i] + dt * (f[j][i]/mass[i]);
    //   }
    // }

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




