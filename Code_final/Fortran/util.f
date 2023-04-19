
      SUBROUTINE vis_forces(N,f,vis,velo)
      IMPLICIT NONE
      INTEGER N, i
      DOUBLE PRECISION f(N),vis(N),velo(N)   
          DO i=1,N
            f(i) = -vis(i) * velo(i)
          END DO
      END

      SUBROUTINE wind_forces(N,f,vis,velo)
      IMPLICIT NONE
      INTEGER N, i
      DOUBLE PRECISION f(N),vis(N),velo
          DO i=1,N
            f(i) = f(i) -vis(i) * velo
          END DO
      END

      SUBROUTINE add_norm(N,r,delta)
      IMPLICIT NONE
      INTEGER N, k
      DOUBLE PRECISION r(N),delta(N)
        DO k=1,N
          r(k) = r(k) + (delta(k) * delta(k))
        END DO
      END

      FUNCTION forces(Wv,deltav,rv)
      IMPLICIT NONE
      DOUBLE PRECISION forces,Wv,deltav,rv
        forces=Wv*deltav/(rv**3)
        RETURN
      END

