c
c solver/mkp.f
c (C) 2004 K. Fukagata 
c
      SUBROUTINE MKP(ALPHA)
      INCLUDE '../par.f'
      INCLUDE '../common.f'
      COMMON /PP/ PHI
      COMMON /PT/ PHIT
      COMMON /PG/ GXPHI,GYPHI,GZPHI
      REAL*8 PHI(0:NX+1,0:NY+1,0:NZ+1)
      REAL*8 PHITMP(0:NX*2+1,0:NY+1,0:NZ+1)
      REAL*8 GXPHI(0:NX+1,0:NY+1,0:NZ+1)
      REAL*8 GYPHI(0:NX+1,0:NY+1,0:NZ+1)
      REAL*8 GZPHI(0:NX+1,0:NY+1,0:NZ+1)
      COMPLEX*16 PHIT(0:NX,0:NY+1,0:NZ+1)
c
c 1. RHS of Poisson eq.
c      
      DO 110 K=1,NZ
      DO 110 J=1,NY         
      DO 110 I=1,NX
         PHI(I,J,K)=(
     +        +(UX(I,J,K)-UX(I-1,J,K))/DX
     +        +(UY(I,J,K)-UY(I,J-1,K))/DY(J)
     +        +(UZ(I,J,K)-UZ(I,J,K-1))/DZ)/ALPHA/DT
 110  CONTINUE
c Reflection
      DO 200 K=1,NZ
      DO 200 J=1,NY
         DO 210 I=1,NX
            PHITMP(I,J,K)=PHI(I,J,K)
 210     CONTINUE
         DO 220 I=NX+1,NX*2
            PHITMP(I,J,K)=PHI(NX*2+1-I,J,K)
 220     CONTINUE
 200  CONTINUE
c
      DO 120 J=1,NY
      DO 120 I=1,NX*2
         PHITMP(I,J,0)=PHITMP(I,J,NZ)
 120  CONTINUE
      DO 130 K=1,NZ
      DO 130 J=1,NY
         PHITMP(0,J,K)=PHITMP(NX*2,J,K)
 130  CONTINUE
      DO 140 J=1,NY
         PHITMP(0,J,0)=PHITMP(NX*2,J,NZ)
 140  CONTINUE
c
c 2. Solving Poisson eq. using trigonometric expansion 
c
c P2S (external): 2D-FFT from physical to spectral spaces
c S2P (external): 2D-FFT from spectral to physical spaces
c SPN (below)   : Solver for non-zero mode
c SP0B (below)  : Solver for zero mode
c PHIBC (below) : Boundaary condition for PHI (Neumann)
c
      CALL P2S(PHITMP,PHIT)
      CALL SPN()
      CALL SP0B()
      CALL S2P(PHIT,PHITMP)
      DO 250 K=1,NZ
      DO 250 J=1,NY
      DO 250 I=1,NX
         PHI(I,J,K)=PHITMP(I,J,K)
 250  CONTINUE
      CALL PHIBC() 
c
c 3. Compute Gradient Phi
c
      DO 310 K=1,NZ
      DO 310 J=1,NY
      DO 310 I=1,NX
         GXPHI(I,J,K)=(PHI(I+1,J,K)-PHI(I,J,K))/DX
         GYPHI(I,J,K)=(PHI(I,J+1,K)-PHI(I,J,K))/DY12(J)
         GZPHI(I,J,K)=(PHI(I,J,K+1)-PHI(I,J,K))/DZ
 310  CONTINUE
c
c 4. Velocity correction
c
      DO 410 K=1,NZ
      DO 410 J=1,NY
      DO 410 I=1,NX
         UX(I,J,K)=UX(I,J,K)-ALPHA*DT*GXPHI(I,J,K) 
         UY(I,J,K)=UY(I,J,K)-ALPHA*DT*GYPHI(I,J,K)
         UZ(I,J,K)=UZ(I,J,K)-ALPHA*DT*GZPHI(I,J,K)
c
c         GXPH(I,J,K)=-ALPHA*DT*GXPHI(I,J,K)
c         GYPH(I,J,K)=-ALPHA*DT*GYPHI(I,J,K)
 410  CONTINUE     
c
c 5. Pressure 
c
      DO 510 K=1,NZ
      DO 510 J=1,NY
      DO 510 I=1,NX
         P(I,J,K)=P(I,J,K)+PHI(I,J,K)
         GXP(I,J,K)=GXP(I,J,K)+GXPHI(I,J,K)
         GYP(I,J,K)=GYP(I,J,K)+GYPHI(I,J,K)
         GZP(I,J,K)=GZP(I,J,K)+GZPHI(I,J,K)
 510  CONTINUE
      RETURN
      END


      SUBROUTINE SPN()
c
c CTRIDIAG (external): TDMA solver for complex matrix
c
      INCLUDE '../par.f'
      INCLUDE '../common.f'
      COMMON /PT/ PHIT
      REAL*8 AA(NY,3)
      COMPLEX*16 PHIT(0:NX,0:NY+1,0:NZ+1),B(NY),D(NY)
*poption parallel,tlocal(I,J,AA,B)
      DO 210 K=0,NZ-1
      DO 210 I=0,NX
         IF((I+1)*(K+1).GT.1)THEN
            DO 220 J=1,NY
               AA(J,1)=AC(J,1)
               AA(J,2)=AC(J,2)-KXZ(I,K)
               AA(J,3)=AC(J,3)
               B(J)=PHIT(I,J,K)
               D(J)=PHIT(I,J,K)
 220        CONTINUE
            AA(1,2)=AA(1,2)+AA(1,1)
            AA(1,1)=0.D0
            AA(NY,2)=AA(NY,2)+AA(NY,3)
            AA(NY,3)=0.D0
            CALL CTRIDIAG(AA,B,NY)
            DO 230 J=1,NY
               PHIT(I,J,K)=B(J)
 230        CONTINUE
         ENDIF
 210  CONTINUE
      RETURN
      END


      SUBROUTINE SP0B()
      INCLUDE '../par.f'
      INCLUDE '../common.f'
      COMMON /PT/ PHIT
      REAL*8 AA(NY,3)
      COMPLEX*16 PHIT(0:NX,0:NY+1,0:NZ+1),D(NY)
      DO 200 J=1,NY
         AA(J,1)=AC(J,1)
         AA(J,2)=AC(J,2)
         AA(J,3)=AC(J,3)
         D(J)=PHIT(0,J,0)
 200  CONTINUE
      PHIT(0,0,0)=(0.D0,0.D0)
      PHIT(0,1,0)=(0.D0,0.D0)
      DO 230 J=1,NY-1
         PHIT(0,J+1,0)=(D(J)-AA(J,1)*PHIT(0,J-1,0)-AA(J,2)*PHIT(0,J,0))
     +        /AA(J,3)
 230  CONTINUE
      PHIT(0,NY+1,0)=PHIT(0,NY,0)
      RETURN
      END


      SUBROUTINE PHIBC()
      INCLUDE '../par.f'
      COMMON /PP/ PHI
      REAL*8 PHI(0:NX+1,0:NY+1,0:NZ+1)
      DO 100 K=1,NZ
      DO 100 I=0,NX
         PHI(I,0,K)=PHI(I,1,K)
         PHI(I,NY+1,K)=PHI(I,NY,K)
 100  CONTINUE
      DO 110 J=1,NY
      DO 110 I=0,NX
         PHI(I,J,NZ+1)=PHI(I,J,1)
         PHI(I,J,0   )=PHI(I,J,NZ)
 110  CONTINUE
      DO 120 K=0,NZ
      DO 120 J=1,NY
         PHI(NX+1,J,K)=PHI(1, J,K)
         PHI(0,   J,K)=PHI(NX,J,K)
 120  CONTINUE
      DO 130 J=1,NY
         PHI(0,J,0)=PHI(1,J,NZ) 
         PHI(0,J,NZ+1)=PHI(1,J,1)
         PHI(NX+1,J,0)=PHI(NX,J,NZ)
         PHI(NX+1,J,NZ+1)=PHI(NX,J,1)
 130  CONTINUE
      RETURN
      END
