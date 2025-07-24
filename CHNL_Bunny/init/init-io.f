c
c init/init-io.f
c (C) 2004-2024 K. Fukagata 
c
      SUBROUTINE INIT
      INCLUDE '../par.f'
      INCLUDE '../common.f'
      WRITE(*,*)'************************************************'
      WRITE(*,*)'*       Basic Cylinder Flow                    *'
      WRITE(*,*)'*       2024-05-15 by K. Fukagata              *'
      WRITE(*,*)'************************************************'
c
c 1. Read parameters 
c
      WRITE(*,*)'init: Reading inp.csv'
      OPEN(28,FILE='inp.csv',STATUS='OLD')
      READ(28,*)NNN
      READ(28,*)IIC
      READ(28,*)IIO
      READ(28,*)DT
      READ(28,*)REB
      READ(28,*)PRL
      CLOSE(28)
c
      WRITE(*,*)'***************** PARAMETERS *******************'
      WRITE(*,*)'           Number of iteration', NNN
      WRITE(*,*)'              Statistics every', IIC
      WRITE(*,*)'       Output field data every', IIO
      WRITE(*,*)'                       Delta t', DT
      WRITE(*,*)'Diameter-based Reynolds number', REB
      WRITE(*,*)'                Prandtl number', PRL
      WRITE(*,*)'************************************************'
c
c 2. Initial velocity 
c     

      DO 100 K=0,NZ+1
      DO 100 J=0,NY+1
      DO 100 I=0,NX+1
         UX(I,J,K)=1.D0 + 1.D-1*DSIN(2.D0*PI*Y(I)/YL)
     +        + 1.D-1*DSIN(2.D0*PI*Z(I)/ZL)
         UY(I,J,K)=0.D0
         UZ(I,J,K)=0.D0
         P(I,J,K)=0.D0
         TT(I,J,K)=0.D0
 100  CONTINUE 
c 
c 3. Initial inlet B.C.
c
      DO 500 K=0,NZ+1
      DO 500 J=0,NY+1
         VXI(J,K)=UX(0,J,K)
         VYI(J,K)=0.5D0*(UY(0,J,K)+UY(1,J,K))
         VZI(J,K)=0.5D0*(UZ(0,J,K)+UZ(1,J,K))
         TTI(J,K)=0.5D0*(TT(0,J,K)+TT(1,J,K))
         VXO(J,K)=UX(NX,J,K)
         VYO(J,K)=0.5D0*(UY(NX,J,K)+UY(NX+1,J,K))
         VZO(J,K)=0.5D0*(UZ(NX,J,K)+UZ(NX+1,J,K))
         TTO(J,K)=0.5D0*(TT(NX,J,K)+TT(NX+1,J,K))
 500  CONTINUE
      CALL UBC()
c
c 5. Zero initialization    
c
      DO 90 K=0,NZ+1
      DO 90 J=0,NY+1
      DO 90 I=0,NX+1
         DUX(I,J,K)=0.D0
         DUY(I,J,K)=0.D0
         DUZ(I,J,K)=0.D0
         DTT(I,J,K)=0.D0
         FX(I,J,K)=0.D0
         FY(I,J,K)=0.D0
         FZ(I,J,K)=0.D0
         FT(I,J,K)=0.D0
c         FFXIB(I,J,K)=0.D0
c         FFYIB(I,J,K)=0.D0
 90   CONTINUE
c
c 4. Pressure gradient
c
      DO 80 K=1,NZ
      DO 80 J=1,NY         
      DO 80 I=1,NX
         GXP(I,J,K)=(P(I+1,J,K)-P(I,J,K))/DX
         GYP(I,J,K)=(P(I,J+1,K)-P(I,J,K))/DY12(J)
         GZP(I,J,K)=(P(I,J,K+1)-P(I,J,K))/DZ
 80   CONTINUE
      RETURN
      END

