c
c init/grid.f
c (C) 2004 K. Fukagata 
c
      SUBROUTINE GRID()
      INCLUDE '../par.f'
      INCLUDE '../common.f'
C------------------------------------------------------------------------------
C 1. Define grid
C------------------------------------------------------------------------------
c x direction
      DO 100 I=0,NX+1
         X  (I)=XL*(DBLE(I)-0.5D0)/DBLE(NX)
         X12(I)=XL*(DBLE(I)      )/DBLE(NX)
 100  CONTINUE
      DX=XL/DBLE(NX)
c y direction
      DO 120 J=0,NY+1
         Y  (J)=YL*(DBLE(J)-0.5D0)/DBLE(NY)
         Y12(J)=YL*(DBLE(J)      )/DBLE(NY)
 120  CONTINUE
      DY=YL/DBLE(NY)
      DY12=DY
c z direction
      DO 150 K=0,NZ+1
         Z  (K)=ZL*(DBLE(K)-0.5D0)/DBLE(NZ)
         Z12(K)=ZL*(DBLE(K)      )/DBLE(NZ)
 150  CONTINUE
      DZ=ZL/DBLE(NZ)
C------------------------------------------------------------------------------
C 2. Constants for FFT
C------------------------------------------------------------------------------
c x direction
      DO 200 I=0,NX
         KX(I)=DBLE(I)/(XL*2.D0)*PI*2.D0
 200  CONTINUE
c z direction
      DO 210 K=0,NZ
         IF (K .LE. NZ/2) THEN
            KZ(K)=DBLE(K)
         ELSE
            KZ(K)=DBLE(K-NZ)
         ENDIF
 210  CONTINUE
      DO 220 K=0,NZ-1
         KZ(K)=KZ(K)/ZL*PI*2.0D0
 220  CONTINUE
c Laplacian
      DO 230 K=0,NZ
      DO 230 I=0,NX
         KXZ(I,K)=2.D0*(1.D0-COS(KX(I)*DX))/DX/DX
     +          +2.D0*(1.D0-COS(KZ(K)*DZ))/DZ/DZ
 230  CONTINUE
C------------------------------------------------------------------------------
C 3. Matrix coefficients
C------------------------------------------------------------------------------
c For centered veriables (p, ux, uz)
      DO 300 J=1,NY
         AC(J,1)=1.D0/DY12(J-1)/DY(J)
         AC(J,3)=1.D0/DY12(J)  /DY(J)
         AC(J,2)=-AC(J,1)-AC(J,3)
 300  CONTINUE
c For edged variables (uy)
      DO 310 J=1,NY
         AE(J,1)=1.D0/DY(J)  /DY12(J)
         AE(J,3)=1.D0/DY(J+1)/DY12(J)
         AE(J,2)=-AE(J,1)-AE(J,3)
 310  CONTINUE      
      RETURN
      END


