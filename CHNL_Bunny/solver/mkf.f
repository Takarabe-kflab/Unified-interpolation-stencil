c
c solver/mkf.f
c (C) 2004 K. Fukagata 
c
      SUBROUTINE MKF()
      INCLUDE '../par.f'
      INCLUDE '../common.f'
      DO 200 K=1,NZ
      DO 200 J=1,NY
      DO 200 I=1,NX
c
c Fx
c
         FX(I,J,K)=
c -d(ux ux)/dx
     +        -(
     +         (UX(I,  J,K)+UX(I+1,J,K))*(UX(I,  J,K)+UX(I+1,J,K))
     +        -(UX(I-1,J,K)+UX(I,  J,K))*(UX(I-1,J,K)+UX(I,  J,K))
     +        )/4.D0/DX
c -d(uy ux)/dy
     +        -(
     +         (UY(I,J,  K)+UY(I+1,J,  K))*(UX(I,J,  K)+UX(I,J+1,K))
     +        -(UY(I,J-1,K)+UY(I+1,J-1,K))*(UX(I,J-1,K)+UX(I,J,  K))
     +        )/4.D0/DY(J)
c -d(uz ux)/dz
     +        -(
     +         (UZ(I,J,K  )+UZ(I+1,J,K  ))*(UX(I,J,K  )+UX(I,J,K+1))
     +        -(UZ(I,J,K-1)+UZ(I+1,J,K-1))*(UX(I,J,K-1)+UX(I,J,K  ))
     +        )/4.D0/DZ
c -dP/dx
     +        +GPM
c
c Fy
c
         FY(I,J,K)=
c -d(ux uy)/dx
     +        -(
     +         (DY(J)*UX(I,  J,K)+DY(J+1)*UX(I,  J+1,K))
     +              *(UY(I,  J,K)        +UY(I+1,J,  K))
     +        -(DY(J)*UX(I-1,J,K)+DY(J+1)*UX(I-1,J+1,K))
     +              *(UY(I-1,J,K)        +UY(I,  J,  K))
     +        )/4.D0/DX/DY12(J)
c -d(uy uy)/dy
     +        -(
     +         (UY(I,J,  K)+UY(I,J+1,K))*(UY(I,J,K  )+UY(I,J+1,K))
     +        -(UY(I,J-1,K)+UY(I,J,  K))*(UY(I,J-1,K)+UY(I,J,K  ))
     +        )/4.D0/DY12(J)
c -d(uz uy)/dz
     +        -(
     +         (DY(J)*UZ(I,J,K  )+DY(J+1)*UZ(I,J+1,K  ))
     +              *(UY(I,J,K  )        +UY(I,J,  K+1))
     +        -(DY(J)*UZ(I,J,K-1)+DY(J+1)*UZ(I,J+1,K-1))
     +              *(UY(I,J,K-1)        +UY(I,J,  K))
     +        )/4.D0/DZ/DY12(J)
c
c Fz
c
         FZ(I,J,K)=
c -d(ux uz)/dx
     +        -(
     +         (UX(I,  J,K)+UX(I,  J,K+1))*(UZ(I,  J,K)+UZ(I+1,J,K))
     +        -(UX(I-1,J,K)+UX(I-1,J,K+1))*(UZ(I-1,J,K)+UZ(I,  J,K))
     +        )/4.D0/DX
c -d(uy uz)/dy
     +        -(
     +         (UY(I,J,  K)+UY(I,J,  K+1))*(UZ(I,J,  K)+UZ(I,J+1,K))
     +        -(UY(I,J-1,K)+UY(I,J-1,K+1))*(UZ(I,J-1,K)+UZ(I,J,  K))
     +        )/4.D0/DY(J)
c -d(uz uz)/dz
     +        -(
     +         (UZ(I,J,K  )+UZ(I,J,K+1))*(UZ(I,J,K  )+UZ(I,J,K+1))
     +        -(UZ(I,J,K-1)+UZ(I,J,K  ))*(UZ(I,J,K-1)+UZ(I,J,K  ))
     +        )/4.D0/DZ
c
c FT
c
         FT(I,J,K)=0.D0
c -d(ux T)/dx
     +        -(
     +         UX(I,  J,K)*(TT(I,  J,K)+TT(I+1,J,K))
     +        -UX(I-1,J,K)*(TT(I-1,J,K)+TT(I,  J,K))
     +        )/2.D0/DX
c -d(uy T)/dy
     +        -(
     +         UY(I,J,  K)*(TT(I,J,  K)+TT(I,J+1,K))
     +        -UY(I,J-1,K)*(TT(I,J-1,K)+TT(I,J,  K))
     +        )/2.D0/DY(J)
c -d(uz T)/dz
     +        -(
     +         UZ(I,J,K  )*(TT(I,J,K  )+TT(I,J,K+1))
     +        -UZ(I,J,K-1)*(TT(I,J,K-1)+TT(I,J,K  ))
     +        )/2.D0/DZ
c
c Diffusion term
c
         DUX(I,J,K)=(
     +         (UX(I-1,J,K)-2.D0*UX(I,J,K)+UX(I+1,J,K))/DX/DX
     +        +(UX(I,J,K-1)-2.D0*UX(I,J,K)+UX(I,J,K+1))/DZ/DZ
     +        +AC(J,1)*UX(I,J-1,K)
     +        +AC(J,2)*UX(I,J,  K)
     +        +AC(J,3)*UX(I,J+1,K)
     +        )/REB
         DUY(I,J,K)=(
     +         (UY(I-1,J,K)-2.D0*UY(I,J,K)+UY(I+1,J,K))/DX/DX
     +        +(UY(I,J,K-1)-2.D0*UY(I,J,K)+UY(I,J,K+1))/DZ/DZ
     +        +AE(J,1)*UY(I,J-1,K)
     +        +AE(J,2)*UY(I,J,  K)
     +        +AE(J,3)*UY(I,J+1,K)
     +        )/REB
         DUZ(I,J,K)=(
     +         (UZ(I-1,J,K)-2.D0*UZ(I,J,K)+UZ(I+1,J,K))/DX/DX
     +        +(UZ(I,J,K-1)-2.D0*UZ(I,J,K)+UZ(I,J,K+1))/DZ/DZ
     +        +AC(J,1)*UZ(I,J-1,K)
     +        +AC(J,2)*UZ(I,J,  K)
     +        +AC(J,3)*UZ(I,J+1,K)
     +        )/REB
         DTT(I,J,K)=(
     +         (TT(I-1,J,K)-2.D0*TT(I,J,K)+TT(I+1,J,K))/DX/DX
     +        +(TT(I,J,K-1)-2.D0*TT(I,J,K)+TT(I,J,K+1))/DZ/DZ
     +        +AC(J,1)*TT(I,J-1,K)
     +        +AC(J,2)*TT(I,J,  K)
     +        +AC(J,3)*TT(I,J+1,K)
     +        )/REB/PRL
 200  CONTINUE
      RETURN
      END




