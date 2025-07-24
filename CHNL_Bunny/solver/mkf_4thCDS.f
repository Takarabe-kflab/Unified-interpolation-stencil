c
c solver/mkf2015.f
c (C) 2015 K. Fukagata 
c     2018 Y. Nabae
c          apply energy-conservative 4th CDS
c          in the streamwise and spanwise direction
c     2022 Y. Nabae
c          Bug fixed
c
      SUBROUTINE MKF()
      INCLUDE '../par.f'
      INCLUDE '../common.f'
      REAL*8 C98, C18, DXI, DZI, DX3I, DZ3I, REBI, REBPRI
      REAL*8 DXSQ576I, DZSQ576I
      REAL*8 DY12I(NY), DY2I(NY), DY4I(NY)
c
c Pre-define divisions
c
      C98=9.D0/8.D0
      C18=1.D0/8.D0
      DXI=1.D0/DX
      DZI=1.D0/DZ
      DX3I=1.D0/3.D0/DX
      DZ3I=1.D0/3.D0/DZ
      REBI=1.D0/REB
      REBPRI=REBI/PRL
      DXSQ576I=1.D0/576.D0/DX/DX
      DZSQ576I=1.D0/576.D0/DZ/DZ
C$OMP PARALLEL
C$OMP DO
      DO J=1,NY
         DY12I(J)=1.D0/DY12(J)
         DY2I(J)=0.5D0/DY(J)
         DY4I(J)=0.25D0/DY(J)
      ENDDO
C$OMP END DO
c     
c Fx
c
C$OMP DO
      DO 200 K=1,NZ
      DO 200 J=1,NY
      DO 200 I=1,NX
c
c Fx
c
         FX(I,J,K)=
c -d(ux ux)/dx
     +        -(
     +          (
     +            (
     +              (UX(I+1,J,K)+UX(I  ,J,K))*C98
     +             -(UX(I+2,J,K)+UX(I-1,J,K))*C18
     +            )*(UX(I+1,J,K)+UX(I  ,J,K))
     +           -(
     +              (UX(I  ,J,K)+UX(I-1,J,K))*C98
     +             -(UX(I+1,J,K)+UX(I-2,J,K))*C18
     +            )*(UX(I  ,J,K)+UX(I-1,J,K))
     +          )*C98*DXI
     +         -(
     +            (  
     +              (UX(I+2,J,K)+UX(I+1,J,K))*C98
     +             -(UX(I+3,J,K)+UX(I  ,J,K))*C18
     +            )*(UX(I+3,J,K)+UX(I  ,J,K))
     +           -(
     +              (UX(I-1,J,K)+UX(I-2,J,K))*C98
     +             -(UX(I  ,J,K)+UX(I-3,J,K))*C18
     +            )*(UX(I  ,J,K)+UX(I-3,J,K))
     +          )*C18*DX3I
     +        )*0.25D0
c -d(uy ux)/dy
     +        -(
     +           ( (UY(I  ,J  ,K)+UY(I+1,J  ,K))*C98
     +            -(UY(I-1,J  ,K)+UY(I+2,J  ,K))*C18 )
     +           *( UX(I,J  ,K)+UX(I,J+1,K) )
     +          -( (UY(I  ,J-1,K)+UY(I+1,J-1,K))*C98
     +            -(UY(I-1,J-1,K)+UY(I+2,J-1,K))*C18 )
     +           *( UX(I,J-1,K)+UX(I,J  ,K) )
     +         )*DY4I(J)
c -d(uz ux)/dz
     +        -(
     +          (
     +            (
     +              (UZ(I+1,J,K  )+UZ(I  ,J,K  ))*C98
     +             -(UZ(I+2,J,K  )+UZ(I-1,J,K  ))*C18
     +            )*(UX(I  ,J,K+1)+UX(I  ,J,K  ))
     +           -(
     +              (UZ(I+1,J,K-1)+UZ(I  ,J,K-1))*C98
     +             -(UZ(I+2,J,K-1)+UZ(I-1,J,K-1))*C18
     +            )*(UX(I  ,J,K  )+UX(I  ,J,K-1))
     +          )*C98*DZI
     +         -(
     +            (
     +              (UZ(I+1,J,K+1)+UZ(I  ,J,K+1))*C98
     +             -(UZ(I+2,J,K+1)+UZ(I-1,J,K+1))*C18
     +            )*(UX(I  ,J,K+3)+UX(I  ,J,K  ))
     +           -(
     +              (UZ(I+1,J,K-2)+UZ(I  ,J,K-2))*C98
     +             -(UZ(I+2,J,K-2)+UZ(I-1,J,K-2))*C18
     +            )*(UX(I  ,J,K  )+UX(I  ,J,K-3))
     +          )*C18*DZ3I
     +        )*0.25D0
c -dP/dx
     +        +GPM
 200  CONTINUE
c
c Fy
c
C$OMP DO
      DO 300 K=1,NZ
      DO 300 J=1,NY
      DO 300 I=1,NX
         FY(I,J,K)=
c -d(ux uy)/dx
     +        -(
     +          (
     +            (DY(J+1)*UX(I  ,J+1,K)+DY(J)*UX(I  ,J,K))
     +                   *(UY(I+1,J  ,K)      +UY(I  ,J,K))
     +           -(DY(J+1)*UX(I-1,J+1,K)+DY(J)*UX(I-1,J,K))
     +                   *(UY(I  ,J  ,K)      +UY(I-1,J,K))
     +          )*C98*DXI
     +         -(
     +            (DY(J+1)*UX(I+1,J+1,K)+DY(J)*UX(I+1,J,K))
     +                   *(UY(I+3,J  ,K)      +UY(I  ,J,K))
     +           -(DY(J+1)*UX(I-2,J+1,K)+DY(J)*UX(I-2,J,K))
     +                   *(UY(I  ,J  ,K)      +UY(I-3,J,K))
     +          )*C18*DX3I
     +        )*0.25D0*DY12I(J)
c -d(uy uy)/dy
     +        -(
     +         (UY(I,J,  K)+UY(I,J+1,K))*(UY(I,J,K  )+UY(I,J+1,K))
     +        -(UY(I,J-1,K)+UY(I,J,  K))*(UY(I,J-1,K)+UY(I,J,K  ))
     +        )*0.25D0*DY12I(J)
c -d(uz uy)/dz
     +        -(
     +          (
     +            (DY(J+1)*UZ(I,J+1,K  )+DY(J)*UZ(I,J,K  ))
     +                   *(UY(I,J  ,K+1)      +UY(I,J,K  ))
     +           -(DY(J+1)*UZ(I,J+1,K-1)+DY(J)*UZ(I,J,K-1))
     +                   *(UY(I,J  ,K  )      +UY(I,J,K-1))
     +          )*C98*DZI
     +         -(
     +            (DY(J+1)*UZ(I,J+1,K+1)+DY(J)*UZ(I,J,K+1))
     +                   *(UY(I,J  ,K+3)      +UY(I,J,K  ))
     +           -(DY(J+1)*UZ(I,J+1,K-2)+DY(J)*UZ(I,J,K-2))
     +                   *(UY(I,J  ,K  )      +UY(I,J,K-3))
     +          )*C18*DZ3I
     +        )*0.25D0*DY12I(J)
 300  CONTINUE
c
c Fz
c
C$OMP DO
      DO 400 K=1,NZ
      DO 400 J=1,NY
      DO 400 I=1,NX
         FZ(I,J,K)=
c -d(ux uz)/dx
     +        -(
     +          (
     +            (
     +              (UX(I  ,J,K+1)+UX(I,J,K  ))*C98
     +             -(UX(I  ,J,K+2)+UX(I,J,K-1))*C18
     +            )*(UZ(I+1,J,K  )+UZ(I,J,K  ))
     +           -(
     +              (UX(I-1,J,K+1)+UX(I-1,J,K  ))*C98
     +             -(UX(I-1,J,K+2)+UX(I-1,J,K-1))*C18
     +            )*(UZ(I  ,J,K  )+UZ(I-1,J,K  ))
     +          )*C98*DXI
     +         -(
     +            (
     +              (UX(I+1,J,K+1)+UX(I+1,J,K  ))*C98
     +             -(UX(I+1,J,K+2)+UX(I+1,J,K-1))*C18
     +            )*(UZ(I+3,J,K  )+UZ(I  ,J,K  ))
     +           -(
     +              (UX(I-2,J,K+1)+UX(I-2,J,K  ))*C98
     +             -(UX(I-2,J,K+2)+UX(I-2,J,K-1))*C18
     +            )*(UZ(I  ,J,K  )+UZ(I-3,J,K  ))
     +          )*C18*DX3I
     +        )*0.25D0
c -d(uy uz)/dy
     +        -(
     +           ( (UY(I,J  ,K  )+UY(I,J  ,K+1))*C98
     +            -(UY(I,J  ,K-1)+UY(I,J  ,K+2))*C18 )
     +           *( UZ(I,J  ,K)+UZ(I,J+1,K) )
     +          -( (UY(I,J-1,K  )+UY(I,J-1,K+1))*C98
     +            -(UY(I,J-1,K-1)+UY(I,J-1,K+2))*C18 )
     +           *( UZ(I,J-1,K)+UZ(I,J  ,K) )
     +        )*DY4I(J)
c -d(uz uz)/dz
     +        -(
     +          (
     +            (
     +              (UZ(I,J,K+1)+UZ(I,J,K  ))*C98
     +             -(UZ(I,J,K+2)+UZ(I,J,K-1))*C18
     +            )*(UZ(I,J,K+1)+UZ(I,J,K  ))
     +           -(
     +              (UZ(I,J,K  )+UZ(I,J,K-1))*C98
     +             -(UZ(I,J,K+1)+UZ(I,J,K-2))*C18
     +            )*(UZ(I,J,K  )+UZ(I,J,K-1))
     +          )*C98*DZI
     +         -(
     +            (
     +              (UZ(I,J,K+2)+UZ(I,J,K+1))*C98
     +             -(UZ(I,J,K+3)+UZ(I,J,K  ))*C18
     +            )*(UZ(I,J,K+3)+UZ(I  ,J,K  ))
     +           -(
     +              (UZ(I,J,K-1)+UZ(I,J,K-2))*C98
     +             -(UZ(I,J,K  )+UZ(I,J,K-3))*C18
     +            )*(UZ(I,J,K  )+UZ(I,J,K-3))
     +          )*C18*DZ3I
     +        )*0.25D0
 400  CONTINUE
c
c FT
c
C$OMP DO
      DO 500 K=1,NZ
      DO 500 J=1,NY
      DO 500 I=1,NX
         FT(I,J,K)=0.D0
c -d(ux T)/dx
     +        -(
     +          (
     +            UX(I  ,J,K)*(TT(I+1,J,K)+TT(I  ,J,K))
     +           -UX(I-1,J,K)*(TT(I  ,J,K)+TT(I-1,J,K))
     +          )*C98*DXI
     +         -(
     +            UX(I+1,J,K)*(TT(I+3,J,K)+TT(I  ,J,K))
     +           -UX(I-2,J,K)*(TT(I  ,J,K)+TT(I-3,J,K))
     +          )*C18*DX3I
     +        )*0.5D0
c -d(uy T)/dy
     +        -(
     +         UY(I,J,  K)*(TT(I,J,  K)+TT(I,J+1,K))
     +        -UY(I,J-1,K)*(TT(I,J-1,K)+TT(I,J,  K))
     +        )*DY2I(J)
c -d(uz T)/dz
     +        -(
     +          (
     +            UZ(I,J,K  )*(TT(I,J,K+1)+TT(I,J,K  ))
     +           -UZ(I,J,K-1)*(TT(I,J,K  )+TT(I,J,K-1))
     +          )*C98*DZI
     +         -(
     +            UZ(I,J,K+1)*(TT(I,J,K+3)+TT(I,J,K  ))
     +           -UZ(I,J,K-2)*(TT(I,J,K  )+TT(I,J,K-3))
     +          )*C18*DZ3I
     +        )*0.5D0
c Q=0
     +        + 0.D0
 500  CONTINUE
c
c Diffusion term
c
C$OMP DO
      DO 600 K=1,NZ
      DO 600 J=1,NY
      DO 600 I=1,NX
         DUX(I,J,K)=(
ccccc 2nd order ccccc
c     +         (UX(I-1,J,K)-2.D0*UX(I,J,K)+UX(I+1,J,K))*DXSQI
c     +        +(UX(I,J,K-1)-2.D0*UX(I,J,K)+UX(I,J,K+1))*DZSQI
ccccc 4th order ccccc
     +         (-( 1460.D0*UX(I,J,K)
     +            +54.D0*(UX(I-2,J,K)+UX(I+2,J,K)) )
     +          +UX(I-3,J,K)+UX(I+3,J,K)
     +          +783.D0*(UX(I-1,J,K)+UX(I+1,J,K))
     +         )*DXSQ576I
     +        +(-( 1460.D0*UX(I,J,K)
     +            +54.D0*(UX(I,J,K-2)+UX(I,J,K+2)) )
     +          +UX(I,J,K-3)+UX(I,J,K+3) 
     +          +783.D0*(UX(I,J,K-1)+UX(I,J,K+1))
     +         )*DZSQ576I
     +        +AC(J,1)*UX(I,J-1,K)
     +        +AC(J,2)*UX(I,J,  K)
     +        +AC(J,3)*UX(I,J+1,K)
     +        )*REBI
         DUY(I,J,K)=(
ccccc 2nd order ccccc
c     +         (UY(I-1,J,K)-2.D0*UY(I,J,K)+UY(I+1,J,K))*DXSQI
c     +        +(UY(I,J,K-1)-2.D0*UY(I,J,K)+UY(I,J,K+1))*DZSQI
ccccc 4th order ccccc
     +         (-( 1460.D0*UY(I,J,K)
     +            +54.D0*(UY(I-2,J,K)+UY(I+2,J,K)) )
     +          +UY(I-3,J,K)+UY(I+3,J,K)
     +          +783.D0*(UY(I-1,J,K)+UY(I+1,J,K))
     +         )*DXSQ576I
     +        +(-( 1460.D0*UY(I,J,K)
     +            +54.D0*(UY(I,J,K-2)+UY(I,J,K+2)) )
     +          +UY(I,J,K-3)+UY(I,J,K+3) 
     +          +783.D0*(UY(I,J,K-1)+UY(I,J,K+1))
     +         )*DZSQ576I
     +        +AE(J,1)*UY(I,J-1,K)
     +        +AE(J,2)*UY(I,J,  K)
     +        +AE(J,3)*UY(I,J+1,K)
     +        )*REBI
         DUZ(I,J,K)=(
ccccc 2nd order ccccc
c     +         (UZ(I-1,J,K)-2.D0*UZ(I,J,K)+UZ(I+1,J,K))*DXSQI
c     +        +(UZ(I,J,K-1)-2.D0*UZ(I,J,K)+UZ(I,J,K+1))*DZSQI
ccccc 4th order ccccc
     +         (-( 1460.D0*UZ(I,J,K)
     +            +54.D0*(UZ(I-2,J,K)+UZ(I+2,J,K)) )
     +          +UZ(I-3,J,K)+UZ(I+3,J,K)
     +          +783.D0*(UZ(I-1,J,K)+UZ(I+1,J,K))
     +         )*DXSQ576I
     +        +(-( 1460.D0*UZ(I,J,K)
     +            +54.D0*(UZ(I,J,K-2)+UZ(I,J,K+2)) )
     +          +UZ(I,J,K-3)+UZ(I,J,K+3) 
     +          +783.D0*(UZ(I,J,K-1)+UZ(I,J,K+1))
     +         )*DZSQ576I
     +        +AC(J,1)*UZ(I,J-1,K)
     +        +AC(J,2)*UZ(I,J,  K)
     +        +AC(J,3)*UZ(I,J+1,K)
     +        )*REBI
         DTT(I,J,K)=(
ccccc 2nd order ccccc
c     +         (TT(I-1,J,K)-2.D0*TT(I,J,K)+TT(I+1,J,K))*DXSQI
c     +        +(TT(I,J,K-1)-2.D0*TT(I,J,K)+TT(I,J,K+1))*DZSQI
ccccc 4th order ccccc
     +         (-( 1460.D0*TT(I,J,K)
     +            +54.D0*(TT(I-2,J,K)+TT(I+2,J,K)) )
     +          +TT(I-3,J,K)+TT(I+3,J,K)
     +          +783.D0*(TT(I-1,J,K)+TT(I+1,J,K))
     +         )*DXSQ576I
     +        +(-( 1460.D0*TT(I,J,K)
     +            +54.D0*(TT(I,J,K-2)+TT(I,J,K+2)) )
     +          +TT(I,J,K-3)+TT(I,J,K+3) 
     +          +783.D0*(TT(I,J,K-1)+TT(I,J,K+1))
     +         )*DZSQ576I
     +        +AC(J,1)*TT(I,J-1,K)
     +        +AC(J,2)*TT(I,J,  K)
     +        +AC(J,3)*TT(I,J+1,K)
     +        )*REBPRI
 600  CONTINUE
C$OMP END PARALLEL
      RETURN
      END




