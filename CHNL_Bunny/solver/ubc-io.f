c
c solver/ubc-io.f
c (C) 2004-2007 K. Fukagata
c (C) 2010      H. Hasebe 
c (C) 2011      N. Tomiyama
c (C) 2011      K. Nakayama
c (C) 2025      I. Takarabe
c      
      SUBROUTINE UBC()
      INCLUDE '../par.f'
      INCLUDE '../common.f'
c
c Inflow/Outflow in x direction
c
c x=0 boundary
      DO 130 K=0,NZ+1
      DO 130 J=0,NY+1
         UX(0,J,K)=VXI(J,K)
         UY(0,J,K)=2.D0*VYI(J,K)-UY(1,J,K)
         UZ(0,J,K)=2.D0*VZI(J,K)-UZ(1,J,K)
         P(0,J,K)=P(1,J,K)
         TT(0,J,K)=2.D0*TTI(J,K)-TT(1,J,K)         
c x=XL boundary
         UX(NX,J,K)=VXO(J,K)
         UY(NX+1,J,K)=2.D0*VYO(J,K)-UY(NX,J,K)
         UZ(NX+1,J,K)=2.D0*VZO(J,K)-UZ(NX,J,K)
         P(NX+1,J,K)=P(NX,J,K)
         TT(NX+1,J,K)=2.D0*TTO(J,K)-TT(NX,J,K)
c for 4th CDS
         UX(-1,J,K)=UX(0,J,K)
         UY(-1,J,K)=UY(0,J,K)
         UZ(-1,J,K)=UZ(0,J,K)
         P(-1,J,K)=P(0,J,K)
         TT(-1,J,K)=TT(0,J,K)
ccccc
         UX(NX+2,J,K)=UX(NX+1,J,K)
         UY(NX+2,J,K)=UY(NX+1,J,K)
         UZ(NX+2,J,K)=UZ(NX+1,J,K)
         P(NX+2,J,K)=P(NX+1,J,K)
         TT(NX+2,J,K)=TT(NX+1,J,K)
 130  CONTINUE
c
c Wall boundary conditions
c
      DO 110 K=0,NZ+1
      DO 110 I=0,NX+1
c Lower wall
         UX(I,0,K)=UX(I,1,K)
         UY(I,0,K)=0.D0
         UZ(I,0,K)=UZ(I,1,K)
         P (I,0,K)=P(I,1,K)
         TT(I,0,K)=TT(I,1,K)
c Upper wall
         UX(I,NY+1,K)=UX(I,NY,K)
         UY(I,NY,K)=0.D0
         UZ(I,NY+1,K)=UZ(I,NY,K)
         P (I,NY+1,K)=P(I,NY,K)
         TT(I,NY+1,K)=TT(I,NY,K)
 110  CONTINUE
c
c Periodic boundary conditions in z direction
c
      DO 120 J=0,NY+1
      DO 120 I=0,NX+1
c z=0 boundary
         UX(I,J,0)=UX(I,J,NZ)
         UY(I,J,0)=UY(I,J,NZ)
         UZ(I,J,0)=UZ(I,J,NZ)
         P (I,J,0)=P (I,J,NZ)
         TT(I,J,0)=TT(I,J,NZ)
c z=ZL boundary
         UX(I,J,NZ+1)=UX(I,J,1)
         UY(I,J,NZ+1)=UY(I,J,1)
         UZ(I,J,NZ+1)=UZ(I,J,1)
         P (I,J,NZ+1)=P (I,J,1)
         TT(I,J,NZ+1)=TT(I,J,1)
c for 4th CDS
ccccc
         UX(I,J,-1)=UX(I,J,NZ-1)
         UY(I,J,-1)=UY(I,J,NZ-1)
         UZ(I,J,-1)=UZ(I,J,NZ-1)
         P (I,J,-1)=P (I,J,NZ-1)
         TT(I,J,-1)=TT(I,J,NZ-1)
ccccc
         UX(I,J,-2)=UX(I,J,NZ-2)
         UY(I,J,-2)=UY(I,J,NZ-2)
         UZ(I,J,-2)=UZ(I,J,NZ-2)
         P (I,J,-2)=P (I,J,NZ-2)
         TT(I,J,-2)=TT(I,J,NZ-2)
ccccc
         UX(I,J,NZ+2)=UX(I,J,2)
         UY(I,J,NZ+2)=UY(I,J,2)
         UZ(I,J,NZ+2)=UZ(I,J,2)
         P (I,J,NZ+2)=P (I,J,2)
         TT(I,J,NZ+2)=TT(I,J,2)
ccccc
         UX(I,J,NZ+3)=UX(I,J,3)
         UY(I,J,NZ+3)=UY(I,J,3)
         UZ(I,J,NZ+3)=UZ(I,J,3)
         P (I,J,NZ+3)=P (I,J,3)
         TT(I,J,NZ+3)=TT(I,J,3)
 120  CONTINUE
      RETURN
      END
