c
c stat/stat0.f
c (C) 2006 K. Fukagata 
c
      SUBROUTINE STAT0()
      INCLUDE '../par.f'
      INCLUDE '../common.f'
      CALL UBC()
c
c Instanteneous bulk and center velocities
c   
      UXB=0.D0
      UZB=0.D0
      UXC=0.D0
      DO 100 K=1,NZ
      DO 100 I=1,NX
         UXC=UXC+0.5D0*(UX(I,NY/2,K)+UX(I,NY/2+1,K))
 100  CONTINUE
      UXC=UXC/DBLE(NX*NZ)
      DO 120 K=1,NZ
      DO 120 J=1,NY
      DO 120 I=1,NX
         UXB=UXB+DY(J)*UX(I,J,K)
         UZB=UZB+DY(J)*UZ(I,J,K)
 120  CONTINUE
      UXB=0.5D0*UXB/DBLE(NX*NZ)
      UZB=0.5D0*UZB/DBLE(NX*NZ)
c
c     Output
c     
      WRITE(11,610)TID,UXB,UXC,GPM,DTDY,UZB
 610  FORMAT(' ',6E15.7)
c     
      RETURN
      END
      
