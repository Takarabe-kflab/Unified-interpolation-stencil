c
c outmat.f
c --- Output for Matlab
c (C) 2012-2024 K. Fukagata 
c
      SUBROUTINE OUTMAT()
      INCLUDE '../par.f'
      INCLUDE '../common.f'
      CHARACTER*6 EXT
      REAL*8 Q(NX,NY,NZ)
C
      WRITE(EXT,620)III
      DO 10 I=1,6
         IF(EXT(I:I).EQ.' ')EXT(I:I)='0'
 10   CONTINUE

      DO 400 K=1,NZ
      DO 400 J=1,NY
      DO 400 I=1,NX
         Q(I,J,K)=((P(I,J+1,K)-P(I,J,K))/(Y(J+1)-Y(J))
     +        -(P(I,J,K)-P(I,J-1,K))/(Y(J)-Y(J-1)))
     +        /(Y12(J)-Y12(J-1))
     +        +(P(I-1,J,K)-2.D0*P(I,J,K)+P(I+1,J,K))/DX**2
     +        +(P(I,J,K-1)-2.D0*P(I,J,K)+P(I,J,K+1))/DZ**2
 400  CONTINUE
c
c surface data
c
      OPEN(10,FILE='Anim3D/uvwp.'//EXT)
      DO 100 K=1,NZ
      DO 100 J=1,NY
      DO 100 I=1,NX
         WRITE(10,600)
     +           X(I),Y(J),Z(K),     
     +           (UX(I,J,K)+UX(I-1,J,K))*0.5D0, 
     +           (UY(I,J,K)+UY(I,J-1,K))*0.5D0,
     +           (UZ(I,J,K)+UZ(I,J,K-1))*0.5D0,
     +           P(I,J,K), Q(I,J,K)
 100  CONTINUE
      CLOSE(10)
 600  FORMAT(8E12.4)
 620  FORMAT(I6)
      RETURN
      END

