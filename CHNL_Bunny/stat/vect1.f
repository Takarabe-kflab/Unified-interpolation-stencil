c
c stat/vect1.f
c (C) 2007 K. Fukagata 
c
      SUBROUTINE VECT1()
      INCLUDE '../par.f'
      INCLUDE '../common.f'
      CHARACTER*6 EXT
      CALL UBCWING()
c
c Set file name
c     
      WRITE(EXT,620)III
 620  FORMAT(I6)
      DO 10 I=1,6
         IF(EXT(I:I).EQ.' ')EXT(I:I)='0'
 10   CONTINUE
      OPEN(10,FILE='Out/vect.'//EXT)
      WRITE(*,*)'vect1: writing Out/vect.'//EXT
c
c Instanteneous vector field and spanwise vorticity in a x-y plane
c   
      K=1
      DO 100 J=1,NY
      DO 110 I=1,NX
         UXTMP=(DY(J+1)*UX(I,J+1,K)+DY(J)*UX(I,J,K))/(DY12(J)*2.D0)
         UYTMP=(UY(I+1,J,K)+UY(I,J,K))*0.5D0
         WZTMP=(UY(I+1,J,K)-UY(I,J,K))/DX
     +        -(UX(I,J+1,K)-UX(I,J,K))/DY12(J)

         WRITE(10,610) FXIB(I,J,K), FYIB(I,J,K),
     +   WZTMP,
     +   UXTMP, UYTMP,
     +   UX(I,J,K), UY(I,J,K)
c
 110  CONTINUE
      WRITE(10,*)
 100  CONTINUE
 610  FORMAT(' ',10E15.7)
c     
      CLOSE(10)
      RETURN
      END
      
