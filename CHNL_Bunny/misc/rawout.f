c
c misc/rawout.f
c (C) 2004 K. Fukagata 
c
      SUBROUTINE RAWOUT()
      INCLUDE '../par.f'
      INCLUDE '../common.f'
      WRITE(*,*)'writing ml3'
      OPEN(28,FILE='ml3',FORM='UNFORMATTED')
      WRITE(28)TID
      DO 100 K=1,NZ
      DO 100 J=1,NY
      DO 100 I=0,NX+1
         WRITE(28)
     +            UX(I,J,K),UY(I,J,K),UZ(I,J,K),P(I,J,K),TT(I,J,K)
 100  CONTINUE
      CLOSE(28)
      RETURN
      END
