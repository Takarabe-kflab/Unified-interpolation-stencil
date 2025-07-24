c
c output/omega.f
c (c)2013 Y.ANZAI
c
      SUBROUTINE OMEGA()
      INCLUDE '../par.f'
      INCLUDE '../common.f'
      CHARACTER*6 EXT
c
      WRITE(EXT,600) III
 600  FORMAT(I6)
      DO 102 I=1,6
         IF(EXT(I:I).EQ.' ')EXT(I:I)='0'
 102     CONTINUE
c
      OPEN(10,FILE='ANM/omegaout.'//EXT)
c
      DO 101 J=40,443
      DO 101 I=315,NX
      DO 101 K=1,NZ
      WRITE(10,601)X(I),Y(J),X12(I),Y12(J),W(I,J,K)
c
 601  FORMAT('',7E15.7)
 101  CONTINUE
      CLOSE(10)
      RETURN
      END
