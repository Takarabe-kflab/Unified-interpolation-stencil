c output/omega1.f
c (c)2012 Y.Anzai
c
      SUBROUTINE OMEGA1()
      INCLUDE '../par.f'
      INCLUDE '../common.f'
c
      OPEN(61,FILE='omega1.out',STATUS='UNKNOWN')
c
      DO 101 J=40,443
      DO 101 I=315,NX
      DO 101 K=1,NZ
      WRITE(61,600)X(I),Y(J),X12(I),Y12(J),W(I,J,K)

 600  FORMAT('',7E15.5)
 101  CONTINUE
      CLOSE(61)
      RETURN
      END
