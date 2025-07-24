c
c output/omega7.f
c (c)2012 Y.Anzai
c
      SUBROUTINE OMEGA7()
      INCLUDE '../par.f'
      INCLUDE '../common.f'
c
      OPEN(10,FILE='omega7.out',STATUS='UNKNOWN')
c
      DO 101 J=40,443
      DO 101 I=315,NX
      DO 101 K=1,NZ
      WRITE(10,600)X(I),Y(J),X12(I),Y12(J)),W(I,J,K)

 600  FORMAT('',7E15.5)
