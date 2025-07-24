c
c stat/stat11.f
c (C) 2012 Y, anzai
c
        SUBROUTINE STAT11()
        INCLUDE '../par.f'
        INCLUDE '../common.f'

        OPEN(6000,FILE='stat11.out',STATUS='UNKNOWN')
        DO 6100 K=1,NZ
        DO 6100 J=1,NY
        DO 6100 I=1,NX
c        WRITE(6000,6200)X(I),Y(J),FPAY6(I,J,K),FPAX6(I,J,K),
c     +                 FPAY5(I,J,K),FPAX5(I,J,K)
c        WRITE(6000,6200)X(I),Y(J),X12(I),Y12(J),UYrms(I,J,K),
c     +                  REYST(I,J,K)
        WRITE(6000,6200)X(I),Y(J),X12(I),PM1(I,J,K),
     +                  P2M(I,J,K),PRrms(I,J,K)
 6200   FORMAT('',6E15.7)
 6100   CONTINUE
        CLOSE(6000)

        RETURN
        END
