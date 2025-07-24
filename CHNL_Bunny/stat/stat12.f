c
c stat/stat12.f
c (C) 2012 Y, anzai
c
        SUBROUTINE STAT12()
        INCLUDE '../par.f'
        INCLUDE '../common.f'

        OPEN(6001,FILE='stat12.out',STATUS='UNKNOWN')
        DO 6400 K=1,NZ
        DO 6400 J=1,NY
        DO 6400 I=1,NX
        WRITE(6001,6300)X(I),Y(J),X12(I),Y12(J),UXrms(I,J,K),
     +                  UYrms(I,J,K) 
 6300   FORMAT('',6E15.6)
 6400   CONTINUE
        CLOSE(6001)

        RETURN
        END
