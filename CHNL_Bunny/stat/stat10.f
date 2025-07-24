c
c stat/stat10.f
c (C) 2012 Y, anzai
c
        SUBROUTINE STAT10()
        INCLUDE '../par.f'
        INCLUDE '../common.f'

        OPEN(5000,FILE='stat10.out',STATUS='UNKNOWN')
        DO 5100 K=1,NZ
        DO 5100 J=0,NY
        DO 5100 I=0,NX
        WRITE(5000,5400)X(I),Y(J),X12(I),Y12(J),UXM(I,J,K),UYM(I,J,K)   
 5400   FORMAT('',6E15.6)
 5100   CONTINUE
        CLOSE(5000)

        RETURN
        END



