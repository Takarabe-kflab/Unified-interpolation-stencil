c
c output/outvelo.f
c  (c)2011 K. Nakayama
c
        SUBROUTINE OUTVELO()
        INCLUDE '../par.f'
        INCLUDE '../common.f'
c
        OPEN(10,FILE='outvelo.out',STATUS='UNKNOWN')
c
        DO 100 J=1,NY
        DO 100 I=1,NX
        DO 100 K=1,NZ
        WRITE(10,600)UX(I,J,K),UY(I,J,K),X(I),Y(J),CD,CL
 600    FORMAT(' ',6E15.7)       
 100  CONTINUE
        CLOSE(10)
        RETURN
        END
