c
c stat/stat2.f
c (C) 2012 K, Nakayama
c
        SUBROUTINE STAT2()
        INCLUDE '../par.f'
        INCLUDE '../common.f'

        NST=NST+1

        DO 1000 I=1,NX
        DO 1000 J=1,NY
        DO 1000 K=1,NZ
c           UXM(I,J,K)=UXM(I,J,K)+UX(I,J,K)
c           UYM(I,J,K)=UYM(I,J,K)+UY(I,J,K)
c           UZM(I,J,K)=UZM(I,J,K)+UZ(I,J,K)
c           PM(I,J,K) =PM(I,J,K) +P(I,J,K)
c           WM(I,J,K) =WM(I,J,K) +W(I,J,K)
           DFM(I,J,K) =DFM(I,J,K) +FDF(I,J,K)
           LFM(I,J,K) =LFM(I,J,K) +FLF(I,J,K)
           DPM(I,J,K) =DPM(I,J,K) +FDP(I,J,K)
           LPM(I,J,K) =LPM(I,J,K) +FLP(I,J,K)
 1000   CONTINUE


        IF (NNN/III.EQ.1) THEN
c           UXM(I,J,K)=UXM(I,J,K)/DBLE(NST)
c           UYM(I,J,K)=UYM(I,J,K)/DBLE(NST)
c           UZM(I,J,K)=UZM(I,J,K)/DBLE(NST)
c           PM(I,J,K)=PM(I,J,K)/DBLE(NST)
c           WM(I,J,K)=WM(I,J,K)/DBLE(NST)
           DFM(I,J,K) =DFM(I,J,K)/DBLE(NST) 
           LFM(I,J,K) =LFM(I,J,K)/DBLE(NST)
           DPM(I,J,K) =DPM(I,J,K)/DBLE(NST)
           LPM(I,J,K) =LPM(I,J,K)/DBLE(NST)
           
        OPEN(10,FILE='Out/stat.out',STATUS='UNKNOWN')
        DO 101 J=1,NY
        DO 101 I=1,NX
        DO 101 K=1,NZ
        WRITE(10,600)X(I),Y(J),DFM(I,J,K),LFM(I,J,K),
     +               DPM(I,J,K),LPM(I,J,K)
c        WRITE(10,600)X(I),Y(J),UXM1(I,J,K),UYM1(I,J,K),
c     +               PM1(I,J,K),WM(I,J,K) 

 600    FORMAT('',7E15.6)
 101    CONTINUE
        CLOSE(10)

        ENDIF


        RETURN
        END
