c
c misc/chkcont.f
c (C) 2004 K. Fukagata 
c
      SUBROUTINE CHKCONT()
      INCLUDE '../par.f'
      INCLUDE '../common.f'
      PHIMAX=0.D0
      PHIAVG=0.D0
      MAXI=0
      MAXJ=0
      DO 100 K=1,NZ
      DO 100 J=1,NY
      DO 100 I=1,NX
         PHITMP=(
     +         (UX(I,J,K)-UX(I-1,J,K))/DX
     +        +(UY(I,J,K)-UY(I,J-1,K))/DY(J))
c     +        +(UZ(I,J,K)-UZ(I,J,K-1))/DZ
c     +        )
         IF(ABS(PHITMP).GT.PHIMAX)THEN
            PHIMAX=ABS(PHITMP)
            MAXI=I
            MAXJ=J
         ENDIF
         PHIAVG=PHIAVG+ABS(PHITMP)
 100  CONTINUE
      PHIAVG=PHIAVG/DBLE(NX*NY*NZ)
      write(*,*)'chkcont: PHIAVG =',PHIAVG
      write(*,*)'chkcont: PHIMAX =',PHIMAX,'(@ I=',MAXI,')'
      write(*,*)'chkcont: PHIMAX =',PHIMAX,'(@ J=',MAXJ,')'
      RETURN
      END
