c
c misc/chkcfl.f
c (C) 2004 K. Fukagata 
c
      SUBROUTINE CHKCFL()
      INCLUDE '../par.f'
      INCLUDE '../common.f'
      CFLLIM=3.0D0
      CFL=0.D0
      ICFL=0
      JCFL=0
      DO 200 K=1,NZ
      DO 200 J=1,NY
      DO 200 I=1,NX
         CFLTMP=ABS(UX(I-1,J,K)+UX(I,J,K))/2.D0*DT/DX
     +         +ABS(UY(I,J-1,K)+UY(I,J,K))/2.D0*DT/DY(J)
     +         +ABS(UZ(I,J,K-1)+UZ(I,J,K))/2.D0*DT/DZ
         IF(CFLTMP.GT.CFL)THEN
            JCFL=J
            CFL=CFLTMP
         ENDIF
 200  CONTINUE
      WRITE(*,*)'chkcfl: CFL =',CFL,'(@ J=',JCFL,')'
      IF(CFL.GT.CFLLIM)STOP 'CFL too large!'
      RETURN
      END
