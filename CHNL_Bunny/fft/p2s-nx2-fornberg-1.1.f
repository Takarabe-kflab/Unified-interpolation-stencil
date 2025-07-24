c
c fft/p2s-nx2-fornberg-1.1.f
c (C) 2001-2009 K. Fukagata 
c
      SUBROUTINE P2S(FP,FPT)
      INCLUDE '../par.f'
      COMMON /F1F/ WCX
      COMMON /F2/ IFFLG
      REAL*8 DATAX(0:NX*2-1),FP(0:NX*2+1,0:NY+1,0:NZ+1)
      COMPLEX*16 DATAXT(0:NX-1), DATAZ(0:NZ-1)
      COMPLEX*16 FPT(0:NX,0:NY+1,0:NZ+1), CI, H0, H1
      COMPLEX*16 WCX(0:NX*2+1)
      EQUIVALENCE(DATAX,DATAXT)
      CI=DCMPLX(0.D0,1.D0)
      IF(IFFLG.EQ.0)THEN
         WRITE(*,*)'Making trigonometric table'
         DO 100 I=0,NX*2-1
            WCX(I)=DCMPLX(DCOS(2.D0*PI*DBLE(I)/DBLE(NX*2)),
     +           DSIN(2.D0*PI*DBLE(I)/DBLE(NX*2)))
 100     CONTINUE
         IFFLG=1
      ENDIF
c
      DO 200 J=1,NY
c x direction
         DO 210 K=0,NZ-1
            DO 220 I=0,NX*2-1
               DATAX(I)=FP(I,J,K)
 220        CONTINUE
            CALL FFT(DATAX,NX,-1)
c symmetry
            H0=DATAXT(0)
            H1=DCONJG(DATAXT(0))
            FPT(0,J,K)=0.5D0*(H0+H1)-0.5D0*CI*(H0-H1)
            DO 230 I=1,NX-1
               H0=DATAXT(I)
               H1=DCONJG(DATAXT(NX-I))
               FPT(I,J,K)=0.5D0*(H0+H1)-0.5D0*CI*(H0-H1)*DCONJG(WCX(I))
 230        CONTINUE
            H0=DATAXT(0)
            H1=DCONJG(DATAXT(0))
            FPT(NX,J,K)=0.5D0*(H0+H1)+0.5D0*CI*(H0-H1)
            DO 235 I=0,NX
               FPT(I,J,K)=FPT(I,J,K)/DBLE(NX*2)
 235        CONTINUE
 210     CONTINUE
c z direction
         IF (NZ.GT.1) THEN
         DO 240 I=0,NX
            DO 250 K=0,NZ-1
               DATAZ(K)=FPT(I,J,K)
 250        CONTINUE
            CALL FFT(DATAZ,NZ,-1)            
            DO 260 K=0,NZ-1
               FPT(I,J,K)=DATAZ(K)/DBLE(NZ)
 260        CONTINUE
 240     CONTINUE
         ENDIF
         DO 300 I=0,NX
           FPT(I,J,NZ)=FPT(I,J,0)
 300     CONTINUE
 200  CONTINUE
      RETURN
      END

      SUBROUTINE S2P(FPT,FP)
      INCLUDE '../par.f'
      COMMON /F1F/ WCX
      COMMON /F2/ IFFLG
      REAL*8 DATAX(0:NX*2-1),FP(0:NX*2+1,0:NY+1,0:NZ+1)
      COMPLEX*16 DATAXT(0:NX-1), DATAZ(0:NZ-1)
      COMPLEX*16 FPT(0:NX,0:NY+1,0:NZ+1), CI, H0, H1
      COMPLEX*16 WCX(0:NX*2+1)
      EQUIVALENCE(DATAX,DATAXT)
      CI=DCMPLX(0.D0,1.D0)
      IF(IFFLG.EQ.0)THEN
         WRITE(*,*)'Making trigonometric table'
         DO 100 I=0,NX*2-1
            WCX(I)=DCMPLX(DCOS(2.D0*PI*DBLE(I)/DBLE(NX*2)),
     +           DSIN(2.D0*PI*DBLE(I)/DBLE(NX*2)))
 100     CONTINUE
         IFFLG=1
      ENDIF
c
      DO 200 J=1,NY
c z direction
         If (NZ.GT.1) THEN
         DO 240 I=0,NX
            DO 250 K=0,NZ-1
               DATAZ(K)=FPT(I,J,K)
 250        CONTINUE
            CALL FFT(DATAZ,NZ,1)            
            DO 260 K=0,NZ-1
               FPT(I,J,K)=DATAZ(K)
 260        CONTINUE
 240     CONTINUE
         ENDIF
c symmetry
         DO 210 K=0,NZ-1
            DO 230 I=0,NX-1
               H0=0.5D0*(FPT(I,J,K)+DCONJG(FPT(NX-I,J,K)))
               H1=0.5D0*(FPT(I,J,K)-DCONJG(FPT(NX-I,J,K)))*WCX(I)
               DATAXT(I)=H0+CI*H1
 230        CONTINUE
c theta direction
            CALL FFT(DATAX,NX,1)
            DO 220 I=0,NX*2-1
               FP(I,J,K)=DATAX(I)*2.D0
 220        CONTINUE
 210     CONTINUE
c b.c.
         DO 300 K=0,NZ
            FP(NX*2,J,K)=FP(0,J,K)
 300     CONTINUE
         DO 310 I=0,NX*2
            FP(I,J,NZ)=FP(I,J,0)
 310     CONTINUE
 200  CONTINUE
      RETURN
      END

