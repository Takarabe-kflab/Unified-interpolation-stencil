c
c stat/stat1.f
c (C) 2004-2006 K. Fukagata
c 2009-05-13 Bug fixed by K. Fukagata 
c
      SUBROUTINE STAT1()
      INCLUDE '../par.f'
      INCLUDE '../common.f'
      CALL UBC()
c
      NST=NST+1
      write(*,*)'STAT:',NST 
c
c 1st order moments
c
c 1: <ux> i+1/2, j, k
c 2: <uy> i, j+1/2, k 
c 3: <uz> i, j, k+1/2
c 4: <p>  i, j, k
c 5: <wx> i, j+1/2, k+1/2
c 6: <wy> i+1/2, j, k+1/2
c 7: <wz> i+1/2, j+1/2, k
c 8: <T> i,j,k
      DO 1000 J=0,NY
         DO 1010 IS=1,8
            ST(J,IS)=ST(J,IS)*DBLE((NST-1)*NX*NZ)
 1010    CONTINUE
         DO 1020 K=1,NZ
         DO 1020 I=1,NX
            ST(J,1)=ST(J,1)+UX(I,J,K)
            ST(J,2)=ST(J,2)+UY(I,J,K)
            ST(J,3)=ST(J,3)+UZ(I,J,K)
            ST(J,4)=ST(J,4)+P (I,J,K)
            ST(J,5)=ST(J,5)
     +           +(UZ(I,J+1,K)-UZ(I,J,K))/DY12(J)
     +           -(UY(I,J,K+1)-UY(I,J,K))/DZ
            ST(J,6)=ST(J,6)
     +           +(UX(I,J,K+1)-UX(I,J,K))/DZ
     +           -(UZ(I+1,J,K)-UZ(I,J,K))/DX
            ST(J,7)=ST(J,7)
     +           +(UY(I+1,J,K)-UY(I,J,K))/DX
     +           -(UX(I,J+1,K)-UX(I,J,K))/DY12(J)
            ST(J,8)=ST(J,8)+TT(I,J,K)
 1020    CONTINUE
         DO 1030 IS=1,8
            ST(J,IS)=ST(J,IS)/DBLE(NST*NX*NZ)
 1030    CONTINUE
 1000 CONTINUE
c
c higher order moments 
c
      DO 2000 J=0,NY
         DO 2010 IS=9,22
            ST(J,IS)=ST(J,IS)*DBLE((NST-1)*NX*NZ)
 2010    CONTINUE
         DO 2020 K=1,NZ
         DO 2020 I=1,NX
c u'
            UXR=UX(I,J,K)-ST(J,1)
            UYR=UY(I,J,K)-ST(J,2)
            UZR=UZ(I,J,K)-ST(J,3)
            PR = P(I,J,K)-ST(J,4)
            WXR=(UZ(I,J+1,K)-UZ(I,J,K))/DY12(J)
     +         -(UY(I,J,K+1)-UY(I,J,K))/DZ
     +         -ST(J,5)
            WYR=(UX(I,J,K+1)-UX(I,J,K))/DZ
     +         -(UZ(I+1,J,K)-UZ(I,J,K))/DX
     +         -ST(J,6)
            WZR=(UY(I+1,J,K)-UY(I,J,K))/DX
     +         -(UX(I,J+1,K)-UX(I,J,K))/DY12(J)
     +         -ST(J,7)
            TR=TT(I,J,K)-ST(J,8)
c    
c 2nd order moments
c
c 9: <ux'^2>
c 10: <uy'^2>
c 11: <uz'^2>
c 12: <p'^2>
c 13: <wx'^2>
c 14: <wy'^2>
c 15: <wz'^2>
c 16: <ux uy> i+1/2, j+1/2, k
c 17: <uy uz> i, j+1/2, k+1/2
c 18: <uz ux> i+1/2, j, k+1/2
c 19: <T'^2>
c 20: <ux'T'> i+1/2,j,k
c 21: <uy'T'> i,j+1/2,k
c 22: <uz'T'> i,j,k+1/2
            ST(J,9 )=ST(J,9 )+UXR**2
            ST(J,10)=ST(J,10)+UYR**2
            ST(J,11)=ST(J,11)+UZR**2           
            ST(J,12)=ST(J,12)+PR**2
            ST(J,13)=ST(J,13)+WXR**2
            ST(J,14)=ST(J,14)+WYR**2
            ST(J,15)=ST(J,15)+WZR**2           
            ST(J,16)=ST(J,16)
     +           +((UX(I,J,K)-ST(J,1)+UX(I,J+1,K)-ST(J+1,1))/2.D0)
     +           *((UY(I,J,K)+UY(I+1,J,K))/2.D0-ST(J,2))
            ST(J,17)=ST(J,17)
     +           +((UY(I,J,K)+UY(I,J,K+1))/2.D0-ST(J,2))
     +           *((UZ(I,J,K)-ST(J,3)+UZ(I,J+1,K)-ST(J+1,3))/2.D0)
            ST(J,18)=ST(J,18)
     +           +((UZ(I,J,K)+UZ(I+1,J,K))/2.D0-ST(J,3))
     +           *((UX(I,J,K)+UX(I,J,K+1))/2.D0-ST(J,1))
            ST(J,19)=ST(J,19)+TR**2
            ST(J,20)=ST(J,20)
     +           +((TT(I,J,K)+TT(I+1,J,K))/2.D0-ST(J,8))
     +           *(UX(I,J,K)-ST(J,1))
            ST(J,21)=ST(J,21)
     +           +(TT(I,J,K)-ST(J,8)+TT(I,J+1,K)-ST(J+1,8))/2.D0
     +           *(UY(I,J,K)-ST(J,2))
            ST(J,22)=ST(J,22)
     +           +((TT(I,J,K)+TT(I,J,K+1))/2.D0-ST(J,8))
     +           *(UZ(I,J,K)-ST(J,3))
 2020    CONTINUE
         DO 2030 IS=9,22
            ST(J,IS)=ST(J,IS)/DBLE(NST*NX*NZ)
 2030     CONTINUE
 2000 CONTINUE
c
c Output
c
      OPEN(10,FILE='ctr.out',STATUS='UNKNOWN')
      WRITE(10,*)NST
      CLOSE(10)
      OPEN(10,FILE='stat.out',STATUS='UNKNOWN')
      DO 9000 J=0,NY
c     WRITE(10,600)Y(J),Y12(J),(ST(J,IS),IS=1,MST)
      WRITE(10,600)Y(J),Y12(J),ST(J,1)
         
 600     FORMAT(' ',77E15.7)
 9000 CONTINUE
      CLOSE(10)
      RETURN
      END
