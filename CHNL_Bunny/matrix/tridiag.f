c
c matrix/tridiag.f
c (C) 2004 K. Fukagata 
c
      SUBROUTINE TRIDIAG(A,B,NN)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 A(NN,3), B(NN)
      DO 20 I=1,NN-1
         A(I,3)=A(I,3)/A(I,2)
         B(I)=B(I)/A(I,2)
         A(I+1,2)=A(I+1,2)-A(I+1,1)*A(I,3)         
         B(I+1)=B(I+1)-A(I+1,1)*B(I)
 20   CONTINUE
      B(NN)=B(NN)/A(NN,2)
      DO 30 I=NN-1,1,-1
         B(I)=B(I)-A(I,3)*B(I+1)
 30   CONTINUE
      RETURN
      END

      SUBROUTINE TRID2(A,B,M,N)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 A(N,3), B(M,N)
      DO 20 I=1,N-1
         A(I,3)=A(I,3)/A(I,2)
         A(I+1,2)=A(I+1,2)-A(I+1,1)*A(I,3)         
         DO 25 J=1,M
            B(J,I)=B(J,I)/A(I,2)
            B(J,I+1)=B(J,I+1)-A(I+1,1)*B(J,I)
 25      CONTINUE
 20   CONTINUE
      DO 27 J=1,M
         B(J,N)=B(J,N)/A(N,2)
 27   CONTINUE
      DO 30 I=N-1,1,-1
         DO 35 J=1,M
            B(J,I)=B(J,I)-A(I,3)*B(J,I+1)
 35      CONTINUE
 30   CONTINUE
      RETURN
      END

      SUBROUTINE CTRIDIAG(A,B,NN)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 A(NN,3)
      COMPLEX*16 B(NN)
      DO 20 I=1,NN-1
         A(I,3)=A(I,3)/A(I,2)
         B(I)=B(I)/A(I,2)
         A(I+1,2)=A(I+1,2)-A(I+1,1)*A(I,3)         
         B(I+1)=B(I+1)-A(I+1,1)*B(I)
 20   CONTINUE
      B(NN)=B(NN)/A(NN,2)
      DO 30 I=NN-1,1,-1
         B(I)=B(I)-A(I,3)*B(I+1)
 30   CONTINUE
      RETURN
      END

      SUBROUTINE CYCLIC(A,ABOT,ATOP,B,N)
      INCLUDE '../par.f'
      REAL*8 A(N,3),B(N)
      REAL*8 U(NMAX)
      GAMMA=-A(1,2)
      DO 100 I=1,N
         U(I)=0.         
 100  CONTINUE
      U(1)=GAMMA
      U(N)=ABOT
      A(1,2)=A(1,2)-GAMMA
      A(N,2)=A(N,2)-ABOT*ATOP/GAMMA  
c 2 tridiag
      DO 110 I=1,N-1
         A(I,3)=A(I,3)/A(I,2)
         B(I)=B(I)/A(I,2)
         U(I)=U(I)/A(I,2)
         A(I+1,2)=A(I+1,2)-A(I+1,1)*A(I,3)         
         B(I+1)=B(I+1)-A(I+1,1)*B(I)
         U(I+1)=U(I+1)-A(I+1,1)*U(I)
 110  CONTINUE
      B(N)=B(N)/A(N,2)
      U(N)=U(N)/A(N,2)
      DO 120 I=N-1,1,-1
         B(I)=B(I)-A(I,3)*B(I+1)
         U(I)=U(I)-A(I,3)*U(I+1)
 120  CONTINUE
c
      FACT=(B(1)+ATOP*B(N)/GAMMA)/(1.D0+U(1)+ATOP*U(N)/GAMMA)
      DO 130 I=1,N
         B(I)=B(I)-FACT*U(I)
 130  CONTINUE
      RETURN
      END

      SUBROUTINE CYCL2(A,ABOT,ATOP,B,M,N)
      INCLUDE '../par.f'
      REAL*8 A(N,3),B(M,N)
      REAL*8 U(NMAX),FACT(NMAX)
      GAMMA=-A(1,2)
      DO 100 I=1,N
         U(I)=0.         
 100  CONTINUE
      U(1)=GAMMA
      U(N)=ABOT
      A(1,2)=A(1,2)-GAMMA
      A(N,2)=A(N,2)-ABOT*ATOP/GAMMA  
c 2 tridiag
      DO 110 I=1,N-1
         A(I,3)=A(I,3)/A(I,2)
         A(I+1,2)=A(I+1,2)-A(I+1,1)*A(I,3)         
         DO 115 J=1,M
            B(J,I)=B(J,I)/A(I,2)
            B(J,I+1)=B(J,I+1)-A(I+1,1)*B(J,I)
 115     CONTINUE
         U(I)=U(I)/A(I,2)
         U(I+1)=U(I+1)-A(I+1,1)*U(I)
 110  CONTINUE
      DO 117 J=1,M
         B(J,N)=B(J,N)/A(N,2)
 117  CONTINUE
      U(N)=U(N)/A(N,2)
      DO 120 I=N-1,1,-1
         DO 125 J=1,M
            B(J,I)=B(J,I)-A(I,3)*B(J,I+1)
 125     CONTINUE
         U(I)=U(I)-A(I,3)*U(I+1)
 120  CONTINUE
c
      DO 127 J=1,M
         FACT(J)=(B(J,1)+ATOP*B(J,N)/GAMMA)/(1.D0+U(1)+ATOP*U(N)/GAMMA)
 127  CONTINUE
      DO 130 I=1,N
      DO 130 J=1,M
         B(J,I)=B(J,I)-FACT(J)*U(I)
 130  CONTINUE
      RETURN
      END

