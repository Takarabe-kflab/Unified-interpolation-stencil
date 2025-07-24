      SUBROUTINE FFT(A,N,ISIGN)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 A(2,0:*)
      J=0
c permutation matrix
      DO 20 I=0,(N-2)
         IF(I.LT.J)THEN
            DO 25 L=1,2
               T=A(L,J)
               A(L,J)=A(L,I)
               A(L,I)=T
 25            CONTINUE
         ENDIF
         K=N/2
 10      IF(K.LE.J)THEN
            J=J-K
            K=K/2
            GOTO 10
         ENDIF
         J=J+K
 20   CONTINUE
c log2N matrix-vector multiplications
      S=0.D0
      C=-1.0D0
      L=1
 30   LH=L
      L=L+L
      UR=1.0D0
      UI=0.0D0
      DO 50 J=0,LH-1
         DO 40 I=J,(N-1),L
            IP=I+LH
            TR=A(1,IP)*UR-A(2,IP)*UI
            TI=A(1,IP)*UI+A(2,IP)*UR
            A(1,IP)=A(1,I)-TR
            A(2,IP)=A(2,I)-TI
            A(1,I)=A(1,I)+TR
            A(2,I)=A(2,I)+TI
 40      CONTINUE
         TI=UR*S+UI*C
         UR=UR*C-UI*S
         UI=TI
 50   CONTINUE
      S=SQRT(0.5D0*(1.0D0-C))*ISIGN
      C=SQRT(0.5D0*(1.0D0+C))
      IF(L.LT.N)GOTO 30
      RETURN
      END
