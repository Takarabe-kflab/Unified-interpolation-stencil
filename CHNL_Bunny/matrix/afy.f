c
c matrix/afy.f
c (C) 2004 K. Fukagata 
c
      SUBROUTINE AFY(ALPHA)
      INCLUDE '../par.f'
      INCLUDE '../common.f'
      REAL*8 A12(NY-1,3),BB12(NX,NY-1)
      REAL*8 A(NY,3),BB(NX,NY)
      BETA=0.5D0*ALPHA*DT/REB    
c
c Ux
c
*poption parallel,tlocal(BB,A,J)
      DO 300 K=1,NZ          
c Set matrix coefficient
         DO 310 J=1,NY
            A(J,1)=    -BETA*AC(J,1)
            A(J,2)=1.D0-BETA*AC(J,2)
            A(J,3)=    -BETA*AC(J,3)
 310     CONTINUE
         DO 315 J=1,NY
         DO 315 I=1,NX
            BB(I,J)=DUX(I,J,K)
 315     CONTINUE
c Dirichlet B.C.
         A(1, 2)=A(1,2)-A(1,1)
         A(1, 1)=0.D0        
         A(NY,2)=A(NY,2)-A(NY,3)
         A(NY,3)=0.D0
c Solve matrix
         CALL TRID2(A,BB,NX,NY)
         DO 320 J=1,NY
         DO 320 I=1,NX
            DUX(I,J,K)=BB(I,J)
 320     CONTINUE
 300  CONTINUE
c
c Uy
c
*poption parallel,tlocal(BB12,A12,J)
      DO 100 K=1,NZ 
c Set matrix coefficient                    
         DO 110 J=1,NY-1
            A12(J,1)=    -BETA*AE(J,1)
            A12(J,2)=1.D0-BETA*AE(J,2)
            A12(J,3)=    -BETA*AE(J,3)
 110     CONTINUE
         DO 115 J=1,NY-1
         DO 115 I=1,NX
            BB12(I,J)=DUY(I,J,K)
 115     CONTINUE
c Dirichlet B.C. (delta Uy is zero)
         A12(1,   1)=0.D0
         A12(NY-1,3)=0.D0
c Solve matrix
         CALL TRID2(A12,BB12,NX,NY-1)
         DO 120 J=1,NY-1
         DO 120 I=1,NX
            DUY(I,J,K)=BB12(I,J)
 120     CONTINUE
         DO 125 I=1,NX
            DUY(I,NY,K)=0.D0
 125     CONTINUE
 100  CONTINUE
c
c Uz
c
*poption parallel,tlocal(BB,A,J)
      DO 200 K=1,NZ 
c Set Matrix coefficient                    
         DO 210 J=1,NY
            A(J,1)=    -BETA*AC(J,1)
            A(J,2)=1.D0-BETA*AC(J,2)
            A(J,3)=    -BETA*AC(J,3)
 210     CONTINUE
         DO 215 J=1,NY
         DO 215 I=1,NX
            BB(I,J)=DUZ(I,J,K)
 215     CONTINUE
c Dirichlet B.C.        
         A(1,2)=A(1,2)-A(1,1)
         A(1,1)=0.D0        
         A(NY,2)=A(NY,2)-A(NY,3)
         A(NY,3)=0.D0
c Solve matrix
         CALL TRID2(A,BB,NX,NY)
         DO 220 J=1,NY
         DO 220 I=1,NX
            DUZ(I,J,K)=BB(I,J)
 220     CONTINUE
 200  CONTINUE
c
c TT
c
*poption parallel,tlocal(BB,A,J)
      DO 400 K=1,NZ 
c Set Matrix coefficient                    
         DO 410 J=1,NY
            A(J,1)=    -BETA*AC(J,1)/PRL
            A(J,2)=1.D0-BETA*AC(J,2)/PRL
            A(J,3)=    -BETA*AC(J,3)/PRL
 410     CONTINUE
         DO 415 J=1,NY
         DO 415 I=1,NX
            BB(I,J)=DTT(I,J,K)
 415     CONTINUE
c Dirichlet B.C.        
         A(1,2)=A(1,2)-A(1,1)
         A(1,1)=0.D0        
         A(NY,2)=A(NY,2)-A(NY,3)
         A(NY,3)=0.D0
c Solve matrix
         CALL TRID2(A,BB,NX,NY)
         DO 420 J=1,NY
         DO 420 I=1,NX
            DTT(I,J,K)=BB(I,J)
 420     CONTINUE
 400  CONTINUE
      RETURN
      END
