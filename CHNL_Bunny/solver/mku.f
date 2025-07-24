c
c solver/mku.f
c (C) 2004 K. Fukagata 
c
      SUBROUTINE MKU(GAMMA,ZETA,ALPHA)
      INCLUDE '../par.f'
      INCLUDE '../common.f'
c------------------------------------------------------------------------------
c 1. Explicit step
c   (DUX in LHS = delta Ux; DUX in RHS = Diffusion term)
c------------------------------------------------------------------------------
      DO 100 K=1,NZ
      DO 100 J=1,NY         
      DO 100 I=1,NX
         DUX(I,J,K)=DT*(GAMMA*FX(I,J,K)
     +                  +ZETA*FX1(I,J,K)
     +                  +ALPHA*(-GXP(I,J,K)+DUX(I,J,K))
     +                 )
         DUY(I,J,K)=DT*(GAMMA*FY(I,J,K)
     +                  +ZETA*FY1(I,J,K)
     +                  +ALPHA*(-GYP(I,J,K)+DUY(I,J,K))
     +                 )
         DUZ(I,J,K)=DT*(GAMMA*FZ(I,J,K)
     +                 +ZETA*FZ1(I,J,K)
     +                 +ALPHA*(-GZP(I,J,K)+DUZ(I,J,K)))
         DTT(I,J,K)=DT*(GAMMA*FT(I,J,K)
     +                 +ZETA*FT1(I,J,K)
     +                 +ALPHA*DTT(I,J,K))
c RHS at previous substep
         FX1(I,J,K)=FX(I,J,K)
         FY1(I,J,K)=FY(I,J,K)
         FZ1(I,J,K)=FZ(I,J,K)
         FT1(I,J,K)=FT(I,J,K)
 100  CONTINUE
c------------------------------------------------------------------------------
c 2. Implicit step (Approximate factrization)
c------------------------------------------------------------------------------
      CALL AFY(ALPHA)
c------------------------------------------------------------------------------
c 3. Update U
c------------------------------------------------------------------------------
      DO 300 K=1,NZ
      DO 300 J=1,NY         
      DO 300 I=1,NX
         UX(I,J,K)=UX(I,J,K)+DUX(I,J,K)
         UY(I,J,K)=UY(I,J,K)+DUY(I,J,K)
         UZ(I,J,K)=UZ(I,J,K)+DUZ(I,J,K)
         TT(I,J,K)=TT(I,J,K)+DTT(I,J,K)
 300  CONTINUE
c------------------------------------------------------------------------------
c 4. Preserve Diffusion term
c------------------------------------------------------------------------------
      DO 400 K=1,NZ
      DO 400 J=1,NY         
      DO 400 I=1,NX
c        DDUX(I,J,K)=DUX(I,J,K)
c        DDUY(I,J,K)=DUY(I,J,K)
 400  CONTINUE
      RETURN
      END

