c
c common.f
c (C) 2004-2024 K. Fukagata 
c
c
c 0. Numerical parameters
c      
      COMMON /A1/ TID, DT
      COMMON /A2/ NNN, INS, NST, IIO, IIC, III
c 
c------------------------------------------------------------------------------
c TID: Time counter
c DT : Time step [input]
c NNN: Number of time steps to be computed [input]
c INS: Flag indicating whether the statistics computation is continued [input] 
c      (0: Initialize. 1: Read "stat" and continue)
c NST: Counter for statistics computation
c IIO: Output field data every IIO step [input]
c IIC: Output "stat" every IIC step [input]
c------------------------------------------------------------------------------
c
c 1. Flow Parameters
c
      COMMON /C11/ REB, PRL, DTDY
c
c 2. Coordinate
c
      COMMON /C21/ X, X12, DX 
      COMMON /C22/ Y, Y12, DY, DY12, YY12
      COMMON /C23/ Z, Z12, DZ
      REAL*8 X(0:NX+2), X12(0:NX+2), DX
      REAL*8 Y(0:NY+2), Y12(0:NY+2), DY(0:NY+2),
     +       DY12(0:NY+2), YY12(0:NY+2)
      REAL*8 Z(0:NZ+2), Z12(0:NZ+2), DZ
c
c 3. Velocity and B.C.
c
      COMMON /C31/ UX, VXI, VXO
      COMMON /C32/ UY, VYI, VYO
      COMMON /C33/ UZ, VZI, VZO
      COMMON /C34/ TT, TTI, TTO
      COMMON /C35/ VXM, VZM
      REAL*8 UX(-2:NX+3,0:NY+2,-2:NZ+3),
     +     VXI(0:NY+2,0:NZ+2), VXO(0:NY+2,0:NZ+2)
      REAL*8 UY(-2:NX+3,0:NY+2,-2:NZ+3),
     +     VYI(0:NY+2,0:NZ+2), VYO(0:NY+2,0:NZ+2)
      REAL*8 UZ(-2:NX+3,0:NY+2,-2:NZ+3), 
     +     VZI(0:NY+2,0:NZ+2), VZO(0:NY+2,0:NZ+2)
      REAL*8 TT(-2:NX+3,0:NY+2,-2:NZ+3), 
     +     TTI(0:NY+2,0:NZ+2), TTO(0:NY+2,0:NZ+2)
c 
c 4. Advection term
c
      COMMON /C41/ FX, FX1
      COMMON /C42/ FY, FY1
      COMMON /C43/ FZ, FZ1
      COMMON /C44/ FT, FT1
      REAL*8 FX(0:NX+2,0:NY+2,0:NZ+2), FX1(0:NX+2,0:NY+2,0:NZ+2)
      REAL*8 FY(0:NX+2,0:NY+2,0:NZ+2), FY1(0:NX+2,0:NY+2,0:NZ+2)
      REAL*8 FZ(0:NX+2,0:NY+2,0:NZ+2), FZ1(0:NX+2,0:NY+2,0:NZ+2)
      REAL*8 FT(0:NX+2,0:NY+2,0:NZ+2), FT1(0:NX+2,0:NY+2,0:NZ+2)
c
c 5. Delta U
c
      COMMON /C51/ DUX, DUY, DUZ, DTT
c, DUDY
      REAL*8 DUX(0:NX+2,0:NY+2,0:NZ+2) 
      REAL*8 DUY(0:NX+2,0:NY+2,0:NZ+2) 
      REAL*8 DUZ(0:NX+2,0:NY+2,0:NZ+2) 
      REAL*8 DTT(0:NX+2,0:NY+2,0:NZ+2)
c
c 6. Pressure
c
      COMMON /C61/ P, GXP, GYP, GZP
      REAL*8   P(-2:NX+3,0:NY+2,-2:NZ+3)
      REAL*8 GXP(0:NX+2,0:NY+2,0:NZ+2)
      REAL*8 GYP(0:NX+2,0:NY+2,0:NZ+2)
      REAL*8 GZP(0:NX+2,0:NY+2,0:NZ+2)
c
c 7. Statistics
c
      COMMON /C71/ ST
      REAL*8 ST(0:NY+2,MST)
c
c 8. FFT
c
      COMMON /C81/ KX,KZ,KXZ
      REAL*8 KX(0:NX), KZ(0:NZ), KXZ(0:NX,0:NZ)
c
c 9. TDMA
c
      COMMON /C91/ AC, AE
      REAL*8 AC(1:NY,3), AE(1:NY,3)
c
c 10. IBM
c
c     GPU, GPV, GPW: table with the convolution kernel for U, V and W
c     UGP: velocity of one ghost point 
c     IGPU, IGPV, IGPW: I of the GP for U, V and W
c     JGPU, JGPV, JGPW: J of the GP for U, V and W
c     NGPU, NGPV, NGPW: number of ghost points for U, V and W
c     NIP: number of inside0 (inside but not ghost) points of the shape
c     INSP: table with coordinates of the inside0 points 
      
      COMMON /C101/ GPU, GPV, GPW, UGP
      REAL*8 GPU(1:NGPMAX,126), GPV(1:NGPMAX,126), GPW(1:NGPMAX,126), UGP
      COMMON /C102/ IGPU, IGPV, IGPW, JGPU, JGPV, JGPW, KGPU, KGPV, KGPW
      INTEGER IGPU(1:NGPMAX), IGPV(1:NGPMAX), IGPW(1:NGPMAX)
      INTEGER JGPU(1:NGPMAX), JGPV(1:NGPMAX), JGPW(1:NGPMAX)
      INTEGER KGPU(1:NGPMAX), KGPV(1:NGPMAX), KGPW(1:NGPMAX)
      COMMON /C103/ NGPU, NGPV, NGPW
      INTEGER NGPU, NGPV, NGPW
c      INTEGER INSP(1:NIPMAX,2), NIP
