c
c par.f
c (C) 2004 K. Fukagata 
c
      IMPLICIT REAL*8(A-H,O-Z)
c
c Note: NX=2^(NFXX), NZ=2^(NFXZ)
c
      PARAMETER(PI=3.14159265358979312D0)
      PARAMETER(NX=512, NY=256, NZ=256, NMAX=512)
      PARAMETER(XL=8.0D0, YL=4.0D0, ZL=4.0D0)
      PARAMETER(NVXE0=3333)
      PARAMETER(MST=22)
c IBM:
c NGPMAX and NIPMAX are used to initialize GP and INS tables
c UOB and VOP are the velocity of the object in x and y direction
      PARAMETER(NGPMAX=100000)
      PARAMETER(NIPMAX=100000)
      PARAMETER(UOB=0.D0)
      PARAMETER(VOB=0.D0)
