c
c ctrl/ctrl0.f (dummy)
c (C) 2004 K. Fukagata 
c
      SUBROUTINE CTRL()
      INCLUDE '../par.f'
      INCLUDE '../common.f'
      COMMON /RAND/ AMP
c
c Perturbation --- needs to be modified!!!
c
      AMP=0.1D0
      OM=4.D0/3.D0
c
      VXM=0.D0
      VXRMS=0.D0
      DO 100 K=1,NZ
      DO 100 J=1,NY
        VXTMP=1.D0
        VXI(J,K)=VXTMP
        VXM=VXM+VXTMP
        VXRMS=VXRMS+VXTMP**2
 100  CONTINUE
      VXM=VXM/DBLE(NY*NZ)
      VXRMS=DSQRT(VXRMS/DBLE(NY*NZ)-VXM**2)
c
      WRITE(*,*)'ctrl0: U_i^f mean, rms: ', VXM, VXRMS
c
      RETURN
      END
