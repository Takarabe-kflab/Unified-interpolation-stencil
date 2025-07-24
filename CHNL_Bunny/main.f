c
c main.f
c (C) 2004 K. Fukagata 
c
      PROGRAM CHNL
      INCLUDE 'par.f'
      INCLUDE 'common.f'
      REAL*8 CG(3), CZ(3), CA(3)
c Coefficients for RK3/CN scheme
      CG(1)=8.D0/15.D0
      CG(2)=5.D0/12.D0
      CG(3)=3.D0/4.D0
      CZ(1)=0.D0
      CZ(2)=-17.D0/60.D0
      CZ(3)=-5.D0/12.D0
      CA(1)=2.D0*4.D0/15.D0
      CA(2)=2.D0*1.D0/15.D0
      CA(3)=2.D0*1.D0/6.D0
c
c 1. Initialization
c
      CALL GRID()
      CALL INIT()
      CALL KER()
      CALL UBC()
      CALL GHP()
c
c 2. Main loop
c
      OPEN(11,FILE='trace.out',STATUS='UNKNOWN')
c
      DO 10000 III=1,NNN
         TID=TID+DT
         write(*,*)'---------------------------------------------------'
         write(*,*)'main: Time step =',III
         write(*,*)'main: Time =',TID
         write(*,*)'---------------------------------------------------'
c Compute control input
         CALL OUTLET()
         CALL UBC()
         CALL GHP()
c
c RK3/CN integtaion
C------------------------------------------------------------------------------
         DO 1000 L=1,3
            write(*,*)'main: Substep',L
c Provisional velocity
            CALL MKF()
            CALL MKU(CG(L),CZ(L),CA(L))
            CALL UBC()
            CALL GHP()
c Velocity correction
            CALL MKP(CA(L)) 
            CALL UBC()
            CALL GHP()
 1000    CONTINUE
C------------------------------------------------------------------------------
c Check
         CALL CHKCONT()
         CALL CHKCFL()
c Output
         IF(MOD(III,10).EQ.0)THEN
            CALL STAT0()
         ENDIF
         IF(MOD(III,IIC).EQ.0)THEN
            CALL STAT1()
         ENDIF
         IF(MOD(III,IIO).EQ.0) CALL RAWOUT()
         IF(MOD(III,240).EQ.0) CALL OUTMAT()
10000 CONTINUE
      CLOSE(11)

      END



