c solver/outlet.f
c 2007/7/20 shimada
c 2009-05-12 K. Fukagata

      SUBROUTINE OUTLET()
      INCLUDE '../par.f'
      INCLUDE '../common.f'
c
c outlet 
c
      DO 700 K=1,NZ
      DO 700 J=1,NY
         VXO(J,K)=VXO(J,K)-DT/DX*1.D0*(VXO(J,K)-UX(NX-1,J,K))
         VYO(J,K)=VYO(J,K)-2.D0*DT/DX*1.D0*(VYO(J,K)-UY(NX,J,K))
         VZO(J,K)=VZO(J,K)-2.D0*DT/DX*1.D0*(VZO(J,K)-UZ(NX,J,K))
         TTO(J,K)=TTO(J,K)-2.D0*DT/DX*1.D0*(TTO(J,K)-TT(NX,J,K))
 700  CONTINUE
c
c Flux correction
c
c inlet
      QI=0.D0
      DO 131 K=1,NZ
      DO 131 J=1,NY
         QI=QI+DZ*DY(J)*VXI(J,K)
 131  CONTINUE
c outlet
      QO=0.D0
      DO 132 K=1,NZ
      DO 132 J=1,NY
         QO=QO+DZ*DY(J)*VXO(J,K)
 132  CONTINUE
      QERR=QI-QO
      write(*,*)'ubc: Qi, Qo, Qerr = ',QI, QO, QERR
c correction at upper/lower boundaries
      DO K=1,NZ
      DO J=1,NY
         VXO(J,K)=VXO(J,K)*(QI+QU-QL)/QO
      END DO
      END DO
      RETURN
      END
