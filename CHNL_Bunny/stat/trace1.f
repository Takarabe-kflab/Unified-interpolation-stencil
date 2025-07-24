c
c stat/trace1.f
c (C) 2011 N.Tomiyama 
c
      SUBROUTINE TRACE1()
      INCLUDE '../par.f'
      INCLUDE '../common.f'
      CHARACTER*6 EXT
      CALL UBCWING()
c
         WRITE(11,610) TID,
c     +   CL, CD, FXIBM, FYIBM, 
     +   UX(16,48,1), UX(48,48,1), UX(96,48,1), UX(128,48,1)
c
 610  FORMAT(' ',10E15.7)
c
      RETURN
      END
      
