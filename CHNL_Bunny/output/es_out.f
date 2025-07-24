c
c output/ga_out.f
c (C)2011 K.Nakayama
c
      SUBROUTINE ES_OUT()
      INCLUDE '../par.f'
      INCLUDE '../common.f'
     
      OPEN(10,FILE='output_test.out',STATUS='UNKNOWN')

      WRITE(10,600)CDM
 600  FORMAT('',1E15.7)
      CLOSE(10)

      RETURN
      END



