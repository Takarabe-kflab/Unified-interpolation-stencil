c
c output/es_in.f
c (C)2011 K.Nakayama
c
      SUBROUTINE ES_IN()
      INCLUDE '../par.f'
      INCLUDE '../common.f'

      OPEN(10,FILE='input1.dat',STATUS='OLD')
      READ(10,*)A1
      CLOSE(10)
      
      OPEN(20,FILE='input2.dat',STATUS='OLD')
      READ(20,*)A2
      CLOSE(20)

      RETURN
      END




