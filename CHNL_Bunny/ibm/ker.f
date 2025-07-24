c ibm/ker.f
c 2024 I. Takarabe
c
	SUBROUTINE KER()
	INCLUDE 'par.f'
	INCLUDE 'common.f'
	WRITE(*,*)'*********************************************'
	WRITE(*,*)'*       Reading convolution kernels         *'

c
	OPEN(28,FILE='kernel/GP_weights_u.txt',STATUS='OLD')
	READ(28,*)NGPU
	DO 100 I=1,NGPU
	   READ(28,*)IGPU(I),JGPU(I),KGPU(I),(GPU(I,J),J=1,126)
 100	CONTINUE
	CLOSE(28)
	WRITE(*,*)'*** GP_weights_u successfully imported  ******'
	
	OPEN(28,FILE='kernel/GP_weights_v.txt',STATUS='OLD')
	READ(28,*)NGPV
	DO 200 I=1,NGPV
	   READ(28,*)IGPV(I),JGPV(I),KGPV(I),(GPV(I,J),J=1,126)
 200	CONTINUE
 	CLOSE(28)
	WRITE(*,*)'*** GP_weights_v successfully imported  ******'

	OPEN(28,FILE='kernel/GP_weights_w.txt',STATUS='OLD')
	READ(28,*)NGPW
	DO 300 I=1,NGPW
	   READ(28,*)IGPW(I),JGPW(I),KGPW(I),(GPW(I,J),J=1,126)
 300	CONTINUE
 	CLOSE(28)
	WRITE(*,*)'*** GP_weights_w successfully imported  ******'
	WRITE(*,*)'*********************************************'
	RETURN
	END
