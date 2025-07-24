c ibm/ghp.f
c 2024 I. Takarabe

	SUBROUTINE GHP()
	INCLUDE 'par.f'
	INCLUDE 'common.f'

	INTEGER I,J,K

	DO 140 N=1,NGPU
    	I=IGPU(N)
	   	J=JGPU(N)
		K=KGPU(N)
	   	UGP=0.D0
	   	DO 150 A=1,5
	   		DO 160 B=1,5
				DO 170 C=1,5
	      			UGP=UGP+UX(I+A-3,J+B-3,K+C-3)*GPU(N,25*(A-1)+5*(B-1)+C)
 170				CONTINUE
 160			CONTINUE
 150 		CONTINUE
	   	UX(I,J,K)=UGP+UOB*GPU(N,126)
 140	CONTINUE

	DO 240 N=1,NGPV
    	I=IGPV(N)
	   	J=JGPV(N)
		K=KGPV(N)
	   	UGP=0.D0
	   	DO 250 A=1,5
	   		DO 260 B=1,5
				DO 270 C=1,5
	      			UGP=UGP+UY(I+A-3,J+B-3,K+C-3)*GPV(N,25*(A-1)+5*(B-1)+C)
 270				CONTINUE
 260			CONTINUE
 250 		CONTINUE
	   	UY(I,J,K)=UGP+VOB*GPV(N,126)
 240	CONTINUE

	DO 340 N=1,NGPW
    	I=IGPW(N)
	   	J=JGPW(N)
		K=KGPW(N)
	   	UGP=0.D0
	   	DO 350 A=1,5
	   		DO 360 B=1,5
				DO 370 C=1,5
	      			UGP=UGP+UZ(I+A-3,J+B-3,K+C-3)*GPW(N,25*(A-1)+5*(B-1)+C)
 370				CONTINUE
 360			CONTINUE
 350 		CONTINUE
	   	UZ(I,J,K)=UGP+VOB*GPW(N,126)
 340	CONTINUE

	RETURN
	END
	      
