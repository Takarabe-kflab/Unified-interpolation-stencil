#   
#	File name			: Makefile
#	Host machine		: 
#	Execution file name	: 
#	Date			: 
#	Version			: 
#	Author			: Koji Fukagata
#

DEST			= .

EXTHDRS			=

HDRS			= 	par.f	\
				common.f

LIBS			=

INPS			= inp.csv

COMPILER		= ifort

OPTFLAG			= -xHost -mcmodel=large

MAKEFILE		= Makefile

PRINT			= a2ps --prologue=color --pretty-print=fortran -o list.ps

PROGRAM			= Sphere3D

SRCS			=	main.f                  \
				init/grid1.f 	        \
				init/init-io.f	        \
            	solver/outlet.f         \
				solver/ubc-io.f      	 \
				solver/mkf.f        	\
				solver/mku.f	        \
				solver/mkp.f	        \
				matrix/afy.f	        \
				matrix/tridiag.f	\
				stat/stat0.f		\
				stat/stat1.f		\
				misc/chkcont.f	        \
				misc/chkcfl.f	        \
				misc/rawout.f	        \
				misc/outmat.f	        \
				fft/fft-fornberg.f	\
				fft/p2s-nx2-fornberg-1.1.f	\
				ibm/ker.f	\
				ibm/ghp.f

OBJS			= $(SRCS:.f=.o) 	

LIBS			=  	

.f.o:
			$(COMPILER) $(OPTFLAG) -c $*.f -o $*.o -lm

all:			$(PROGRAM)

$(PROGRAM):		$(OBJS) $(LIBS)
			@echo -n "Loading $(PROGRAM) ... "
			@$(COMPILER) $(OPTFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)
			@echo "done"

clean:;			@rm -f $(SRCS:.f=.o)

depend:;		@mkmf -f $(MAKEFILE) PROGRAM=$(PROGRAM) DEST=$(DEST)

index:;			@ctags -wx $(HDRS) $(SRCS)

install:		$(PROGRAM)
			@echo Installing $(PROGRAM) in $(DEST)
			@install -s $(PROGRAM) $(DEST)

print:;			@cat Makefile $(HDRS) $(SRCS) |$(PRINT)  

log:;			@cat *.log > $(PROGRAM).LOG

pkg:;			@tar cvf $(PROGRAM).tar $(HDRS) $(SRCS) $(MAKEFILE)

data:;                  @tar cvf data.tar inp.csv out Out/stat Out/trace 
 
scalar:;                @f77 -loglist -pmfunc -listlang=euc -W0,'MP(P(0)),LANGLVL(SAVE(0))' -c matrix/tridiag.f -o matrix/tridiag.o -lm;

program:		$(PROGRAM)

tags:			$(HDRS) $(SRCS); @ctags $(HDRS) $(SRCS)

update:			$(DEST)/$(PROGRAM)

#$(DEST)/$(PROGRAM):	$(SRCS) $(LIBS) $(HDRS) $(EXTHDRS)
#			@make -f $(MAKEFILE) DEST=$(DEST) install
#
#       Makefile  : End of File
#














