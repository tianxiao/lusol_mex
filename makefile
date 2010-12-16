# general settings
MEXFLAGS = -g -v
FFLAGS = -fPIC -g

# Ubuntu settings with gfortran
FC = gfortran-4.3
FCLIB = -L/usr/lib/gcc/x86_64-linux-gnu/4.3/ -lgfortran
MEX = mex
OUTMEX = lusol_mex.mexa64

# Mac OS X 10.6 settings
#FC = gfortran
#FCLIB = -lgfortran
#MEX = /Applications/MATLAB_R2009aSV.app/bin/mex
#OUTMEX = lusol_mex.mexmaci

OBJECTS = lusol.o

$(OUTMEX) : lusol_mex.c $(OBJECTS)
	$(MEX) $(MEXFLAGS) $(FCLIB) lusol_mex.c $(OBJECTS)

# cleanup
.PHONY: clean clean_mex clean_all
clean : 
	rm -f *.o

clean_mex :
	rm -f *.mex*

clean_all : clean clean_mex
