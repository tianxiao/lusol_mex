# general settings
#MEXFLAGS = -g -v
MEXFLAGS = -g

# Ubuntu settings with gfortran
# FC = gfortran-4.3
# FFLAGS = -fPIC -g
# FCLIB = -L/usr/lib/gcc/x86_64-linux-gnu/4.3/ -lgfortran
# MEX = mex
# OUTMEX = lusol_mex.mexa64

# Mac OS X 10.6 settings
FC = gfortran
FFLAGS = -arch x86_64 -fPIC -g
FCLIB = -L/usr/local/Cellar/gfortran/4.2.4-5664/lib/gcc/i686-apple-darwin10/4.2.1/x86_64/ -lgfortran
MEX = /Applications/MATLAB_R2011a_Student.app/bin/mex
OUTMEX = lusol_mex.mexmaci64

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
