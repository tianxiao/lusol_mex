# general settings
MEX = mex
OUTMEX = lusol_mex.mexa64
#MEXFLAGS = -g -v
MEXFLAGS = 

# Ubuntu settings with gfortran
#FC = gfortran
#FFLAGS = -fPIC
#FCLIB = -L/usr/lib/gcc/x86_64-unknown-linux-gnu/4.3.4/ -lgfortran
#FCLIB = -lgfortran

# Mac OS X 10.6 settings
#FC = gfortran
#FFLAGS = -arch x86_64 -fPIC -g
#FCLIB = -L/usr/local/Cellar/gfortran/4.2.4-5664/lib/gcc/i686-apple-darwin10/4.2.1/x86_64/ -lgfortran
#MEX = /Applications/MATLAB_R2011a_Student.app/bin/mex
#OUTMEX = lusol_mex.mexmaci64

$(OUTMEX) : lusol_mex.c lusol.f
	$(MEX) $(MEXFLAGS) $(FCLIB) lusol_mex.c lusol.f -o $(OUTMEX)

# cleanup
.PHONY: clean_mex
clean_mex :
	rm -f *.mex*
