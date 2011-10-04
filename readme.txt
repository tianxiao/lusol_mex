lusol_mex provides a Matlab interface to Michael Saunders' LUSOL.

LUSOL maintains LU factors of a square or rectangular sparse matrix A.

The software for LUSOL is provided by SOL, Stanford University under the terms
of the Common Public License (CPL):
http://oss.software.ibm.com/developerworks/opensource/license-cpl.html

The website for LUSOL is:
http://www.stanford.edu/group/SOL/software/lusol.html

2010-11-12: First version of the readme file.
2010-12-15: Added xunit tests.
2011-01-17: added tutorial in doc/

Please send comments regarding the mex interface to:
  Nick Henderson
  nwh@stanford.edu
  Institute for Computational and Mathematical Engineering
  Stanford University

--------------------------------------------------------------------------------

Directories

doc      some documentation and notes
fspec    files and code to automatically generate lusol_mex.c
html     ftagshtml documentation of lusol.f
run      example matlab scripts to drive lusol

Files

lusol.f       fortran 77 LUSOL code 
lusol.m       matlab class to drive LUSOL
lusol_mex.c   automatically generated mex interface to lusol
makefile      makefile to build the interface
test_lusol.m  matlab xunit test cases for lusol.m, not complete coverage
