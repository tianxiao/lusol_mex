clear
RandStream.setDefaultStream(RandStream('mt19937ar','seed',1));

% set up matrix
m = 20;
n = 20;
density = 0.4;
A = sprandn(m,n,density);

% run lusol
%options = lusol.luset();
lu = lusol();

[inform nsing depcol] = lu.factorize(A)