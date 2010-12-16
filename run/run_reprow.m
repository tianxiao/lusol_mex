clear
RandStream.setDefaultStream(RandStream('mt19937ar','seed',1));

%% set up matrix
m = 20;
n = 20;
density = 0.4;
A = sprandn(m,n,density);

%% run lusol
%options = lusol.luset();
lu = lusol(A);

%% choose a row for replacement
irep = 3;
r = ones(n,1);

%% replace the row
inform = lu.reprow(r,irep)

%% test with product
x = ones(n,1);
A1 = A;
A1(irep,:) = r;
b1 = A1*x;
b2 = lu.mulA(x);
norm(b1-b2)

%% test with solve
b3 = A1\x;
b4 = lu.solveA(x);
norm(b3 - b4)