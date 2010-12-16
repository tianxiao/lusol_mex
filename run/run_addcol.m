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

%% generate colum to add
v = ones(m,1);

%% add the column
[inform diag vnorm] = lu.addcol(v)

%% test it
A1 = [A v];
x = ones(n+1,1);

b1 = A1*x;
b2 = lu.mulA(x);

norm(b1-b2)
