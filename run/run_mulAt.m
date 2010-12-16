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

%% generate x vector
x = randn(n,1);

%% test
y1 = A'*x;
y2 = lu.mulAt(x);

norm(y1-y2)