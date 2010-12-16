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

%% generate row to add
w = ones(n,1);

%% add the row
[inform diag] = lu.addrow(w)

%% test
A1 = [A; w'];
x = rand(n,1);
b = A1*x;
blu = lu.mulA(x);

norm(b-blu)
