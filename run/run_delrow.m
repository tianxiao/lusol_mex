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

%% choose row to delete
i = 8;

%% delete the row
inform = lu.delrow(i)

%% test with a multiply
A1 = [A(1:i-1,:); A(i+1:end,:)];
x = ones(n,1);
b1 = A1*x;
b2 = lu.mulA(x);
norm(b1-b2(1:end-1))