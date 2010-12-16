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

%% choose column to delete
j = 10;

%% delete the column
inform = lu.delcol(j);

%% test with product
A1 = [A(:,1:j-1) A(:,j+1:end)];
x = ones(n-1,1);
b1 = A1*x;
b2 = lu.mulA(x);
norm(b1-b2)