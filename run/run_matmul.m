% test matrix multiply methods

%% select rng seed
rand('twister',6548);

%% generate matrix
m = 20;
n = 30;
d = .2;
A = sprand(m,n,d);
%spy(A)

%% get lusol object
mylu = lusol(A);

%% get factors
L = mylu.L0();
U = mylu.U();

%% mulL
Xc1 = 5;
X1 = rand(m,Xc1);
Y1 = mylu.mulL(X1);
norm(L*X1 - Y1,1)

%% mulLt
Xc2 = 5;
X2 = rand(m,Xc2);
Y2  = mylu.mulLt(X2);
norm(L'*X2 - Y2,1)

%% mulU
Xc3 = 5;
X3 = rand(n,Xc3);
Y3   = mylu.mulU(X3);
norm(U*X3 - Y3,1)

%% mulUt
Xc4 = 5;
X4 = rand(m,Xc4);
Y4   = mylu.mulUt(X4);
norm(U'*X4 - Y4,1)

%% mulA
Xc5 = 5;
X5 = rand(n,Xc5);
Y5   = mylu.mulA(X5);
norm(A*X5 - Y5,1)

%% mulAt
Xc6 = 5;
X6 = rand(m,Xc6);
Y6   = mylu.mulAt(X6);
norm(A'*X6 - Y6,1)
