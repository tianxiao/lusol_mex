% test matrix solve methods

%% select rng seed
rand('twister',6548);

%% generate matrix
m = 20;
n = 30;
d = .2;
A = sprand(m,n,d);
spy(A)

%% get lusol object
mylu = lusol(A);

%% get factors
L = mylu.L0();
U = mylu.U();

%% solveL
Bc1 = 5;
B1 = rand(m,Bc1);
[X1 inform1] = mylu.solveL(B1);
norm(L*X1 - B1,1)

%% solveLt
Bc2 = 5;
B2 = rand(m,Bc2);
[X2 inform2] = mylu.solveLt(B2);
norm(L'*X2 - B2,1)

%% solveU
Bc3 = 5;
B3 = rand(m,Bc3);
[X3 inform3 resid3] = mylu.solveU(B3);
norm(U*X3 - B3,1)

%% solveUt
Bc4 = 5;
B4 = rand(n,Bc4);
[X4 inform4 resid4] = mylu.solveUt(B4);
norm(U'*X4 - B4,1)

%% solveA
Bc5 = 5;
B5 = rand(m,Bc5);
[X5 inform5 resid5] = mylu.solveA(B5);
norm(A*X5 - B5,1)

%% solveAt
Bc6 = 5;
B6 = rand(n,Bc6);
[X6 inform6 resid6] = mylu.solveAt(B6);
norm(A'*X6 - B6,1)
