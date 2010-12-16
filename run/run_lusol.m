clear
RandStream.setDefaultStream(RandStream('mt19937ar','seed',111111));

% set up matrix
m = 20;
n = 20;
density = 0.4;
A = sprandn(m,n,density);
%A(:,1) = A(:,2);

% initialize lusol
options = lusol.luset();
lu = lusol(1,options); % here the matrix 1 is factorized

% factorize A
[inform nsing depcol] = lu.factorize(A);
fprintf('\nafter factorize...\n')
fprintf('inform = %d\n',inform)
fprintf('nsing = %d\n',nsing)
fprintf('sum(depcol) = %d\n',sum(depcol))

% test multiply methods
xn = ones(n,1);
xm = ones(m,1);
fprintf('\ntesting multiply methods...\n')
fprintf('norm(A*xn - lu.mulA(xn),1) = %g\n',norm(A*xn - lu.mulA(xn),1))
fprintf('norm(A*xn - lu.mulL(lu.mulU(xn)),1) = %g\n',norm(A*xn - lu.mulL(lu.mulU(xn)),1))
fprintf('norm(A''*xm - lu.mulAt(xm),1) = %g\n',norm(A'*xm - lu.mulAt(xm),1))
fprintf('norm(A''*xm - lu.mulUt(lu.mulLt(xm)),1) = %g\n',norm(A''*xm - lu.mulUt(lu.mulLt(xm)),1))

% test solve
fprintf('\ntesting solve methods...\n')
fprintf('norm(A\\xn - lu.solveA(xn),1) = %g\n',norm(A*xn - lu.mulA(xn),1))
fprintf('norm(A\\xn - lu.solveU(lu.solveL(xn)),1) = %g\n',norm(A\xn - lu.solveU(lu.solveL(xn)),1))
fprintf('norm(A''\\xn - lu.solveAt(xn),1) = %g\n',norm(A'\xn - lu.solveAt(xn),1))
fprintf('norm(A''\\xn - lu.solveLt(lu.solveUt(xn)),1) = %g\n',norm(A'\xn - lu.solveLt(lu.solveUt(xn)),1))

% test repcol by negating column jrep
jrep = 4;
crep = -A(:,4);
B = A;
B(:,jrep) = -A(:,4);
fprintf('\ntesting repcol...\n')
[inform diag vnorm] = lu.repcol(crep,jrep);
fprintf('inform = %d\n',inform)
fprintf('diag = %d\n',diag)
fprintf('vnorm = %d\n',vnorm)
fprintf('norm(B*xn - lu.mulA(xn)) = %g\n',norm(B*xn - lu.mulA(xn)))
fprintf('norm(B\\xn - lu.solveA(xn)) = %g\n',norm(B\xn - lu.solveA(xn)))

% test repcol by replacing with dependent column
jrep = 8;
crep = A(:,1);
%B = A;
%B(:,jrep) = -A(:,4);
fprintf('\ntesting repcol, replacing with dependent column...\n')
[inform diag vnorm] = lu.repcol(crep,jrep);
fprintf('inform = %d\n',inform)
fprintf('diag = %d\n',diag)
fprintf('vnorm = %d\n',vnorm)
%fprintf('norm(B*xn - lu.mulA(xn)) = %g\n',norm(B*xn - lu.mulA(xn)))
%fprintf('norm(B\\xn - lu.solveA(xn)) = %g\n',norm(B\xn - lu.solveA(xn)))

