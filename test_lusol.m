% test_lusol simple LUSOL test cases
% 
% This test suite uses matlab xunit:
% http://www.mathworks.com/matlabcentral/fileexchange/22846-matlab-xunit-test-framework
%
% It does not provide complete coverage of all lusol.m methods.
%
% 2010-12-15 (nwh) first version, mulLt and mulAt fail at this point
%

function test_suite = test_lusol
initTestSuite;
end

function s = setup
RandStream.setDefaultStream(RandStream('mt19937ar','seed',9999));
s.n = 50;
s.d = .2;
s.A = sprand(s.n,s.n,s.d);
s.y = randn(s.n,1);
s.tol = 1e-8;
end

function test_solveA(s)
mylu = lusol(s.A);
x = mylu.solveA(s.y);
assertVectorsAlmostEqual(s.A*x,s.y);
end

function test_solveAt(s)
mylu = lusol(s.A);
x = mylu.solveAt(s.y);
assertVectorsAlmostEqual(s.A'*x,s.y);
end

function test_mulA(s)
mylu = lusol(s.A);
x = mylu.mulA(s.y);
assertVectorsAlmostEqual(s.A*s.y,x);
end

function test_LU(s)
mylu = lusol(s.A);
L = mylu.L0();
U = mylu.U();
err = norm(s.A - L*U,'inf');
assertElementsAlmostEqual(err,0)
end

function test_LU_perm(s)
mylu = lusol(s.A);
[L p1] = mylu.L0();
[U p2 q] = mylu.U();
assertEqual(p1,p2);
err = norm(s.A(p1,q) - L*U,'inf');
assertElementsAlmostEqual(err,0)
end

function test_solveL(s)
mylu = lusol(s.A);
L = mylu.L0();
x1 = L\s.y;
x2 = mylu.solveL(s.y);
assertVectorsAlmostEqual(x1,x2);
end

function test_solveLt(s)
mylu = lusol(s.A);
L = mylu.L0();
x1 = L'\s.y;
x2 = mylu.solveLt(s.y);
assertVectorsAlmostEqual(x1,x2);
end

function test_solveU(s)
mylu = lusol(s.A);
U = mylu.U();
x1 = U\s.y;
x2 = mylu.solveU(s.y);
assertVectorsAlmostEqual(x1,x2);
end

function test_solveUt(s)
mylu = lusol(s.A);
U = mylu.U();
x1 = U'\s.y;
x2 = mylu.solveUt(s.y);
assertVectorsAlmostEqual(x1,x2);
end

function test_mulL(s)
mylu = lusol(s.A);
L = mylu.L0();
x1 = L*s.y;
x2 = mylu.mulL(s.y);
assertVectorsAlmostEqual(x1,x2);
end

function test_mulU(s)
mylu = lusol(s.A);
U = mylu.U();
x1 = U*s.y;
x2 = mylu.mulU(s.y);
assertVectorsAlmostEqual(x1,x2);
end

function test_mulUt(s)
mylu = lusol(s.A);
U = mylu.U();
x1 = U'*s.y;
x2 = mylu.mulUt(s.y);
assertVectorsAlmostEqual(x1,x2);
end

function test_repcol(s)
j = 1; % column index to replace
v = ones(s.n,1); % column to insert

% replace in mylu
mylu = lusol(s.A);
inform = mylu.repcol(v,j);
assertEqual(inform,0);

% replace in s.A and check with a solve
s.A(:,j) = v;
x1 = s.A\s.y;
x2 = mylu.solveA(s.y);
assertElementsAlmostEqual(x1,x2)
end

% tests that fail

% function test_mulAt(s)
% mylu = lusol(s.A);
% x = mylu.mulAt(s.y);
% assertVectorsAlmostEqual(s.A'*s.y,x);
% end

% function test_mulLt(s)
% mylu = lusol(s.A);
% L = mylu.L0();
% x1 = L'*s.y;
% x2 = mylu.mulLt(s.y);
% assertVectorsAlmostEqual(x1,x2);
% end