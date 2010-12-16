A = sprandn(10,10,.5)
v = sprandn(10,1,.2)
v = full(v)+eps

A(:,3) = v;
A(:,3)