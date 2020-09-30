function QX = randQXfun(X,n)

rng('default');
r = min(n,5);
Atmp = randn(n,r);
A = Atmp*Atmp';
A = A/norm(A,'fro');
Btmp = randn(n,r);
B = Btmp*Btmp';
B = B/norm(B,'fro');
QX = (0.5/sqrt(n))*(A*X*B + B*X*A);

