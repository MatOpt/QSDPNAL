%% function stein transformation
%% Q(X) = X - AXA where A = D + L'*L;
%%************************************************

function QX = Qfun_stein(L,D,X,A)
if nargin > 3
   AX = A*X;
   QX = AX*A;
else
   AX = D*X + L'*(L*X);
   AXA = AX*D + (AX*(L'))*L;
   QX = AXA;
   QX = 0.5*(QX + QX');
end
