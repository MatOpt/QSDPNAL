   
function axb=QXfun_corr(L1,D1,L2,D2,x)
   
   ax=D1*x+L1'*(L1*x);
   axb=ax*D2+(ax*(L2'))*L2;
   axb=(axb+axb')/2;
   
end