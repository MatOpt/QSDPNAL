%%*****************************************************
%% NCM: nearest correlation matrix
%% type1 : diag(X) = e
%% type2 : Xe = e
%% type3 : diag(X) = e, Xe = e
%%
%% W = a nonnegative weight matrix. 
%%*****************************************************

   function [blk,At,C,b,L,U] = NCM_Hnorm(G,H,type)

   n = length(G); 
   HG = H.*G;
   normHG = norm(HG,'fro');
   
   C = randn(n);
   C = C*normHG^2/norm(C,'fro');
   C = (C+C')/2;

   blk{1,1} = 's'; blk{1,2} = n;  

   Acelltmp = cell(1,2*n); JJ = [1:n]';
   for k = 1:n; 
           Acelltmp{k} = spconvert([k,k,1;n,n,0]);
           II = k*ones(n,1);
           tmp = spconvert([II, JJ, ones(n,1); n, n, 0]);
           tmp = 0.5*(tmp + tmp');
           tmp(k,k) = 1;
           Acelltmp{k+n} = tmp;
   end
   if type == 1
       Acell = Acelltmp(1,1:n);
       b = ones(n,1);
   elseif type == 2
       Acell = Acelltmp(1,n+1:2*n);
       b = ones(n,1);
   elseif type == 3
       Acell = Acelltmp;
       b = ones(2*n,1);
   end
   
   Atmp = svec(blk(1,:),Acell,1); 
   At = Atmp{1}; 
      
      
      
    
   L = 0.7*min(min(G));
   U = 1.2*max(max(G));

%%*****************************************************
