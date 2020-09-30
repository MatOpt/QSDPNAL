%% demo for QSDPNAL for solving QSDP-BIQ, QSDP-QAP, LS, LSI, NCM and NCMP problems

clear all
startup

BIQ = 1;
QAP = 1;
LS  = 1;
LSI = 1;
NCM = 1;
NCMP= 1;

%% Solve QSDP-BIQ
%% P): min 1/2 <X, QX> + <C,X>, s.t., AE X = bE, B X >= d, X >= 0, X \in K
%% where K = psd cone 
if BIQ
   G = biqread(['be250.2.sparse']);
   [blk,AEt,C,bE,AIt,bI] = biq_sqsdp(G,3);
   Q.QXfun = @(X) randQXfun(X,blk{1,2});
   options.nonnegative = 1;
   fprintf('\n #######################################');
   fprintf('\n ****************** %s ************','be250.2');
   fprintf('\n #######################################');   
   [obj,Z,W,QW,S,yE,yI,X,runhist,info] = ... 
          qsdpnal(blk,Q,AEt,bE,C,AIt,bI,options);
   %[obj,Z,W,QW,S,yE,yI,X,runhist,info] = qsdpnal_BIQ(blk,Q,AEt,bE,C,AIt,bI,options);
   clear blk Q AEt bE C AIt bI options
   clear obj Z W QW S yE yI X runhist info
end

%% Solve QSDP-QAP
%% P): min 1/2 <X, QX> + <C,X>, s.t., AE X = bE, B X >= d, X >= 0, X \in K
%% where K = psd cone 
if QAP
   [AA,BB] = qapread(['chr15c.dat']);
   [blk,AEt,C,bE] = qapAW(AA,BB,2);
   AEt = AEt{1}; C = C{1};
   AIt = []; bI = [];
   options.nonnegative = 1;
   Q.QXfun = @(X) randQXfun(X,blk{1,2});
   fprintf('\n ######################################');
   fprintf('\n ****************** %s ************','chr15c');
   fprintf('\n ######################################');
   [obj,Z,W,QW,S,yE,yI,X,runhist,info] = qsdpnal(blk,Q,AEt,bE,C,AIt,bI,options);
   %[obj,Z,W,QW,S,yE,yI,X,runhist,info] = qsdpnal_QAP(blk,Q,AEt,bE,C,AIt,bI,options);
   clear blk Q AEt bE C AIt bI options
   clear obj Z W QW S yE yI X runhist info
end

%% Solve least squares
%% (P) min {1/2||H(X)||_F^2 + <C,X> | A(X)=b , X \in K}
%% where K = psd cone 

if LS
   load('2d250n10nf.mat');
   Hinput = Bt;
   fprintf('\n ######################################');
   fprintf('\n ****************** %s ************','2d250n10nf');
   fprintf('\n ######################################');
   options.Qsolver = 'pcg';
   [obj,xi,S,yE,X,runhist,info] = ...
         LSnal(blk,Hinput,At,b,C,options);
   clear blk Hinput At b C options
   clear obj xi S yE X runhist info
end

%% Solve least squares with inequality constraints
%% (P) min {1/2||H(X)||_F^2 + <C,X> | A(X)=b, B(X)>=d, X \in K}
%% where K = psd cone 

if LSI
   load('3d250n10nf.mat');
   Hinput = Bt;
   fprintf('\n ######################################');
   fprintf('\n ****************** %s ************','3d250n10nf');
   fprintf('\n ######################################');
   options.Qsolver = 'pcg';
   [obj,xi,S,yE,yI,X,runhist,info] = ...
         LSnalI(blk,Hinput,At,b,AIt,bI,C,options);
   clear blk Hinput At b AIt bI C options
   clear obj xi S yE yI X runhist info
end

%% Solve the NCM problem
%% (P) min {1/2||H o X||_F^2 + <C,X> | A(X)=b, X \in K}
%% where K = psd cone 

if NCM
   options = [];
   load divergentGH.mat;
   tmp0 = ones(110,110);
   H0 = kron(tmp0,H);
   alpha = 0.1;
   load('Lymph.mat','S');
   n = length(S);
   H = H0([1:n],[1:n]);
   H = (H+H')/2;
   G = S([1:n],[1:n]);
   clear S;
   E = 2*rand(n)-1;
   E = triu(E) + triu(E,1)';
   E = 0.5*(E + E');
   G = (1-alpha)*G + alpha*E;
   G = 0.5*(G+G');
   G = G - diag(diag(G)) + eye(n);          
   Hscale = sqrt(norm(H.*H.*G,'fro'));
   H = H/(Hscale);
   H2 = H.*H;             
   type = 1;
   [blk,At,~,b,~,~] = NCM_Hnorm(G,H,type);
   C = zeros(blk{1,2});
   C = C - H2.*G;   
   X0 = G; 
   fprintf('\n #######################################');
   fprintf('\n ****************** %s ************','Lymph');
   fprintf('\n #######################################');   

   [obj,Xi,S,yE,X,runhist,info] = qsdpnal_NCM(blk,H,At,b,C,options,X0);
   clear blk H At b C options X0
   clear obj Xi S yE X runhist info
end

%% Solve the NCM problem
%% (P) min {1/2||H o X||_F^2 + <C,X> | A(X)=b, X \in K, X\in P}
%% where K = psd cone, P is ployhedral set e.g., P = {X | L<= X <= U}.
if NCMP
   load divergentGH.mat;
   tmp0 = ones(110,110);
   H0 = kron(tmp0,H);
   alpha = 0.05;
   load('ER.mat','S');
   n = length(S);
   H = H0([1:n],[1:n]);
   H = (H+H')/2;
   G = S([1:n],[1:n]);
   clear S;
   E = 2*rand(n)-1;
   E = triu(E) + triu(E,1)';
   E = 0.5*(E + E');
   G = (1-alpha)*G + alpha*E;
   G = 0.5*(G+G');
   G = G - diag(diag(G)) + eye(n);          
   Hscale = sqrt(norm(H.*H.*G,'fro'));
   H = H/(Hscale);
   H2 = H.*H;             
   type = 1;
   [blk,At,~,b,~,~] = NCM_Hnorm(G,H,type);
   C = zeros(blk{1,2});
   C = C - H2.*G;   
   options.lowerbound = -.5; % Lower bound indices, i.e., X >= -0.5
   SignP.L = 1*ones(n); 
   options.SignP = SignP; 
   X0 = G; 
   fprintf('\n #######################################');
   fprintf('\n ****************** %s ************','ER');
   fprintf('\n #######################################');   
   [obj,Z,Xi,S,yE,X,runhist,info] = qsdpnal_NCMP(blk,H,At,b,C,options,X0);
   clear blk H At b C options X0
   clear obj Z Xi S yE X runhist info
end


