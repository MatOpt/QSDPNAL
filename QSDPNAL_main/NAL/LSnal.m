%%*************************************************************************
%% QSDPNAL: Phase II 
%% An inexact semismooth Netwon based ALM for solving the following Least Squre problem
%% (P) min {1/2||H(X)||_F^2 + <C,X> | A(X)=b, X \in K}
%% where K = psd cone 
%% (P') min {1/2 ||Y||_F^2 + <C,X> | A(X) = b, H(X) = Y, X\in K}
%%
%% (D') max -1/2 ||Xi||^2  + <b,y>
%%      s.t. H*(Xi) + S + A*y = C,
%%       S \in PSD.
%%*************************************************************************
%% QSDPNAL: 
%% Copyright (c) 2015 by
%% Xudong Li, Defeng Sun, and Kim-Chuan Toh
%%*************************************************************************
   function [obj,xi,S,yE,X,runhist,info] = LSnal(blk,H,A,b,C,options,X,xi,S,yE)
%% General options
   rng('default');
   maxiter = 5000;
   stoptol = 1e-6;
   printyes = 1;
   printminoryes=1;
   startfun = @SCB_LS; % Algorithm QSDPNAL-Phase I for initial point
   sigma0 = 1;
   sGS.iter = 0;
   sGS.time = 0;
   use_proximal = 1;
   % options from users
   if isfield(options,'maxiter');  maxiter  = options.maxiter; end
   if isfield(options,'stoptol');  stoptol  = options.stoptol; end
   if isfield(options,'printyes'); printyes = options.printyes; end
   if isfield(options,'printminoryes'); printminoryes = options.printminoryes; end
   if isfield(options,'sigma');    sigma0 = options.sigma; end
   
   
   fprintf('\n');
   fprintf('\n============================================================================================================');
   fprintf('\n************************************************************************************************************');
   fprintf('\n*  QSDPNAL: A two-phase augmented Lagrangian method for convex quadratic semidefinite programming          *');
   fprintf('\n*  Authors: Xudong Li, Defeng Sun, and Kim-Chuan Toh                                                       *');
   fprintf('\n************************************************************************************************************');
   fprintf('\n============================================================================================================');
   fprintf('\n************************************************************************************************************');
   fprintf('\n*                                             Reference                                                    *');
   fprintf('\n*  Xudong Li, Defeng Sun, and Kim-Chuan Toh, QSDPNAL: A two-phase augmented Lagrangian method for convex   *');
   fprintf('\n*  quadratic semidefinite programming, Mathematical Programming Computation, 10 (2018), pp. 703--743.      *');
   fprintf('\n************************************************************************************************************');
   fprintf('\n============================================================================================================');   
   fprintf('\n');
   
%% problem property
   existA = ~isempty(A);
   existH = ~isempty(H);
   tstart = clock;
   n = blk{1,2};
   nbar = n*(n+1)/2;
   Zero = sparse(n,n);
   Fnorm = @(X) mexFnorm(X);
%% AEmap, AETmap, Qmap
   if existA
       if isstruct(A)
          if isfield(A,'At'); At = A.At; end
          if isfield(A,'Amap'); Amap0 = A.Amap; end
          if isfield(A,'ATmap'); ATmap0 = A.ATmap; end    
       else
          At = A;
          Amap0  = @(X) sdpnalAXfun(blk,At,X);
          ATmap0 = @(y) sdpnalAtyfun(blk,At,y);    
       end
   else
      Amap0 = @(X) []; ATmap0 = @(y) Zero;
   end
   if existH
      if isstruct(H)
         if isfield(H,'Ht'); Ht = H.Ht; end
         if isfield(H,'Hmap'); Hmap0 = H.Hmap; end
         if isfield(H,'HTmap'); HTmap0 = H.HTmap; end
      else
         Ht = H;
         Hmap0 = @(X) sdpnalAXfun(blk,Ht,X);
         HTmap0 = @(y) sdpnalAtyfun(blk,Ht,y);
      end
   else
      Hmap0 = @(X) []; HTmap0 = @(y) Zero;
   end
%%
%%
   if ~exist('X','var'); X = Zero; end
   if ~exist('yE','var'); yE = zeros(length(b),1); end
   if ~exist('S','var'); S = Zero; end
   if ~existA; yE = []; b = []; end
   if ~exist('Xi','var'); xi = zeros(size(Ht,2),1); end
%% initial scaling
%% 
   Atorg = At;  borg = b;
   Corg = C; Htorg = Ht;
   normHorg = Fnorm(Htorg);
   normCorg = 1 + Fnorm(Corg);
   normborg = 1 + norm(borg);
%%
%% generate initial point through QSDPNAL-Phase I (SCB_qsdp)
   sGS.op.stoptol = 1e-3;
   sGS.op.maxiter = 1000;
   sGS.op.phase2 = 1;
   fprintf('\n********************************************************');
   fprintf('\n                  QSDPNAL Phase I                       ');
   fprintf('\n********************************************************');
   [obj,xi,S,yE,X,runhist_sGS,info_sGS] = startfun(blk,Htorg,A,borg,Corg,sGS.op,X,xi,S,yE);
   fprintf('\n -------------------------------------------------------');
   sGS.xi0 = xi; sGS.S0 = S; 
   sGS.yE0= yE;  sGS.X0=X;
   nE = info_sGS.nE; 
   nxi = length(xi);
   bscale = info_sGS.bscale; Cscale = info_sGS.Cscale;
   sigma = info_sGS.sigma;
   normA = info_sGS.normA; DA = info_sGS.DA;
   Amap = info_sGS.Amap; ATmap = info_sGS.ATmap;
   b = info_sGS.b; 
   C = info_sGS.C; 
   sGS.iter = sGS.iter + info_sGS.iter;
   sGS.time = sGS.time + info_sGS.totaltime;  
%% 
   At = At*DA;
   X = X/bscale; yE = (normA.*yE)/Cscale; 
   S = S/Cscale; Htxi = HTmap0(xi)/Cscale;
   xi = xi/sqrt(bscale*Cscale);
   objscale = bscale*Cscale;
%%-------------------------------------------------------------------
%% end of intial stage
%%-------------------------------------------------------------------
%% start QSDPNAL-Phase II
   parmain.existH = existH;  parmain.existA = existA;
   parmain.bscale = bscale; parmain.Cscale = Cscale; parmain.objscale = objscale;
   parmain.nE = nE;   parmain.n = n; parmain.nxi = nxi;
   parmain.normA = normA;  
   parmain.Htorg = Htorg;
   parmain.Htxi  = Htxi;
   parmain.normHorg = normHorg;
   parmain.normborg = normborg;
   parmain.normCorg = normCorg;
   parmain.use_proximal = use_proximal;
   
   parmain.sigma = sigma;
   parmain.tstart= tstart;
   parmain.maxiter = maxiter;
   parmain.stoptol = stoptol;
   parmain.printyes = printyes;
   parmain.sGS = sGS;
   
   parmain.info_sGS = info_sGS;
   parmain.obj = obj;
   
   Ainput.Amap =  Amap; Ainput.ATmap =  ATmap; Ainput.At =  At;
   % main solver QSDPNAL-Phase II with restart technique
   [X,yE,S,xi,info_main,runhist] = ...
      nal_LSmain(blk,Hmap0,HTmap0,Ainput,b,C,parmain,X,yE,S,xi);
   iter = info_main.iter;
   msg  = info_main.msg;
   bscale = info_main.bscale;
   Cscale = info_main.Cscale;
%%-----------------------------------------------------------------
%% recover orignal variables
%%-----------------------------------------------------------------
   if (iter == maxiter)
      msg = ' maximum iteration reached';
      info.termcode = 3;
   end
   X = X*bscale; normX = norm(X,'fro');
   HX = Hmap0(X);
   yE = Cscale*(DA*yE); 
   S = Cscale*S; normS = norm(S,'fro');
   xi = xi*sqrt(bscale*Cscale);
   Htxi = HTmap0(xi);
   AX = Amap0(X); AtyE = ATmap0(yE);
   if existA; bTyE = borg'*yE; end
   primobj = 0.5*norm(HX)^2 + sum(sum(Corg.*X));
   dualobj = bTyE - 0.5*Fnorm(xi)^2;
   obj = [primobj,dualobj];
   etaQ = Fnorm(xi + Hmap0(X))^2/(1+normHorg^2);
   if (true)
      etaSp = compnormXn(X) /(1+normX);
      etaSc = full(abs(sum(sum(X.*S))) /(1+normX +normS));
      etaS = max(etaSp,etaSc);
   end
   primfeasorg = norm(AX-borg)/normborg;
   dualfeasorg = norm(S+AtyE+Htxi-Corg,'fro')/normCorg;
   ttime = etime(clock,tstart);   
   runhist.nE = nE;
   runhist.ns = n;
   runhist.iter = iter;
   runhist.totaltime = ttime;
   runhist.primobjorg = primobj; 
   runhist.dualobjorg = dualobj;    
   runhist.relgaporg = (primobj-dualobj)/(1+abs(primobj)+abs(dualobj));  
   runhist.primfeasorg(iter) = primfeasorg;    
   runhist.dualfeasorg(iter) = dualfeasorg; 
   runhist.relerr = max([dualfeasorg,primfeasorg,etaS]);
   info.nE = nE;
   info.ns = n;
   info.iter = iter;
   info.totaltime = ttime;
   info.primobjorg = primobj;
   info.dualobjorg = dualobj;
   info.relgaporg = (primobj-dualobj)/(1+abs(primobj)+abs(dualobj)); 
   info.relerrS = etaS;
   info.primfeasorg = primfeasorg; info.dualfeasorg = dualfeasorg;
   info.relerrQ = etaQ;
   info.relerr = max([dualfeasorg,primfeasorg,etaS,etaQ]);
   info.minX = min(min(X));
   info.maxX = max(max(X));
   info.sGS_iter = sGS.iter;
   info.sGS_time = sGS.time;
   info.totalNT = sum(runhist.innerNT);
   info.totalCG = sum(runhist.innerCG);
   if (printminoryes) 
      if ~isempty(msg); fprintf('\n %s',msg); end
      fprintf('\n--------------------------------------------------------------');
      fprintf('------------------');
      fprintf('\n  number iter = %2.0d',iter);      
      fprintf('\n  time = %3.2f',ttime);       
      %fprintf('\n  time per iter = %5.4f',ttime/iter); 
      fprintf('\n  time of Phase I = %3.2f', sGS.time);
      fprintf('\n  iter of Phase I = %2.0d', sGS.iter);
      fprintf('\n  primobj = %9.8e',runhist.primobjorg);       
      fprintf('\n  dualobj = %9.8e',runhist.dualobjorg);            
      fprintf('\n  primfeasorg = %3.2e, dualfeasorg = %3.2e, relgaporg = %3.2e',...
	      runhist.primfeasorg(end),runhist.dualfeasorg(end),runhist.relgaporg(end));
      fprintf('\n  relerrQ = %3.2e', etaQ);
      fprintf('\n  min(X)    = %3.2e, max(X)    = %3.2e',...
          info.minX,info.maxX); 
      fprintf('\n  total innerNT steps = %3d, total innerCG = %3d', info.totalNT, info.totalCG);
      fprintf('\n--------------------------------------------------------------');
      fprintf('------------------\n');
   end
%%**********************************************************************
    function normXn = compnormXn(X)
    
    d = eig(full(X));   
    normXn = norm(d(d<0));     
%%**********************************************************************
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
