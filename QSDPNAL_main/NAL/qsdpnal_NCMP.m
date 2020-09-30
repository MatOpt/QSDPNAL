%%*************************************************************************
%% QSDPNAL: Phase II
%% An inexact semismooth Newton based ALM for solving NCM problems with 
%% Polyhderal constraints:
%% (P) min {1/2||H o X||_F^2 + <C,X> | A(X)=b, X \in K, X\in P}
%% where K = psd cone, P is ployhedral set e.g., P = {X | L<= X <= U}.
%% (P') min {1/2 ||Y||_F^2 + <C,X> | A(X) = b, H o X = Y, X\in K, X\in P}
%%
%% (D') max -\delta_{P}^*(-Z) - 1/2 ||Xi||^2  + <b,y>
%%      s.t. Z + H o Xi + S + A*y = C,
%%       S \in PSD.
%%*************************************************************************
%% QSDPNAL: 
%% Copyright (c) 2015 by
%% Xudong Li, Defeng Sun, and Kim-Chuan Toh
%%*************************************************************************

   function [obj,Z,Xi,S,yE,X,runhist,info] = qsdpnal_NCMP(blk,H,A,b,C,options,X,Z,Xi,S,yE)
%% General options
   rng('default');
   maxiter = 5000;
   stoptol = 1e-6;
   printyes = 1;
   printminoryes=1;
   startfun = @SCB_NCMP; % Algorithm QSDPNAL-Phase I for initial point
   sigma0 = 1;
   sGS.iter = 0;
   sGS.time = 0;
   use_proximal = 1;
   restart = 1;
%% options from users
   if isfield(options,'maxiter');  maxiter  = options.maxiter; end
   if isfield(options,'stoptol');  stoptol  = options.stoptol; end
   if isfield(options,'printyes'); printyes = options.printyes; end
   if isfield(options,'printminoryes'); printminoryes = options.printminoryes; end
   if isfield(options,'sigma');    sigma0 = options.sigma; end
%% information for \cP
   nonnegative = 0; L = []; U = []; SignP = []; fV = [];
   if isfield(options,'upperbound'); U = options.upperbound; end
   if isfield(options,'lowerbound'); L = options.lowerbound; end
   if isfield(options,'SignP');      SignP = options.SignP; end
   if isfield(options,'fV');         fV = options.fV; end
   if isfield(options,'nonnegative');  nonnegative = options.nonnegative; end

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
   tstart = clock;
   n = blk{1,2};
   Zero = sparse(n,n);
   Fnorm = @(X) mexFnorm(X);
%% AEmap, AETmap, AImap, AITmap, Qmap
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
%%
%% initial scaling
%% 
   Atorg = At;  borg = b; Lorg = L; Uorg = U; fVorg = fV;
   Corg = C; Horg = H;
   normHorg = Fnorm(Horg);
   normCorg = 1 + Fnorm(Corg);
   normborg = 1 + norm(borg);
%%
   if ~exist('X','var'); X = Zero; end
   if ~exist('yE','var'); yE = zeros(length(b),1); end
   if ~exist('S','var'); S = Zero; end
   if ~existA; yE = []; b = []; end
   if ~exist('Xi','var'); Xi = Zero; end
   if ~exist('Z','var'); Z = Zero; end
%%
%% generate initial point through QSDPNAL-Phase I (SCB_qsdp)
   sGS.op.stoptol = 6e-5; % tol of initial stage
   sGS.op.maxiter = 250; % maxiter of initial stage
   sGS.op.lowerbound = Lorg; sGS.op.SignP = SignP;
   sGS.op.upperbound = Uorg; sGS.op.fV = fVorg;
   sGS.op.nonnegative = nonnegative;
   sGS.op.phase2 = 1;
   fprintf('\n********************************************************');
   fprintf('\n                  QSDPNAL Phase I                       ');
   fprintf('\n********************************************************');
   [obj,Z,Xi,S,yE,X,runhist_sGS,info_sGS] = startfun(blk,Horg,A,borg,Corg,sGS.op,X,Z,Xi,S,yE);
   fprintf('\n -------------------------------------------------------');
   sGS.Xi0 = Xi; sGS.S0 = S; sGS.Z0 = Z;
   sGS.yE0= yE;  sGS.X0=X;
   nE = info_sGS.nE; nXi = info_sGS.nXi;
   bscale = info_sGS.bscale; Cscale = info_sGS.Cscale;
   sigma = info_sGS.sigma;
   normA = info_sGS.normA; DA = info_sGS.DA;
   Amap = info_sGS.Amap; ATmap = info_sGS.ATmap;
   b = info_sGS.b; 
   C = info_sGS.C; 
   sGS.iter = sGS.iter + info_sGS.iter;
   sGS.time = sGS.time + info_sGS.totaltime;  
%% 
   X = X/bscale; yE = (normA.*yE)/Cscale; 
   S = S/Cscale; HXi = Horg.*Xi/Cscale;
   Xi = Xi/sqrt(bscale*Cscale); Z=Z/Cscale;
   L = L/bscale; U = U/bscale; fV = fV/bscale;
   H = sqrt(bscale/Cscale)*Horg; 
   objscale = bscale*Cscale;
   At = At*DA;
%%-------------------------------------------------------------------
%% end of intial stage
%%-------------------------------------------------------------------
%% start QSDPNAL-Phase II
   parmain.existA = existA;
   parmain.bscale = bscale; parmain.Cscale = Cscale; parmain.objscale = objscale;
   parmain.nE = nE;   parmain.n = n;  
   parmain.normA = normA;  parmain.DA = DA;
   parmain.HXi  = HXi;
   parmain.normHorg = normHorg;
   parmain.normborg = normborg;
   parmain.normCorg = normCorg;
   parmain.use_proximal = use_proximal;

   parmain.nonnegative = nonnegative; parmain.L = L; parmain.U = U;
   parmain.SignP = SignP; parmain.fV = fV;
   
   parmain.sigma = sigma;
   parmain.tstart= tstart;
   parmain.maxiter = maxiter;
   parmain.stoptol = stoptol;
   parmain.printyes = printyes;
   parmain.sGS = sGS;
   parmain.restart = restart;

   parmain.Lorg = Lorg; parmain.Uorg = Uorg;  parmain.fVorg = fVorg;
   parmain.borg = borg; parmain.Corg = Corg; 
   parmain.Horg = Horg; parmain.A  = A;
   parmain.startfun = startfun;
   
   parmain.info_sGS = info_sGS;
   parmain.obj = obj;
   
   Ainput.Amap =  Amap; Ainput.ATmap =  ATmap; 
   % main solver QSDPNAL-Phase II with restart technique   
   [X,Z,yE,S,Xi,info_main,runhist] = ...
    nal_NCMPmain(blk,H,At,Ainput,b,C,parmain,X,Z,yE,S,Xi);
   iter = info_main.iter;
   msg  = info_main.msg;
   bscale = info_main.bscale;
   Cscale = info_main.Cscale;
   Ztmp  = info_main.Ztmp;
   sGS = info_main.sGS;
%%-----------------------------------------------------------------
%% recover orignal variables
%%-----------------------------------------------------------------
   if (iter == maxiter)
      msg = ' maximum iteration reached';
      info.termcode = 3;
   end
   X = X*bscale; normX = norm(X,'fro');
   HX = Horg.*X;
   yE = Cscale*(DA*yE); 
   Z = Cscale*Z; normZ = norm(Z,'fro'); Ztmp = Ztmp*bscale;
   S = Cscale*S; normS = norm(S,'fro');
   Xi = Xi*sqrt(bscale*Cscale);
   HXi = Horg.*Xi;
   AX = Amap0(X); AtyE = ATmap0(yE);
   if existA; bTyE = borg'*yE; end
   primobj = 0.5*Fnorm(HX)^2 + sum(sum(Corg.*X));
   dualobj = bTyE + sum(sum(Z.*Ztmp)) - 0.5*Fnorm(Xi)^2;
   obj = [primobj,dualobj];
   etaQ = Fnorm(Xi + Horg.*X)^2/(1+normHorg^2);
   if (true)
      etaSp = compnormXn(X) /(1+normX);
      etaSc = full(abs(sum(sum(X.*S)))) /(1+normX +normS);
      etaS = max(etaSp,etaSc);
      projPX = proj_LU_new(X,nonnegative,SignP,Lorg,Uorg,fVorg);
      etaZp = Fnorm(X-projPX)/(1 + normX);
      ZtmX = Ztmp - X;
      normZtmX = norm(ZtmX,'fro');
      trZtmX = full(abs(sum(sum(ZtmX.*Z))));       
      etaZc = trZtmX/(1 + normZ + normZtmX);            
      etaZ = max(etaZp,etaZc);
   end
   etab = norm(AX-borg)/normborg;
   primfeasorg = etab;
   dualfeasorg = norm(Z + S+AtyE+HXi-Corg,'fro')/normCorg;
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
   info.etaZ = etaZ;
   info.etaS = etaS;
   info.etab = etab;
   info.primfeasorg = primfeasorg; info.dualfeasorg = dualfeasorg;
   info.etaQ = etaQ;
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
      fprintf('\n  etaQ = %3.2e', etaQ); fprintf('  etaS = %3.2e', etaS);
      fprintf('\n  etab = %3.2e', etab); fprintf('  etaZ = %3.2e', etaZ);
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
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
