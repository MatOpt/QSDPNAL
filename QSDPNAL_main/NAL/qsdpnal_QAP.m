%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% QSDPNAL: Phase II
%% An inexact semismooth Netwon based ALM for solving the following QSDP problem
%% QSDP-QAP problem
%% D'): max -(\delta_{\cP}^*(-Z) + 1/2 <W, QW> - <b,y> - <d, yI>)
%%     s.t. Z + S - QW + AE^*yE + B^* yI = C,
%%          alpha(u - y_I) = 0, S \in PSD, u >= 0, W \in Range(Q)
%%                                                   
%% P): min 1/2 <X, QX> + <C,X>, s.t., AE X = bE, B X >= d, X\in \cP
%% P'):  min 1/2 <X, QX> + <C,X>, s.t., AE X = bE, B X -alpha v = b_I, X\in \cP, v>=0.
%%*************************************************************************
%% QSDPNAL: 
%% Copyright (c) 2015 by
%% Xudong Li, Defeng Sun, and Kim-Chuan Toh
%%*************************************************************************

   function [obj,Z,W,QW,S,yE,yI,X,runhist,info] = qsdpnal_QAP(blk,Q,A,b,C,B,d,options)
%% General options
   rng('default');
   maxiter = 5000;
   stoptol = 1e-6;
   printyes = 1;
   printminoryes=1;
   startfun = @SCB_qsdp; % Algorithm QSDPNAL-Phase I for initial point
   sigma0 = 1;
   sGS.iter = 0;
   sGS.time = 0;
   use_proximal = 1;
   restart = 1;
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
   existQ = ~isempty(Q);
   existA = ~isempty(A);
   existB = ~isempty(B);
   tstart = clock;
   n = blk{1,2};
   nbar = n*(n+1)/2;
   Zero = sparse(n,n);
   Fnorm = @(X) mexFnorm(X);
   if nonnegative || ~isempty(L) || ~isempty(U) || ~isempty(fV)
      existZ = 1;
   else
      existZ = 0;
   end
   Qsolver = 'pcg'; %solve linear system corresponding to W by PCG
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
   if existB
      if isstruct(B)
         if isfield(B,'Bt'); Bt = B.Bt; end
         if isfield(B,'Bmap'); Bmap0 = B.Bmap; end
         if isfield(B,'BTmap'); BTmap0 = B.BTmap; end
      else
         Bt = B;
         Bmap0  = @(X) sdpnalAXfun(blk,Bt,X);
         BTmap0 = @(y) sdpnalAtyfun(blk,Bt,y);
      end
   else
      Bmap0 = @(X) []; BTmap0 = @(y) Zero;
   end
   if existQ
      Qmap0 = @(X) feval(Q.QXfun,X);
   else
      Qmap0 = @(X) sparse(n,n);
   end
%%
%% initial scaling
%% 
   Atorg = At;  borg = b;
   Corg = C; Lorg = L; Uorg = U; fVorg = fV;
   normCorg = 1 + Fnorm(Corg);
   normborg = 1 + norm(borg);
   if existB; 
      Btorg = Bt; dorg = d; 
      normdorg = 1 + norm(dorg); 
   else
      dorg = []; d = [];
   end
%%
%% generate initial point through QSDPNAL-Phase I (SCB_qsdp)
   sGS.iter = 0;
   sGS.time = 0;
   sGS.op.stoptol = stoptol;
   sGS.op.maxiter = 50000;
   sGS.op.lowerbound = Lorg; sGS.op.SignP = SignP;
   sGS.op.upperbound = Uorg; sGS.op.fV = fVorg;
   sGS.op.nonnegative = nonnegative;
   sGS.op.stopoption = 2;
   fprintf('\n********************************************************');
   fprintf('\n                  QSDPNAL Phase I                       ');
   fprintf('\n********************************************************');
   [obj,Z,W,QW,S,yE,yI,X,runhist_sGS,info_sGS] = startfun(blk,Q,A,borg,Corg,B,dorg,sGS.op);
   fprintf('\n -------------------------------------------------------');
   sGS.Z0 = Z; sGS.W0 = W; sGS.QW0 = QW; sGS.S0 = S; 
   sGS.yE0= yE; sGS.yI0=yI; sGS.X0=X;
   nE = info_sGS.nE; nI = info_sGS.nI;
   bscale = info_sGS.bscale; Cscale = info_sGS.Cscale;
   sigma = info_sGS.sigma;
   sGS.iter = sGS.iter + info_sGS.iter;
   sGS.time = sGS.time + info_sGS.totaltime;
   if existQ; 
      Qmap = info_sGS.Qmap; parQ = info_sGS.parQ; 
      normQorg = info_sGS.normQorg;
   else
      Qmap = Qmap0;
   end
   normA = info_sGS.normA; DA = info_sGS.DA;
   Amap = info_sGS.Amap; ATmap = info_sGS.ATmap;
   Bmap = info_sGS.Bmap; BTmap = info_sGS.BTmap;
   b = info_sGS.b; 
   C = info_sGS.C; 
   if existB; 
      normB = info_sGS.normB; DB = info_sGS.DB;
      d = info_sGS.d;
      u = info_sGS.u; v = info_sGS.v;
      alpha = 1*info_sGS.alpha;
      yI = (normB.*yI)/Cscale; 
      Bt = Bt*DB;
   end
   At = At*DA;
%% 
   X = X/bscale; yE = (normA.*yE)/Cscale; 
   Z = Z/Cscale; S = S/Cscale;
   QW = QW/Cscale; W = W/bscale; 
   objscale = bscale*Cscale;
%%-------------------------------------------------------------------
%% end of intial stage
%%-------------------------------------------------------------------
%% start QSDPNAL-Phase II
   parmain.existQ = existQ; parmain.existB = existB; parmain.existA = existA;
   parmain.bscale = bscale; parmain.Cscale = Cscale; parmain.objscale = objscale;
 
   parmain.Qsolver = Qsolver;
   parmain.DA = DA; 
   parmain.nE = nE; parmain.nI = nI; parmain.n = n;
   
   parmain.nonnegative = nonnegative; parmain.L = L; parmain.U = U;
   parmain.SignP = SignP; parmain.fV = fV;
   
   parmain.normA = normA;
   parmain.normborg = normborg;
   parmain.normCorg = normCorg;
   parmain.normQorg = normQorg;
   
   parmain.sigma = sigma;
  
   parmain.tstart= tstart;
   parmain.maxiter = maxiter;
   parmain.stoptol = stoptol;
   parmain.printyes = printyes;
   parmain.sGS = sGS;
   parmain.use_proximal = use_proximal;
   parmain.restart = restart;
   
   if existB
      parmain.alpha = alpha;
      parmain.normB = normB;
      parmain.Bt = Bt;
      parmain.normB =  normB; parmain.DB = DB;
      parmain.d =  d;
      parmain.u =  u; parmain.v = v;
      parmain.alpha = alpha;
   end
   
   parmain.info_sGS = info_sGS;
   parmain.obj = obj;
   
   parmain.parQ = parQ;
   Ainput.Amap =  Amap; Ainput.ATmap =  ATmap; Ainput.At =  At;
   Binput.Bmap =  Bmap; Binput.BTmap =  BTmap;
   
   parmain.Lorg = Lorg; parmain.Uorg = Uorg;  parmain.fVorg = fVorg;
   parmain.borg = borg; parmain.Corg = Corg;  parmain.dorg  = dorg;
   parmain.startfun = startfun;

   % main solver QSDPNAL-Phase II with restart technique
   [Z,Wvec,QW,S,yE,yI,X,info_main,runhist] = ...
    nal_qsdpQAPmain(blk,Q,Qmap,A,Ainput,B,Binput,b,C,d,parmain,Z,W,QW,S,yE,yI,X);
   iter = info_main.iter;
   msg  = info_main.msg;
   Ztmp = info_main.Ztmp;
   bscale = info_main.bscale;
   Cscale = info_main.Cscale;
   sGS = info_main.sGS;
%%-----------------------------------------------------------------
%% recover orignal variables
%%-----------------------------------------------------------------
   if (iter == maxiter)
      msg = ' maximum iteration reached';
      info.termcode = 3;
   end
   X = X*bscale; normX = norm(X,'fro');
   QX = Qmap0(X);
   yE = Cscale*(DA*yE); 
   Z = Cscale*Z; normZ = norm(Z,'fro'); Ztmp = Ztmp*bscale;
   S = Cscale*S; normS = norm(S,'fro');
   QW = Cscale*QW; Wvec = bscale*Wvec; 
   if existQ; W = smat(blk,Wvec); else W = [];end;
   if existB; 
       yI = Cscale*(DB*yI); 
       BtyI = BTmap0(yI); 
       normyI = norm(yI);
       RpI = dorg - Bmap0(X);
       dTyI = dorg'*yI;
   else
       dTyI = 0; BtyI = Zero;
   end
   AX = Amap0(X); AtyE = ATmap0(yE);
   if existA; bTyE = borg'*yE; end
   primobj = 0.5*sum(sum(X.*QX)) + sum(sum(Corg.*X));
   if existQ; WQW = sum(sum(W.*QW)); else WQW = 0; end
   dualobj = bTyE + sum(sum(Z.*Ztmp)) - 0.5*WQW + dTyI;   
   obj = [primobj,dualobj];
   if ~isempty(L)
      minLX = min(min(X(SignP.L == 1)));
   end
   if ~isempty(U)
      maxUX = max(max(X(SignP.U == 3)));
   end
   if existQ
       etaQ = Fnorm(QX - QW)/normQorg;
   else
       etaQ = 0;
   end
   if (true)
      etaSp = compnormXn(X) /(1+normX);
      etaSc = full(abs(sum(sum(X.*S))) /(1+normX +normS));
      etaS = max(etaSp,etaSc);
      projPX = proj_LU_new(X,nonnegative,SignP,Lorg,Uorg,fVorg);
      etaZp = Fnorm(X-projPX)/(1 + normX);
      ZtmX = Ztmp - X;
      normZtmX = norm(ZtmX,'fro');
      trZtmX = full(abs(sum(sum(ZtmX.*Z))));       
      etaZc = trZtmX/(1 + normZ + normZtmX);            
      etaZ = max(etaZp,etaZc);
      if existB
         etayId = norm(min(yI,0))/(1+normyI);
         normRpI = norm(RpI);
         etayIp = norm(max(RpI,0))/(1+normRpI);
         etayIc = full(abs(yI'*RpI)/(1+normyI+normRpI));
         etayI = max([etayId,etayIp,etayIc]);
      else
         etayI = 0; etayId = 0; etayIp = 0; etayIc = 0;
      end
   end
   primfeasorg = max(norm(AX-borg)/normborg,etayIp);
   dualfeasorg = max(norm(Z+S+AtyE+BtyI-QW-Corg,'fro')/normCorg,etayId);
   ttime = etime(clock,tstart);   
   runhist.nE = nE;
   runhist.ns = n;
   runhist.nI = nI;
   runhist.iter = iter;
   runhist.totaltime = ttime;
   runhist.primobjorg = primobj; 
   runhist.dualobjorg = dualobj;    
   runhist.relgaporg = (primobj-dualobj)/(1+abs(primobj)+abs(dualobj));  
   runhist.primfeasorg(iter) = primfeasorg;    
   runhist.dualfeasorg(iter) = dualfeasorg; 
   runhist.relerr = max([dualfeasorg,primfeasorg,etaS,etaZ,etayIc]);
   info.nE = nE;
   info.nI = nI;
   info.ns = n;
   info.iter = iter;
   info.totaltime = ttime;
   info.primobjorg = primobj;
   info.dualobjorg = dualobj;
   info.relgaporg = (primobj-dualobj)/(1+abs(primobj)+abs(dualobj)); 
   info.relerrS = etaS;
   info.primfeasorg = primfeasorg; info.dualfeasorg = dualfeasorg;
   info.relerrZ = etaZ;
   info.relerrQ = etaQ;
   info.etayI = etayI;
   info.etayI1 = etayId;
   info.etayI2 = etayIp;
   info.etayI3 = etayIc;
   info.relerr = max([dualfeasorg,primfeasorg,etaS,etaQ,etaZ,etayIc]);
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
	       primfeasorg,dualfeasorg,info.relgaporg);
      fprintf('\n  etaQ = %3.2e', etaQ);
      fprintf('\n  etaZ = %3.2e', etaZ);
      fprintf('\n  etaS = %3.2e', etaS);
      fprintf('\n  min(X)    = %3.2e, max(X)    = %3.2e',...
          info.minX,info.maxX); 
      fprintf('\n  total innerNT steps = %3d, total innerCG = %3d', info.totalNT, info.totalCG); 
      if ~isempty(L)
         fprintf('\n Lower bound for bounded part   = %3.2e',minLX);
      end
      if ~isempty(U)
         fprintf('\n Upper bound for bounded part   = %3.2e',maxUX);
      end
      fprintf('\n--------------------------------------------------------------');
      fprintf('------------------\n');
   end
%%**********************************************************************
    function normXn = compnormXn(X)
    
    d = eig(full(X));   
    normXn = norm(d(d<0));     
%%**********************************************************************
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
