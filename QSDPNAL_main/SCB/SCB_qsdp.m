%%*************************************************************************
%% QSDPNAL: Phase I
%% An inexact SCB based semi-Proximal ADMM for solving QSDP problems:
%% (P) min {1/2<X,QX> + <C,X> | A(X)=b, B(X)>=d, X \in K and X \in P}
%% where K = psd cone and P is a polyhedron convex set {X : L<= X <= U}
%% (D) max -(\delta_{\cP}^*(-Z) + 1/2 <W, QW> - <b,y> - <d, yI>)
%%     s.t. Z + S - QW + AE^*yE + B^* yI = C, 
%%          S \in psd, yI>=0, W \in Range(Q),
%% (D') max -(\delta_{\cP}^*(-Z) + 1/2 <W, QW> - <b,y> - <d, yI>)
%%      s.t. Z + S - QW + AE^*yE + B^* yI = C,
%%       alpha(u - y_I) = 0, S \in PSD, u >= 0, W \in Range(Q)
%% (P')  min 1/2 <X, QX> + <C,X>, 
%%       s.t., A(X) = b, B(X)-alpha v = d, 
%%             X\in psd, X\in P , v>=0.
%%*************************************************************************
%% QSDPNAL: 
%% Copyright (c) 2015 by
%% Xudong Li, Defeng Sun, and Kim-Chuan Toh
%%*************************************************************************

   function [obj,Z,W,QW,S,yE,yI,X,runhist,info] = SCB_qsdp(blk,Q,A,b,C,B,d,options,X,Z,W,S,yE,yI)
   rng('default');
%% For inexact criteria
   rhsconst  =0.09;
   feasiconst=0.9;
   tolconst  =100;
%% General options
   maxiter = 5000;
   stoptol = 1e-6;
   printyes = 1;
   printminoryes=1;
   scale_data = 1;
   rescale = 1;
   gamma = 1.618; %step length in SCB ADMM
   sigma = 1;
   gamma_reset_start = 200;
   sGS = 1;
   sig_fix = 0;
   stopoption = 0;
   
   % options for linear system solvers
   % for AEsolver: `pcg' and `sparse_chol'
   % for Qsolver and AIsolver only `pcg'
   Qsolver  = 'pcg';
   AEsolver = 'sparse_chol';%'pcg';%
   AIsolver = 'pcg';
   
   % options from users
   if isfield(options,'maxiter');  maxiter  = options.maxiter; end
   if isfield(options,'stoptol');  stoptol  = options.stoptol; end
   if isfield(options,'printyes'); printyes = options.printyes; end
   if isfield(options,'printminoryes'); printminoryes = options.printminoryes; end
   if isfield(options,'scale_data');    scale_data = options.scale_data; end
   if isfield(options,'rescale');  rescale = options.rescale; end
   if isfield(options,'sigma');    sigma = options.sigma; end
   if isfield(options,'gamma');    gamma = options.gamma; end
   if isfield(options,'stopoption'); stopoption = options.stopoption; end
%% solvers for Q, AE, AI
   if isfield(options,'Qsolver');    Qsolver  = options.Qsolver; end
   if isfield(options,'AEsolver');   AEsolver = options.AEsolver; end
   if isfield(options,'AIsolver');   AIsolver = options.AIsolver; end
   if isfield(options,'parB');  parB = options.parB; else parB = []; end
   if isfield(options,'parQ');  parQ = options.parQ; else parQ = []; end
   if isfield(options,'sGS');        sGS = options.sGS; end
%% information for \cP
   nonnegative = 0; L = []; U = []; SignP = []; fV = [];
   if isfield(options,'upperbound'); U = options.upperbound; end
   if isfield(options,'lowerbound'); L = options.lowerbound; end
   if isfield(options,'SignP');      SignP = options.SignP; end
   if isfield(options,'fV');         fV = options.fV; end
   if isfield(options,'nonnegative');  nonnegative = options.nonnegative; end
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
      Qmap0 =@(X) Zero;
   end
   %%%
   if ~exist('X','var') || ~exist('yE','var') || ~exist('S','var')
      X = Zero; S = Zero; yE = zeros(length(b),1);
   end
   if ~existA; yE = []; b = []; end
   if ~existB; d = []; Bt = []; yI = []; end
   if ~exist('W','var'); W = Zero; QW = Zero; end
   if ~exist('yI','var'); yI = zeros(length(d),1); end
   if ~exist('Z','var'); Z = Zero; end
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
      dorg = 0;
   end
   eigsopt.issym = 1;
   
   if existA
      nE = length(b);
      normA = sum(At.*At)';
      normA = max(1,sqrt(normA));
      DA = spdiags(1./normA,0,nE,nE);
      At = At*DA;
      b = DA*b;
      normAt = Fnorm(At);
      Amap  = @(X) Amap0(X)./normA;
      ATmap = @(y) ATmap0(y./normA);
      AATmap = @(y) Amap(ATmap(y));
      if strcmp(AEsolver,'sparse_chol')
         [LAAt,AAt] = mycholAAt(At,nE);
         dd = diag(AAt);
         if (Fnorm(AAt-speye(nE)) < 1e-12)
            if printyes; fprintf('\n AEAEt is identity ');  end  %% For BIQ, FAP and Theta(+) problems
            LAAt.isidentity = true;
         elseif  (Fnorm(AAt-spdiags(dd,0,nE,nE)) < 1e-15*norm(dd)) 
            fprintf('\n AAt is diagonal '); 
            LAAt.isidentity = false;
         else
            LAAt.isidentity = false;
         end      
         AEsolve = @(rhs) mylinsysolve(LAAt,rhs);
         if printyes; fprintf('\n AEAEt solved by chol decompostion'); end
      elseif strcmp(AEsolver,'pcg')
         rr = min(6,nE);
         rng('default');
         eigsopt.v0 = randn(nE,1);
         [VA,dA,flagA] = eigs(AATmap,nE,rr,'LA',eigsopt);
         dA = diag(dA); rr = sum(dA>0);
         pcgA = min(3,rr);
         parA.V = VA(:,1:pcgA); parA.Vt = parA.V';
         parA.d = dA(1:pcgA);
         parA.precond = 2;
         if printyes; fprintf('\n using pcg to solve AEAEt'); end
      end
   else
      Amap = @(X) []; ATmap = @(y) Zero;
      normA = []; DA = []; nE = 0;
   end
   if existB
      nI = length(d);
      normB = sum(Bt.*Bt)';
      normB = max(1,sqrt(normB));
      DB = spdiags(1./normB,0,nI,nI);
      Bt = Bt*DB;
      d = DB*d;
      normBt = Fnorm(Bt);
      Bmap = @(X) Bmap0(X)./normB;
      BTmap = @(y) BTmap0(y./normB);
      BBTmap = @(y) Bmap(BTmap(y));
      rr = min(6,nI);
      if isempty(parB)
         rng('default');
         eigsopt.v0 = randn(nI,1);
         [VB,dB,flagB] = eigs(BBTmap,nI,rr,'LA',eigsopt);
         dB = diag(dB); rr = sum(dB>0);
         if rr>0; alpha = max(dB)^(0.25)/2; else alpha = 1; end
      else
         alpha = parB.alpha;
      end      
      if isfield(options,'alpha'); alpha = options.alpha; end
      if strcmp(AIsolver,'pcg') && isempty(parB)
         pcgB = min(1,rr);
         parB.V = VB(:,1:pcgB);  parB.Vt = parB.V';
         parB.d = dB(1:pcgB)+alpha^2;
         parB.precond = 2;
         parB.alpha = alpha;
         if printyes; fprintf('\n using pcg to solve AIAIt'); end
      end    
   else
      Bmap = @(X) []; BTmap = @(y) Zero; BBTmap = @(y) [];
      normB = []; DB = []; nI = 0; normBt = []; alpha = [];
   end
   
   if (scale_data == 1)
       if existB; bscale = max(1,norm([b;d])); else; bscale = max(1,norm(b)); end
       Cscale = max(1,Fnorm(C));
   elseif scale_data == 0
       bscale = 1;
       Cscale = 1;
   end
   b = b/bscale; d  = d/bscale; 
   X = X/bscale; L = L/bscale; U = U/bscale; fV = fV/bscale;
   C = C/Cscale;
   yE = (normA.*yE)/Cscale; Z = Z/Cscale; S = S/Cscale; 
   W = W/bscale;
   yI = (normB.*yI)/Cscale;
   objscale = bscale*Cscale;
   
   if printyes
        fprintf('\n *******************************************************');
        fprintf('******************************************');
        fprintf('\n \t\t   SCB_qsdp with  beginning gamma = %6.3f', gamma);
        fprintf('\n ******************************************************');
        fprintf('*******************************************\n');
        if printminoryes
            fprintf('\n n = %3.0f, mE = %3.0f',n, nE);
            if existB; fprintf( ' mI = %3.0d', nI); end
            fprintf('\n scale_data = %2.0f', scale_data);
            fprintf('\n bscale = %3.2e, cscale = %3.2e', bscale, Cscale);
            fprintf('\n ---------------------------------------------------');
        end
        fprintf('\n  iter|  pinforg  dinforg   relgaporg|   pobj         dobj    |');
        fprintf(' time |   sigma |rankS, minspeed|gamma|');
   end
%%-------------------------------------------------------------------
%% start main solver
%  parameters for main solver
   parmain.existQ = existQ; parmain.existB = existB; parmain.existA = existA;
   parmain.existZ = existZ;
   parmain.bscale = bscale; parmain.Cscale = Cscale; parmain.objscale = objscale;
   parmain.AEsolver = AEsolver; parmain.AIsolver = AIsolver; parmain.Qsolver = Qsolver;
   parmain.DA = DA; parmain.DB = DB;
   parmain.nE = nE; parmain.nI = nI; parmain.n = n;
   
   parmain.nonnegative = nonnegative; parmain.L = L; parmain.U = U;
   parmain.SignP = SignP; parmain.fV = fV;
   
   parmain.normA = normA; parmain.normB = normB;
   parmain.normAt = normAt;
   parmain.normBt = normBt;
   parmain.normborg = normborg;
   parmain.normCorg = normCorg;
   
   parmain.sigma = sigma;
   parmain.alpha = alpha;
   parmain.gamma = gamma;
   parmain.tstart= tstart;
   parmain.maxiter = maxiter;
   parmain.rescale = rescale;
   parmain.stoptol = stoptol;
   parmain.printyes = printyes;
   parmain.sGS = sGS;
   parmain.gamma_reset_start = gamma_reset_start;
   parmain.sig_fix = sig_fix;
   parmain.stopoption = stopoption;
   
   parmain.feasiconst = feasiconst;
   parmain.rhsconst = rhsconst;
   parmain.tolconst = tolconst;
   
   if strcmp(AEsolver,'pcg'); 
      parmain.parA = parA; 
   elseif strcmp(AEsolver,'sparse_chol');
      parmain.AEsolve = AEsolve;
   end
   parmain.parB = parB; parmain.parQ = parQ;
   
   Ainput.Amap =  Amap; Ainput.ATmap =  ATmap; Ainput.AATmap =  AATmap;
   Binput.Bmap =  Bmap; Binput.BTmap =  BTmap; Binput.BBTmap =  BBTmap; 

   % main file of SCB ADMM for solving QSDP problem
   [Z,Wvec,QW,S,yE,yI,X,info_main,runhist] = ...
    SCB_qsdpmain(blk,Qmap0,Ainput,Binput,b,C,d,parmain,Z,W,S,yE,yI,X);
   
   % information 
   iter = info_main.iter;
   bscale = info_main.bscale;
   Cscale = info_main.Cscale;
   breakyes = info_main.breakyes;
   msg = info_main.msg;
   iterProj = info_main.iterProj;
   repeatyE = info_main.repeatyE;
   repeatyI = info_main.repeatyI;
   repeatQ  = info_main.repeatQ;
   Ztmp = info_main.Ztmp;
   if existQ
      normQorg = info_main.normQorg;
      Qmap = info_main.Qmap;
      tolQ = info_main.tolQ;
      parQ = info_main.parQ;
   end
   if existB
      u = info_main.u; v = info_main.v;
      d = info_main.d;  normd = info_main.normd;
   end
    b = info_main.b;  normb = info_main.normb;  
    C = info_main.C;  normC = info_main.normC;
    sigma = info_main.sigma;
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
   QW = Cscale*QW;
   Wvec = bscale*Wvec; 
   if existQ; W = smat(blk,Wvec); else W = []; end
   if existB; 
       yI = Cscale*(DB*yI); 
       BtyI = BTmap0(yI); 
       normyI = norm(yI);
       RpI = dorg - Bmap0(X);
       dTyI = dorg'*yI;
   else
       dTyI = 0; BtyI = Zero;
   end
   if existA; bTyE = borg'*yE; end
   AX = Amap0(X); AtyE = ATmap0(yE);
   primobj = full(0.5*sum(sum(X.*QX)) + sum(sum(Corg.*X)));
   if existQ; WQW = sum(sum(W.*QW)); else WQW = 0; end
   dualobj = full(bTyE + sum(sum(Z.*Ztmp)) - 0.5*WQW + dTyI);   
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
   etab = norm(AX-borg)/normborg;
   primfeasorg = max(etab,etayIp);
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
   runhist.relerr = max([dualfeasorg,primfeasorg,etaS,etaQ,etaZ,etayIc]);
   info.nE = nE;
   info.nI = nI;
   info.ns = n;
   info.iter = iter;
   info.totaltime = ttime;
   info.primobjorg = primobj;
   info.dualobjorg = dualobj;
   info.relgaporg = (primobj-dualobj)/(1+abs(primobj)+abs(dualobj)); 
   info.etaS = etaS;
   info.primfeasorg = primfeasorg; info.dualfeasorg = dualfeasorg;
   info.etaZ = etaZ;
   info.etaQ = etaQ;
   info.etayI = etayI;
   info.etayId = etayId;
   info.etayIp = etayIp;
   info.etayIc = etayIc;
   info.etab = etab;
   info.relerr = max([dualfeasorg,primfeasorg,etaS,etaQ,etaZ,etayIc]);
   info.minX = min(min(X));
   info.maxX = max(max(X));
   info.psqmrQ1 = sum(runhist.psqmrQ1);
   info.psqmrQ2 = sum(runhist.psqmrQ2);
   info.psqmrA1 = sum(runhist.psqmrA1);
   info.psqmrA2 = sum(runhist.psqmrA2);
   info.psqmrB1 = sum(runhist.psqmrB1);
   info.psqmrB2 = sum(runhist.psqmrB2);
   info.breakyes = breakyes;
%% added for Phase II
   info.Amap = Amap; info.ATmap = ATmap;
   info.Bmap = Bmap; info.BTmap = BTmap;
   info.bscale = bscale; info.Cscale = Cscale;
   info.sigma = sigma;
   info.Qmap = Qmap;
   info.AAt = AAt;
   info.AEsolve = AEsolve;
   if exist('At','var'); info.At = At; end
   if exist('Bt','var'); info.Bt = Bt; end
   if exist('parA','var');info.parA = parA; end
   if existQ; 
      info.tolQ = tolQ; 
      info.normQorg = normQorg; 
      info.parQ = parQ; 
   end  
   if existB; 
      info.u = u; info.v = v; info.DB = DB;
      info.normB = normB; info.parB = parB;
      info.normBt = normBt;
      info.d = d; info.normd = normd;
      info.alpha = alpha;
      %info.tolyI = tolyI;
      info.BBTmap = BBTmap;
   end
   info.normA = normA; info.DA = DA; info.normAt = normAt;
   info.b = b; info.normb = normb; info.C = C; info.normC = normC;
   info.Ztmp = Ztmp;
   %for output infomation
   info.out.relerr = info.relerr;
   info.out.nE = nE;
   info.out.nI = nI;
   info.out.ns = n;
   info.out.iter = iter;
   info.out.totaltime = ttime;
   info.out.primobjorg = primobj;
   info.out.dualobjorg = dualobj;
   info.out.relerrZ = etaZ;
   info.out.relerrQ = etaQ;
   info.out.relerrS = etaS;
   info.out.etayI = etayI;
   info.out.etayId = etayId;
   info.out.etayIp = etayIp;
   info.out.etayIc = etayIc;
   info.out.relgaporg = info.relgaporg; 
   info.out.primfeasorg = primfeasorg; 
   info.out.dualfeasorg = dualfeasorg;
   info.out.minX = min(min(X));
   info.out.maxX = max(max(X));
   
   if (printminoryes) 
      if ~isempty(msg); fprintf('\n %s',msg); end
      fprintf('\n--------------------------------------------------------------');
      fprintf('------------------');
      fprintf('\n  number iter = %2.0d',iter);      
      fprintf('\n  number of iterProj = %2.0f',iterProj); 
      fprintf('\n  number iter = %2.0d \t repeatQW = %d \t repeaty= %d \t repeatyi = %d',...
         iter,repeatQ,repeatyE,repeatyI);
      if existQ && strcmp(Qsolver,'pcg')
         fprintf('\n  number of PSQMRiterQ for 1st system = %2.0f',info.psqmrQ1);
         fprintf('\n  number of PSQMRiterQ for 2nd system = %2.0f',info.psqmrQ2);
      end
      if existA && strcmp(AEsolver,'pcg')
         fprintf('\n  number of PSQMRiterA for 1st system = %2.0f',info.psqmrA1);
         fprintf('\n  number of PSQMRiterA for 2nd system = %2.0f',info.psqmrA2);
      end
      if existB
         fprintf('\n  number of PSQMRiterB for 1st system = %2.0f',info.psqmrB1);
         fprintf('\n  number of PSQMRiterB for 2nd system = %2.0f',info.psqmrB2);
      end
      fprintf('\n  time = %3.2f',ttime);       
      fprintf('\n  time per iter = %5.4f',ttime/iter);       
      fprintf('\n  primobj = %9.8e',runhist.primobjorg);       
      fprintf('\n  dualobj = %9.8e',runhist.dualobjorg);            
      fprintf('\n  primfeas    = %3.2e, dualfeas    = %3.2e, relgap    = %3.2e',...
	      runhist.primfeas(end),runhist.dualfeas(end),runhist.relgap(end)); 
      fprintf('\n  primfeasorg = %3.2e, dualfeasorg = %3.2e, relgaporg = %3.2e',...
	      runhist.primfeasorg(end),runhist.dualfeasorg(end),runhist.relgaporg(end));
      fprintf('\n  etaQ = %3.2e', etaQ); fprintf('  etaS = %3.2e', etaS);
      fprintf('\n  etaZ = %3.2e', etaZ);
      fprintf('\n  etab = %3.2e', etab);
      if existB
          fprintf('\n etayIp = %3.2e, etayId = %3.2e, etayIc = %3.2e', etayIp, etayId, etayIc);
      end
      fprintf('\n  min(X)    = %3.2e, max(X)    = %3.2e',...
          info.minX,info.maxX); 
      if ~isempty(L)
         fprintf('\n Lower bound for bounded part   = %3.2e',minLX);
      end
      if ~isempty(U)
         fprintf('\n Upper bound for bounded part   = %3.2e',maxUX);
      end
      fprintf('\n sigma = %3.1e',sigma*bscale/Cscale);
      fprintf('\n--------------------------------------------------------------');
      fprintf('------------------\n');
   end
%%**********************************************************************
%% solve A*q = r.
%%**********************************************************************
    function q = mylinsysolve(L,r) 

    if (L.isidentity)
       q = r; 
    else
       if strcmp(L.matfct_options,'chol')
          q(L.perm,1) = mextriang(L.R, mextriang(L.R,r(L.perm),2) ,1);
       elseif strcmp(L.matfct_options,'spcholmatlab')
          q(L.perm,1) = mexbwsolve(L.Rt,mexfwsolve(L.R,r(L.perm,1)));
       end
    end
%%**********************************************************************
%%**********************************************************************
    function normXn = compnormXn(X)
    
    d = eig(full(X));   
    normXn = norm(d(d<0));     
%%**********************************************************************
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
