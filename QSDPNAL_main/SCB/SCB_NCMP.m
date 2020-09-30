%%*************************************************************************
%% QSDPNAL: Phase I
%% An inexact SCB based semi-Proximal ADMM for solving NCM problems with 
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

   function [obj,Z,Xi,S,yE,X,runhist,info] = SCB_NCMP(blk,H,A,b,C,options,X,Z,Xi,S,yE)
%% For inexact criteria
   rng('default');
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
   phase2 = 0;
   % options for linear system solvers
   % for AEsolver: `pcg' and `sparse_chol'
   AEsolver = 'sparse_chol';%'pcg';%

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
   if isfield(options,'phase2'); phase2 = options.phase2; end
%% solvers for AE
   if isfield(options,'AEsolver');   AEsolver = options.AEsolver; end
   if isfield(options,'sGS');        sGS = options.sGS; end
%% information for \cP
   nonnegative = 0; L = []; U = []; SignP = []; fV = [];
   if isfield(options,'upperbound'); U = options.upperbound; end
   if isfield(options,'lowerbound'); L = options.lowerbound; end
   if isfield(options,'SignP');      SignP = options.SignP; end
   if isfield(options,'fV');         fV = options.fV; end
   if isfield(options,'nonnegative');  nonnegative = options.nonnegative; end
%% problem property
   existA = ~isempty(A);
   tstart = clock;
   n = blk{1,2};
   nbar = n*(n+1)/2;
   Zero = sparse(n,n);
   Fnorm = @(X) mexFnorm(X);
%% AEmap, AETmap
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
   %%%
   if ~exist('X','var'); X = Zero; end
   if ~exist('yE','var'); yE = zeros(length(b),1); end
   if ~exist('S','var'); S = Zero; end
   if ~existA; yE = []; b = []; end
   if ~exist('Xi','var'); Xi = Zero; end
   if ~exist('Z','var'); Z = Zero; end
%%
%% initial scaling
%% 
   Atorg = At;  borg = b; Horg = H;
   Corg = C; Lorg = L; Uorg = U; fVorg = fV;
   normCorg = 1 + Fnorm(Corg);
   normborg = 1 + norm(borg);
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
   
   if (scale_data == 1)
       bscale = max(1,norm(b));
       Cscale = max(1,Fnorm(C));
   elseif scale_data == 0
       bscale = 1;
       Cscale = 1;
   end
   b = b/bscale; C = C/Cscale;
   X = X/bscale; Z = Z/Cscale;
   L = L/bscale; U = U/bscale; fV = fV/bscale;
   yE = (normA.*yE)/Cscale; 
   S = S/Cscale;
   objscale = bscale*Cscale;
   Xi = Xi/sqrt(objscale);
   
   if printyes
        fprintf('\n *******************************************************');
        fprintf('******************************************');
        fprintf('\n \t\t   SCB_NCMP with  beginning gamma = %6.3f', gamma);
        fprintf('\n ******************************************************');
        fprintf('*******************************************\n');
        if printminoryes
            fprintf('\n n = %3.0f, mE = %3.0f',n, nE);
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
   parmain.existA = existA;
   parmain.bscale = bscale; parmain.Cscale = Cscale; parmain.objscale = objscale;
   parmain.AEsolver = AEsolver;   
   parmain.DA = DA;  
   parmain.nE = nE;   parmain.n = n;
   
   parmain.nonnegative = nonnegative; parmain.L = L; parmain.U = U;
   parmain.SignP = SignP; parmain.fV = fV;
   
   parmain.normA = normA;  
   parmain.normAt = normAt;
   parmain.normborg = normborg;
   parmain.normCorg = normCorg;
   
   parmain.sigma = sigma;
   parmain.gamma = gamma;
   parmain.tstart= tstart;
   parmain.maxiter = maxiter;
   parmain.rescale = rescale;
   parmain.feasiconst = feasiconst;
   parmain.rhsconst = rhsconst;
   parmain.stoptol = stoptol;
   parmain.printyes = printyes;
   parmain.sGS = sGS;
   parmain.gamma_reset_start = gamma_reset_start;
   parmain.sig_fix = sig_fix;
   parmain.stopoption = stopoption;
   parmain.tolconst = tolconst;
   
   if strcmp(AEsolver,'pcg'); 
      parmain.parA = parA; 
   elseif strcmp(AEsolver,'sparse_chol');
      parmain.AEsolve = AEsolve;
   end
   
   Ainput.Amap =  Amap; Ainput.ATmap =  ATmap; Ainput.AATmap =  AATmap;
   % main file of SCB ADMM for solving NCM problem with \cP
   [X,Z,yE,S,Xi,info_main,runhist] = ...
     SCB_NCMPmain(Horg,Ainput,b,C,parmain,X,Z,yE,S,Xi);
  
   % information
   iter = info_main.iter;
   msg  = info_main.msg;
   bscale = info_main.bscale;
   Cscale = info_main.Cscale;
   breakyes = info_main.breakyes;
   iterProj = info_main.iterProj;
   repeatyE = info_main.repeatyE;
   repeatQ  = info_main.repeatQ;
   sigma = info_main.sigma;
   Ztmp = info_main.Ztmp;
   if phase2
      b = info_main.b;  normb = info_main.normb;  
      C = info_main.C;  normC = info_main.normC;
      H = info_main.H;
   end
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
   Xi = sqrt(bscale*Cscale)*Xi;
   HXi = Horg.*Xi;
   if existA; bTyE = borg'*yE; end
   AX = Amap0(X); AtyE = ATmap0(yE);
   primobj = 0.5*Fnorm(HX)^2 + sum(sum(Corg.*X));
   dualobj = sum(sum(Z.*Ztmp)) + bTyE - 0.5*Fnorm(Xi)^2;   
   obj = [primobj,dualobj];
   etaQ = Fnorm(Horg.*(Xi+HX))/(1+Fnorm(Horg)^2);
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
   end
   etab = norm(AX-borg)/normborg;
   primfeasorg = max(etab);
   dualfeasorg = max(norm(Z+HXi+S+AtyE-Corg,'fro')/normCorg);
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
   runhist.relerr = max([dualfeasorg,primfeasorg,etaS,etaQ]);
   info.nE = nE;
   info.ns = n;
   info.iter = iter;
   info.totaltime = ttime;
   info.primobjorg = primobj;
   info.dualobjorg = dualobj;
   info.relgaporg = (primobj-dualobj)/(1+abs(primobj)+abs(dualobj));
   info.etaZ = etaZ;
   info.etaS = etaS;
   info.primfeasorg = primfeasorg; info.dualfeasorg = dualfeasorg;
   info.etaQ = etaQ;
   info.etab = etab;
   info.relerr = max([dualfeasorg,primfeasorg,etaS,etaQ]);
   info.minX = min(min(X));
   info.maxX = max(max(X));
   info.psqmrA1 = sum(runhist.psqmrA1);
   info.psqmrA2 = sum(runhist.psqmrA2);
%% added for Phase II
   if phase2
      info.nXi = nbar;
      info.Amap = Amap; info.ATmap = ATmap;
      info.bscale = bscale; info.Cscale = Cscale;
      info.sigma = sigma;
      info.AAt = AAt;
      info.AEsolve = AEsolve;
      if exist('At','var'); info.At = At; end
      info.normA = normA; info.DA = DA; info.normAt = normAt;
      info.b = b; info.normb = normb; info.C = C; info.normC = normC;
      info.H = H; info.Ztmp = Ztmp;
      info.breakyes = breakyes;
   end
   if (printminoryes) 
      if ~isempty(msg); fprintf('\n %s',msg); end
      fprintf('\n--------------------------------------------------------------');
      fprintf('------------------');
      fprintf('\n  number iter = %2.0d',iter);      
      fprintf('\n  number of iterProj = %2.0f',iterProj); 
      fprintf('\n  number iter = %2.0d \t repeaty= %d ',...
         iter,repeatyE);
      fprintf('\t repeatQ= %d ',...
         iter,repeatQ);
      if existA && strcmp(AEsolver,'pcg')
         fprintf('\n  number of PSQMRiterA for 1st system = %2.0f',info.psqmrA1);
         fprintf('\n  number of PSQMRiterA for 2nd system = %2.0f',info.psqmrA2);
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
      fprintf('\n  etab = %3.2e', etab); fprintf('  etaZ = %3.2e', etaZ);
      fprintf('\n  min(X)    = %3.2e, max(X)    = %3.2e',...
          info.minX,info.maxX); 
      fprintf('\n sigma = %3.1e',sigma);
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
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
