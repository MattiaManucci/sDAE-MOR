% Fig3Fig4.m
% Purpose:  Reproduce Figure 3 and Figure 4 from the manuscript [1]
% Author:   Mattia Manucci
% Date:     2026-05-06
%
% ------------------------------------------------------------------------------
% Bibliography
% [1] M. Manucci, B. Unger. Balancing-Based model reduction for switched
% descriptor systems, arXiv XXXXXX, 2026
% ------------------------------------------------------------------------------
% Summary of operations:
%  1) Load the data
%  2) Preprocess / compute quantities to plot: solve M pairs of generalized Lyapunov
%  equations (GLEs) where M are the number of modes of the switched DAE
%  3) Compute solutions of the FOM and ROM for two different input
%  functions and switching signals
%  4) Create Figure 3 
%  5) Create Figure 4 
%
% Notes:
%  - This script runs in the base workspace (no input/output arguments).
%  - Modify paths, filenames, and parameters below as needed.
clc
clearvars
close all
addpath(genpath('./Functions/'))
%-------------------------------------------------------------
FS = 15;       % Fontsize
FN = 'times';  % Fontname
LW = 2;      % Linewidth
MS = 7.8;    % Markersize
%% 1) Define parameters and load data
tol=1e-10;
% Constrained mass-spring-dumper system
g=579; % number of masses 
n=2*g+1; % dimension of the system
nf(1)=n-3; % dimension of the differential part
DAE_Index(1)=3; % DAE index
dimKerE(1)=1; % dimension of the kernel of E
mas=ones(g,1);
k1=1.5*ones(g-1,1); k2=2*ones(g,1);
d1=0.7*ones(g-1,1); d2=0.7*ones(g,1);
[E1, A1, B1, C1, M1, D1, K1, G1] = msd_ind3(g, mas, k1, k2, d1, d2);
A{1}=A1; E{1}=E1; B{1}=B1; C{1}=C1;
si(1)=size(B1,2); so(1)=size(C1,1);
% Semi-discretized Stokes problem
nx = 20; ny = 20;
DAE_Index(2)=2;
n = (nx-1)*ny+(ny-1)*nx+nx*ny-1; I=speye(n);
dimKerE(2)=nx*ny-1; opts=[];
[E1, A1, B1, C1, nf(2)] = stokes_ind2(si(1), so(1), nx, ny,opts);
A{2}=A1; E{2}=E1; B{2}=B1/normest(B1); C{2}=C1; BStokesNorm=normest(B1);
%% 2) Preprocess / compute quantities to plot: decouple the DAE modes and solve M pairs of generalized Lyapunov equations (GLEs) where M are the number of modes of the switched DAE
%% Decoupling impulsive and differential components through Quasi-Weierstrass form
nm = 2;
V=cell(nm,1); W=cell(nm,1);
T=cell(nm,1); invS=cell(nm,1);
DS=cell(nm,1); IS=cell(nm,1);
n1=cell(nm,1); Iimp=cell(nm,1);
Iaug=cell(nm,1);
A_decoupled=cell(nm,1); E_decoupled=cell(nm,1); J=cell(nm,1); 
Bdiff=cell(nm,1); Cdiff=cell(nm,1);
Map = cell(nm,nm); Bimp = cell(nm,nm); Cimp = cell(nm,nm); 
N = cell(nm,1); Nnn = cell(nm,1); NnnID = cell(nm,1); time_QWF = zeros(nm,1);
for i=1:nm
    tic
    [V{i},W{i}] = QWF(A{i},E{i},I(:,(n-dimKerE(i)+1):end),DAE_Index(i));
    time_QWF(i)=toc;
    T{i}=[V{i},W{i}]; T{i}(abs(T{i})<1e-12)=0; T{i}=sparse(T{i});
    invS{i}=([E{i}*V{i},A{i}*W{i}]); invS{i}(abs(invS{i})<1e-12)=0; invS{i}=sparse(invS{i});
    n1{i}=size(V{i},2);

    A_decoupled{i} = invS{i}\(A{i}*T{i}); A_decoupled{i}(abs(A_decoupled{i})<2e-12)=0;
    E_decoupled{i} = invS{i}\(E{i}*T{i}); E_decoupled{i}(abs(E_decoupled{i})<2e-12)=0;
    J{i} = A_decoupled{i}(1:nf(i),1:nf(i));
    N{i} = E_decoupled{i}((nf(i)+1):end,(nf(i)+1):end);
    Nnn{i} = sparse(n,n); Nnn{i}((nf(i)+1):end,(nf(i)+1):end) =  N{i}; 
    NnnID{i} = sparse(n,n); NnnID{i}((nf(i)+1):end,(nf(i)+1):end) =  N{i}^0;
    for j=1:(DAE_Index(i)-1)
        Nnn{i,j} = sparse(n,n); Nnn{i,j}((nf(i)+1):end,(nf(i)+1):end) =  N{i}^(j); 
    end
    Iaug{i} = sparse(nf(i),n); Iaug{i}(1:(size(Iaug{i},1)+1):nf(i)^2) = 1;
    Bdiff{i} = Iaug{i}*(invS{i}\B{i});  Cdiff{i} = C{i}*T{i}*Iaug{i}'; 

end
% Impulsive components
for i=1:nm
  for j=1:nm
      if i==j
         Map{i,j} = speye(nf(i));
         Bimp{i,j} = sparse(nf(i),si);
         Cimp{i,j} = sparse(so,nf(i));
      else
         Map{i,j} = Iaug{i}*(T{i}\(T{j}*Iaug{j}'));
         Bimp{i,j} =  (NnnID{j})*(invS{j}\B{j});
         for k = 1:(DAE_Index(j)-1)
             Bimp{i,j} = [Bimp{i,j},Nnn{j,k}*(invS{j}\B{j})];
         end
         Bimp{i,j} = Iaug{i}*(T{i}\Bimp{i,j});
         Cimp{i,j} = C{i}*T{i}*Nnn{i,1}*(T{i}\(T{j}*(Iaug{j}')));
         for k = 2:(DAE_Index(i)-1)
             Cimp{i,j} = [Cimp{i,j};C{i}*T{i}*Nnn{i,k}*(T{i}\(T{j}*(Iaug{j}')))];
         end
      end
   end
end
%% Approximating the GLEs solutions
opts_GLE.tol = tol; 
opts_GLE.tol2 = tol;
options.problem = 1;
opts_GLE.RelativeErr = 1;
if options.problem==0
    opts_GLE.rescale = 0;
else
    opts_GLE.rescale = 1;
end
Z = cell(nm,1); S = cell(nm,1); 
Z_m = cell(nm,1); S_m = cell(nm,1);
rho_P = zeros(nm,1); rho_Q = zeros(nm,1); JJ = J;
LMI_Phat_Z = zeros(nm,2);  LMI_Phat_S = zeros(nm,2);
Z2 = cell(nm,1); S2 = cell(nm,1); 
Z2_m = cell(nm,1); S2_m = cell(nm,1);
rho_P2 = zeros(nm,1); rho_Q2 = zeros(nm,1);
LMI_Phat_Z2 = zeros(nm,2);  LMI_Phat_S2 = zeros(nm,2);
BBimp=[] ;CCimp=[]; EE = cell(nm,1); JJt = cell(nm,1);
sigma_min = [0.4725,16.24]; %Precomputed values of sigma_min of the Lyap operator
for i=1:nm
    for j=1:nm
        EE{j} = speye(nf(i),nf(i));
    end
    for j = 1:nm
        if j~=i
            JJ{j} = Map{i,j}*J{j}*Map{j,i}; JJ{j}((abs(JJ{j})<1e-12))=0;
        end
    end
    J_c=JJ{i};
    JJ{i}=JJ{1};
    JJ{1}=J_c;
    
    opts_GLE.sigma_min = sigma_min(i);
    
    % For differential input-output
    for j=1:nm
        for k=1:si
            if j==i
                BB(:,(j-1)*si+k) = Bdiff{j}(:,k);
            else
                BB(:,(j-1)*si+k) = Map{i,j}*Bdiff{j}(:,k);
            end
        end
        for k=1:so
            if j==i
                CC((j-1)*so+k,:) = Cdiff{j}(k,:);
            else
                CC((j-1)*so+k,:) = Cdiff{j}(k,:)*Map{j,i};
            end
        end
    end
    % For impulsive input-output
    %------------------------------
    for j=1:nm
        BBimp = [BBimp,Bimp{i,j}]; %Index before is i-j
        CCimp = [CCimp;Cimp{j,i}];
    end
    %------------------------------
    BB_2 = [BB,BBimp]; CC_2 = [CC;CCimp];
    % Approximate Solution of GLE for P
    [Z{i},rho_P(i)]=Solve_LSS_GLE_2(JJ,EE,BB_2,opts_GLE,1);
    Z_m{i} = real(Z_m{i}); %Correcting some numerical errors that originates small complex entries
    % Approximate Solution of GLE for P without impulsive components on BB
    [Z2{i},rho_P2(i)]=Solve_LSS_GLE_2(JJ,EE,BB,opts_GLE,1);
    Z2_m{i} = real(Z2_m{i}); %Correcting some numerical errors that originates small complex entries
    % Enforcing LMI for P
    [Z_m{i},LMI_Phat_Z(i,:)] = enforcing_LMI_sDAE(JJ,EE{i},Z{i},BB_2,rho_P(i)); 
    Z_m{i} = real(Z_m{i}); 
    [Z2_m{i},LMI_Phat_Z2(i,:)] = enforcing_LMI_sDAE(JJ,EE{i},Z2{i},BB,rho_P2(i));
    Z2_m{i} = real(Z2_m{i}); 
    % Approximate Solution of GLE for Q
    [S{i},rho_Q(i)]=Solve_LSS_GLE_2(JJ,EE,CC_2,opts_GLE,2);
    % Approximate Solution of GLE for Q without impulsive components on CC
    [S2{i},rho_Q2(i)]=Solve_LSS_GLE_2(JJ,EE,CC,opts_GLE,2);
    for j=1:nm
       JJt{j} = JJ{j}';
    end
    % Enforcing LMI for Q
    [S_m{i},LMI_Phat_S(i,:)] = enforcing_LMI_sDAE(JJt,EE{i},S{i},CC_2',rho_Q(i));
    [S2_m{i},LMI_Phat_S2(i,:)] = enforcing_LMI_sDAE(JJt,EE{i},S2{i},CC',rho_Q2(i));
    JJ = J;
    clear BB  CC
    BBimp=[] ; CCimp=[];
end
%% Mass Matrix
for j=1:nm
        EE{j} = speye(nf(j),nf(j));
end
%% 3) Testing the ROM: first switching sequence and input function
rv=1:1:30; %You can also inset a vector here
nr=numel(rv);
ns=40; % Up to 40 switching times
ns=ns+1;
seed = 123; rng(seed); % for the first switching path
% Choose random order of switches
rr = 1:ns; rr=mod(rr,nm+1); rr(rr==0)=nm;
random_indices = randperm(length(rr));
rr = rr(random_indices); 
% Time scale and time step
scale=0.1; dt=scale*1e-1;
% First input function
u=@(t) 1*BStokesNorm;
u_ex = cell(nm,1);
u_ex{1} = @(t) BStokesNorm*[1;0;0];
u_ex{2} = @(t) BStokesNorm*[1;0];
%% Evaluate the FOM
dtt=cell(ns,1); t=zeros(ns,1); tspan=cell(ns,1); y=cell(ns,1); x=cell(ns,1); imp_input = 0;
tic
for i=1:ns
    if i==1
        x_0=zeros(nf(rr(i)),1);
        t(i)=scale*randi([1 10],1,1);
        dtt{i}=0:dt:t(i);
    else
        t(i)=t(i-1)+scale*randi([1 10],1,1);
        dtt{i}=t(i-1):dt:t(i);
        if rr(i)==rr(i-1)
            x_0=Map{rr(i),rr(i-1)}*x{i-1}(end,:)';
        else
            imp_input = imp_input + norm(u_ex{rr(i-1)}(t(i-1)))^2;
            x_0=Map{rr(i),rr(i-1)}*x{i-1}(end,:)'+Bimp{rr(i),rr(i-1)}*u_ex{rr(i-1)}(t(i-1));
        end
    end
    odefun=@(t,y) (J{rr(i)}*y)+Bdiff{rr(i)}*u(t);

    options = odeset('RelTol',max([tol*1e-2,1e-12]),'AbsTol',max([tol*1e-2,1e-12]),'Mass',EE{rr(i)});
    [tspan{i},x{i}]=ode15s(odefun,dtt{i},x_0,options);
    if i<ns
         y{i}=[Cdiff{rr(i)}*(x{i}');Cimp{rr(i+1),rr(i)}*(x{i}')];
    else
         y{i}=[Cdiff{rr(i)}*(x{i}');Cimp{rr(i),rr(i)}*(x{i}')];
    end
end
time_FOM(1)=toc;
%% Evaluate the ROM
yred=cell(ns,nr); Xplus = cell(nr,ns); Xminus = cell(nr,ns);
xred = cell(nr,ns); MAX_EIG_c = zeros(nm,nr); MAX_EIG_o = zeros(nm,nr);
time_ROM = zeros(nr,1);
for jj=1:nr
    r=rv(jj);
    % Reduced Operators
    Jred=cell(nm,1); Bred=cell(nm,1); Cred=cell(nm,1); Ered=cell(nm,1); BimpRed = cell(nm,nm); CimpRed = cell(nm,nm);
    U=cell(nm,1); H=cell(nm,1); Ss=cell(nm,1); Vs=cell(nm,1); VV=cell(nm,1); W=cell(nm,1);
    r=r*ones(nm,1);
    for i=1:nm
        
        H{i}=Z_m{i}'*S_m{i}; [U{i},Ss{i},Vs{i}]=svd(H{i}); Ss{i}=diag(Ss{i})+tol;
        % Projection Matrices VV and W
        if (r(i)>size(Vs{i},2))||(r(i)>size(U{i},2))
           r(i)=min(size(Vs{i},2),size(U{i},2));
        end
        VV{i}=Z_m{i}*U{i}(:,1:r(i))*diag((Ss{i}(1:r(i))).^(-0.5));
        W{i}=S_m{i}*(Vs{i}(:,1:r(i)))*diag((Ss{i}(1:r(i))).^(-0.5)); 
        Jred{i}=(W{i}'*((J{i}*VV{i})));
        G = eig(Jred{i}); G=max(G); 
        if G>0
            disp(G)
        end
        Ered{i}=(W{i}'*((EE{i}*VV{i})));
        Bred{i}=(W{i}'*Bdiff{i});
        Cred{i}=Cdiff{i}*VV{i};
        for j = 1:nm
            BimpRed{i,j} = W{i}'*Bimp{i,j};
            CimpRed{j,i} = Cimp{j,i}*VV{i};
        end

        % Store conditions
        Lyap_c = Jred{i}*diag((Ss{i}(1:r(i))))+diag((Ss{i}(1:r(i))))*(Jred{i}') + Bred{i}*Bred{i}';
        Lyap_o = (Jred{i}')*diag((Ss{i}(1:r(i))))+diag((Ss{i}(1:r(i))))*Jred{i} + Cred{i}'*Cred{i};
       
        MAX_EIG_c(i,jj) = max(eig(Lyap_c)); MAX_EIG_o(i,jj) = max(eig(Lyap_o));
        
    end
    % Reduced solution size r
    tic
    for i=1:ns

        if i==1
            x_0=zeros(r(rr(i)),1);
        else
            if rr(i)==rr(i-1)
                x_0=(W{rr(i)}'*Map{rr(i),rr(i-1)}*VV{rr(i-1)})*(xred{jj,i-1}(end,:)');
            else
                x_0=(W{rr(i)}'*Map{rr(i),rr(i-1)}*VV{rr(i-1)})*(xred{jj,i-1}(end,:)')+BimpRed{rr(i),rr(i-1)}*u_ex{rr(i-1)}(t(i-1));
            end
            Xminus{jj,i-1} = xred{jj,i-1}(end,:)';
            Xplus{jj,i-1} = x_0;
        end
        odefun=@(t,y) Jred{rr(i)}*y+Bred{rr(i)}*u(t);
        options = odeset('RelTol',max([tol*1e-2,1e-12]),'AbsTol',max([tol*1e-2,1e-12]),'Mass',Ered{rr(i)});
        [~,xred{jj,i}]=ode15s(odefun,dtt{i},x_0,options);
        if i<ns
            yred{i,jj}=[Cred{rr(i)}*(xred{jj,i}');CimpRed{rr(i+1),rr(i)}*(xred{jj,i}')];
        else
            yred{i,jj}=[Cred{rr(i)}*(xred{jj,i}');CimpRed{rr(i),rr(i)}*(xred{jj,i}')];
        end
    end
    time_ROM(jj)=toc;
    Xminus{jj,ns} = xred{jj,i}(end,:)'; 
end
time_ROM_20(1) = time_ROM(20);
Gc_GLE = zeros(nr,ns); Go_GLE = zeros(nr,ns);
for jj=2:nr
    for i = 1:(ns-1)
        if (jj<=numel(Ss{rr(i)}))&&((jj<=numel(Ss{rr(i+1)})))
    
        Gc_GLE(jj-1,i)  = ((Ss{rr(i)}(jj))^2)*norm((Xminus{jj,i}+[Xminus{jj-1,i};0])'*(diag((Ss{rr(i)}(1:jj)).^(-0.5))))^2 ...
                        - ((Ss{rr(i+1)}(jj))^2)*norm((Xplus{jj,i}+[Xplus{jj-1,i};0])'*(diag((Ss{rr(i+1)}(1:jj)).^(-0.5))))^2;
        
        Go_GLE(jj-1,i) =  norm((Xminus{jj,i}-[Xminus{jj-1,i};0])'*(diag((Ss{rr(i)}(1:jj)).^(0.5))))^2 ...
                       -  norm((Xplus{jj,i}-[Xplus{jj-1,i};0])'*(diag((Ss{rr(i+1)}(1:jj)).^(0.5))))^2;

        end
    end
    if (jj<=numel(Ss{rr(i)}))&&((jj<=numel(Ss{rr(i+1)})))
        Go_GLE(jj,ns) =  norm((Xminus{jj,ns}-[Xminus{jj-1,ns};0])'*(diag((Ss{rr(ns)}(1:jj)).^(0.5))))^2;
        Gc_GLE(jj,ns) =  ((Ss{rr(ns)}(jj))^2)*norm((Xminus{jj,ns}+[Xminus{jj-1,ns};0])'*(diag((Ss{rr(ns)}(1:jj)).^(-0.5))))^2;
    end
end
%% 4) Figure 3a
rP=20;
Out_to_Plot =3;
figure
for i=1:ns

    if i==1
        plot(tspan{i},(y{i}(Out_to_Plot,:)'),'--b','LineWidth',LW)
        hold on
        plot(tspan{i},(yred{i,rP}(Out_to_Plot,:)'),'-r','LineWidth',LW)

    end
    plot(tspan{i},(y{i}(Out_to_Plot,:)'),'-b','LineWidth',LW)
    hold on

    plot(tspan{i},(yred{i,rP}(Out_to_Plot,:)'),'--r','LineWidth',LW)

end
xlabel('t','Interpreter','Latex')
ylabel('$\mathbf{y}(t)$','Interpreter','Latex')
lgd=legend(['FOM $n=' num2str(n) '$'],['ROM $r=' num2str(rv(rP)) '$']);
set(lgd, 'Interpreter','Latex','Location','best');

set(gca,'Fontname',FN,'Fontsize',FS);
set(gcf, 'Color', 'w');

% cleanfigure
% matlab2tikz('CMSS_Plot_3a.tex','showInfo', false);
%% Evaluate the norm of X_o, X_c for Standard and Modified Gramians
norm_Xo=zeros(nr,1); norm_Xo_r=zeros(nr,1);
for jj=2:nr
    for j=1:ns

        if (jj<=numel(Ss{rr(j)}))
            norm_Xo(jj)=norm_Xo(jj)+MAX_EIG_o(rr(j),jj)*(((norm(xred{jj,j}(2:(end-1),:)-[xred{jj-1,j}(2:(end-1),:), zeros(numel(xred{jj-1,j}(2:(end-1),1)),1)],"fro").^2))...
                +0.5*(norm(xred{jj,j}(1,:)-[xred{jj-1,j}(1,:),0]).^2)...
                +0.5*(norm(xred{jj,j}(end,:)-[xred{jj-1,j}(end,:), 0]).^2));

             norm_Xo_r(jj)=norm_Xo_r(jj)+MAX_EIG_c(rr(j),jj)*(((norm(xred{jj,j}(2:(end-1),:)+[xred{jj-1,j}(2:(end-1),:), zeros(numel(xred{jj-1,j}(2:(end-1),1)),1)],"fro").^2))...
                +0.5*(norm(xred{jj,j}(1,:)+[xred{jj-1,j}(1,:),0]).^2)...
                +0.5*(norm(xred{jj,j}(end,:)+[xred{jj-1,j}(end,:), 0]).^2));

        end
    end

    if norm_Xo(jj)>0
        norm_Xo(jj)=sqrt(dt)*sqrt(norm_Xo(jj));
    else
        norm_Xo(jj) = 0;
    end
    if norm_Xo_r(jj)>0
        norm_Xo_r(jj)=sqrt(dt)*sqrt(norm_Xo_r(jj));
    else
        norm_Xo_r(jj) = 0;
    end

end
%% Computation of the jumps error estimate component
Gc_GLE = sum(Gc_GLE,2); Go_GLE = sum(Go_GLE,2); 
Jump_Err = Gc_GLE + Go_GLE;
Jump_Err(Jump_Err>0) = 0;
Jump_Err = cumsum((-Jump_Err).^(0.5),'reverse');

%% Computation of the inexact LMI error estimate component
LMI_ERR_COMP = norm_Xo_r+norm_Xo; 
LMI_ERR_COMP(LMI_ERR_COMP<1e-13) = 0; 
LMI_ERR_COMP = cumsum(LMI_ERR_COMP,'reverse');
%% Computation of the SV error component
Dim = zeros(nm,1); Ss_pad = cell(nm,1); 
for i=1:nm
  Dim(i) = numel(Ss{i}); 
end
Dmax = max([max(Dim),rv(end)]); SS = zeros(Dmax,nm); 
for i=1:nm
   Ss_pad{i} = [Ss{i}; zeros( Dmax - length(Ss{i}),1)];
   SS(:,i) = Ss_pad{i};
end

SS = max(SS, [], 2); 
ind = find(SS<tol+0.01*tol); 
SS(ind) = 0; 
% Smallest singular value set to LMI_Phat_Z(2,1) rather than 0, normally 2*tol*SS(1) is preferred but some numerical issue here with the LMI
SS(ind(1)) = LMI_Phat_Z(2,1);
SS = 2*cumsum(SS,'reverse'); 
%% Computing L_2 errors using trapezoidal quadrature rule
final_error_L2_r=zeros(nr,4); 
normU=0; 
for jj=1:nr
final_error_L2=0; normy=0; final_error_L2_m=0;

    for j=1:ns
        final_error_L2=final_error_L2+((norm(y{j}(:,2:(end-1))-yred{j,jj}(:,2:(end-1)),'fro').^2))...
            +0.5*(norm(y{j}(:,1)-yred{j,jj}(:,1)).^2)...
            +0.5*(norm(y{j}(:,end)-yred{j,jj}(:,end)).^2);


        normy=normy+(sum(norm(y{j}(:,2:(end-1)),'fro').^2))...
            +0.5*(norm(y{j}(:,1)).^2)...
            +0.5*(norm(y{j}(:,end)).^2);

        if jj==1
            normU=normU+sum(norm(u(dtt{j}(2:(end-1))'),'fro')^2)+0.5*norm(u(dtt{j}(1)'))^2+...
                0.5*norm(u(dtt{j}(end)'))^2;

            normU=sqrt(dt*normU+imp_input);
        end


    end
    final_error_L2=sqrt(dt)*sqrt(final_error_L2);
    normy = sqrt(dt)*sqrt(normy);

    final_error_L2_r(jj,1)=final_error_L2/normU;
end
%% 5) Data Figure 4a
Err_Estimate(:,1) = SS(rv)+(Jump_Err+LMI_ERR_COMP)/normU;
%% 3) Evaluate ROM without the impulsive part
yred=cell(ns,nr); Xplus = cell(nr,ns); Xminus = cell(nr,ns);
xred = cell(nr,ns); MAX_EIG_c = zeros(nm,nr); MAX_EIG_o = zeros(nm,nr);
time_ROM = zeros(nr,1);
for jj=1:nr
    r=rv(jj);
    % Reduced Operators
    Jred=cell(nm,1); Bred=cell(nm,1); Cred=cell(nm,1); Ered=cell(nm,1); BimpRed = cell(nm,nm); CimpRed = cell(nm,nm);
    U=cell(nm,1); H=cell(nm,1); Ss=cell(nm,1); Vs=cell(nm,1); VV=cell(nm,1); W=cell(nm,1);
    r=r*ones(nm,1);
    for i=1:nm
        
        H{i}=Z2_m{i}'*S2_m{i}; [U{i},Ss{i},Vs{i}]=svd(H{i}); Ss{i}=diag(Ss{i})+tol;
        % Projection Matrices VV and W
        if (r(i)>size(Vs{i},2))||(r(i)>size(U{i},2))
           r(i)=min(size(Vs{i},2),size(U{i},2));
        end
        VV{i}=Z2_m{i}*U{i}(:,1:r(i))*diag((Ss{i}(1:r(i))).^(-0.5));
        W{i}=S2_m{i}*(Vs{i}(:,1:r(i)))*diag((Ss{i}(1:r(i))).^(-0.5));
        Jred{i}=(W{i}'*((J{i}*VV{i})));
        G = eig(Jred{i}); G=max(G); 
        if G>0
            disp(G)
        end
        Ered{i}=(W{i}'*((EE{i}*VV{i})));
        Bred{i}=(W{i}'*Bdiff{i});
        Cred{i}=Cdiff{i}*VV{i};
        for j = 1:nm
            BimpRed{i,j} = W{i}'*Bimp{i,j};
            CimpRed{j,i} = Cimp{j,i}*VV{i};
        end

        % Store conditions
        Lyap_c = Jred{i}*diag((Ss{i}(1:r(i))))+diag((Ss{i}(1:r(i))))*(Jred{i}') + Bred{i}*Bred{i}';
        Lyap_o = (Jred{i}')*diag((Ss{i}(1:r(i))))+diag((Ss{i}(1:r(i))))*Jred{i} + Cred{i}'*Cred{i};
       
        MAX_EIG_c(i,jj) = max(eig(Lyap_c)); MAX_EIG_o(i,jj) = max(eig(Lyap_o));
        
     
    end
    % Reduced solution size r
    tic
    for i=1:ns

        if i==1
            x_0=zeros(r(rr(i)),1);
        else
            if rr(i)==rr(i-1)
                x_0=(W{rr(i)}'*Map{rr(i),rr(i-1)}*VV{rr(i-1)})*(xred{jj,i-1}(end,:)');
            else
                x_0=(W{rr(i)}'*Map{rr(i),rr(i-1)}*VV{rr(i-1)})*(xred{jj,i-1}(end,:)')+BimpRed{rr(i),rr(i-1)}*u_ex{rr(i-1)}(t(i-1));
            end
            Xminus{jj,i-1} = xred{jj,i-1}(end,:)';
            Xplus{jj,i-1} = x_0;
        end
        odefun=@(t,y) Jred{rr(i)}*y+Bred{rr(i)}*u(t);
        options = odeset('RelTol',max([tol*1e-2,1e-12]),'AbsTol',max([tol*1e-2,1e-12]),'Mass',Ered{rr(i)});
        [~,xred{jj,i}]=ode15s(odefun,dtt{i},x_0,options);
        if i<ns
            yred{i,jj}=[Cred{rr(i)}*(xred{jj,i}');CimpRed{rr(i+1),rr(i)}*(xred{jj,i}')];
        else
            yred{i,jj}=[Cred{rr(i)}*(xred{jj,i}');CimpRed{rr(i),rr(i)}*(xred{jj,i}')];
        end
    end
    time_ROM(jj)=toc;
    Xminus{jj,ns} = xred{jj,i}(end,:)'; 
end
time_ROM_20(1) = time_ROM(20);
Gc_GLE = zeros(nr,ns); Go_GLE = zeros(nr,ns);
for jj=2:nr
    for i = 1:(ns-1)
        if (jj<=numel(Ss{rr(i)}))&&((jj<=numel(Ss{rr(i+1)})))
    
        Gc_GLE(jj-1,i)  = ((Ss{rr(i)}(jj))^2)*norm((Xminus{jj,i}+[Xminus{jj-1,i};0])'*(diag((Ss{rr(i)}(1:jj)).^(-0.5))))^2 ...
                        - ((Ss{rr(i+1)}(jj))^2)*norm((Xplus{jj,i}+[Xplus{jj-1,i};0])'*(diag((Ss{rr(i+1)}(1:jj)).^(-0.5))))^2;
        
        Go_GLE(jj-1,i) =  norm((Xminus{jj,i}-[Xminus{jj-1,i};0])'*(diag((Ss{rr(i)}(1:jj)).^(0.5))))^2 ...
                       -  norm((Xplus{jj,i}-[Xplus{jj-1,i};0])'*(diag((Ss{rr(i+1)}(1:jj)).^(0.5))))^2;

        end
    end
    if (jj<=numel(Ss{rr(i)}))&&((jj<=numel(Ss{rr(i+1)})))
        Go_GLE(jj,ns) =  norm((Xminus{jj,ns}-[Xminus{jj-1,ns};0])'*(diag((Ss{rr(ns)}(1:jj)).^(0.5))))^2;
        Gc_GLE(jj,ns) =  ((Ss{rr(ns)}(jj))^2)*norm((Xminus{jj,ns}+[Xminus{jj-1,ns};0])'*(diag((Ss{rr(ns)}(1:jj)).^(-0.5))))^2;
    end
end
%% Evaluate the norm of X_o, X_c for standard and modified GLEs approximate solution
norm_Xo=zeros(nr,1); norm_Xo_r=zeros(nr,1);
for jj=2:nr
    for j=1:ns

        if (jj<=numel(Ss{rr(j)}))
            norm_Xo(jj)=norm_Xo(jj)+MAX_EIG_o(rr(j),jj)*(((norm(xred{jj,j}(2:(end-1),:)-[xred{jj-1,j}(2:(end-1),:), zeros(numel(xred{jj-1,j}(2:(end-1),1)),1)],"fro").^2))...
                +0.5*(norm(xred{jj,j}(1,:)-[xred{jj-1,j}(1,:),0]).^2)...
                +0.5*(norm(xred{jj,j}(end,:)-[xred{jj-1,j}(end,:), 0]).^2));

             norm_Xo_r(jj)=norm_Xo_r(jj)+MAX_EIG_c(rr(j),jj)*(((norm(xred{jj,j}(2:(end-1),:)+[xred{jj-1,j}(2:(end-1),:), zeros(numel(xred{jj-1,j}(2:(end-1),1)),1)],"fro").^2))...
                +0.5*(norm(xred{jj,j}(1,:)+[xred{jj-1,j}(1,:),0]).^2)...
                +0.5*(norm(xred{jj,j}(end,:)+[xred{jj-1,j}(end,:), 0]).^2));

        end
    end

    if norm_Xo(jj)>0
        norm_Xo(jj)=sqrt(dt)*sqrt(norm_Xo(jj));
    else
        norm_Xo(jj) = 0;
    end
    if norm_Xo_r(jj)>0
        norm_Xo_r(jj)=sqrt(dt)*sqrt(norm_Xo_r(jj));
    else
        norm_Xo_r(jj) = 0;
    end

end
%% Computation of the jumps error estimate component
Gc_GLE = sum(Gc_GLE,2); Go_GLE = sum(Go_GLE,2); 
Jump_Err = Gc_GLE + Go_GLE;
Jump_Err(Jump_Err>0) = 0;
Jump_Err = cumsum((-Jump_Err).^(0.5),'reverse');
%% Computation of the LMI error estimate component
LMI_ERR_COMP = norm_Xo_r+norm_Xo; 
LMI_ERR_COMP(LMI_ERR_COMP<1e-13) = 0; 
LMI_ERR_COMP = cumsum(LMI_ERR_COMP,'reverse');
%% Computation of the SV error component
Dim = zeros(nm,1); Ss_pad = cell(nm,1); 
for i=1:nm
  Dim(i) = numel(Ss{i}); 
end
Dmax = max([max(Dim),rv(end)]); SS = zeros(Dmax,nm); 
for i=1:nm
   Ss_pad{i} = [Ss{i}; zeros( Dmax - length(Ss{i}),1)];
   SS(:,i) = Ss_pad{i};
end

SS = max(SS, [], 2); 
ind = find(SS<tol+0.01*tol); 
SS(ind) = 0; 
% Smallest singular value set to LMI_Phat_Z(2,1) rather than 0, normally 2*tol*SS(1) is preferred but some numerical issue here with the LMI
SS(ind(1)) = LMI_Phat_Z(2,1);
SS = 2*cumsum(SS,'reverse'); 
%% Computing L_2 errors using trapezoidal quadrature rule
normU=0; 
for jj=1:nr
final_error_L2=0; normy=0; final_error_L2_m=0;

    for j=1:ns
        final_error_L2=final_error_L2+((norm(y{j}(:,2:(end-1))-yred{j,jj}(:,2:(end-1)),'fro').^2))...
            +0.5*(norm(y{j}(:,1)-yred{j,jj}(:,1)).^2)...
            +0.5*(norm(y{j}(:,end)-yred{j,jj}(:,end)).^2);


        normy=normy+(sum(norm(y{j}(:,2:(end-1)),'fro').^2))...
            +0.5*(norm(y{j}(:,1)).^2)...
            +0.5*(norm(y{j}(:,end)).^2);

        if jj==1
            normU=normU+sum(norm(u(dtt{j}(2:(end-1))'),'fro')^2)+0.5*norm(u(dtt{j}(1)'))^2+...
                0.5*norm(u(dtt{j}(end)'))^2;

            normU=sqrt(dt*normU+imp_input);
        end


    end
    final_error_L2=sqrt(dt)*sqrt(final_error_L2);
    normy = sqrt(dt)*sqrt(normy);

    final_error_L2_r(jj,2)=final_error_L2/normU;
end
Err_Estimate(:,2) = SS(rv)+(Jump_Err+LMI_ERR_COMP)/normU;
%% Numerics for the second switching signal and input function
ns=20;
ns=ns+1;
seed = 123; rng(seed); % For the first switching path
% Choose random order of switches
rr = 1:ns; rr=mod(rr,nm+1); rr(rr==0)=nm;
random_indices = randperm(length(rr));
rr = rr(random_indices); 
% Time scale and time steps
scale=0.2; dt=scale*1e-1;
% Second input function
u=@(t) BStokesNorm*[sin(2*pi*exp(t/8))];
u_ex = cell(nm,1);
u_ex{1} = @(t) BStokesNorm*[sin(2*pi*exp(t/8));(pi*exp(t/8).*cos(2*pi*exp(t/8)))/4;...
    (pi*exp(t/8).*cos(2*pi*exp(t/8)))/32 - ...
    ((pi^2)*exp(t/4).*sin(2*pi*exp(t/8)))/16];
u_ex{2} = @(t) BStokesNorm*[sin(2*pi*exp(t/8));(pi*exp(t/8).*cos(2*pi*exp(t/8)))/4];
% Mass Matrix
for i=1:nm
    EE{i} = speye(nf(i),nf(i));
end
%% Evaluate the FOM
dtt=cell(ns,1); t=zeros(ns,1); tspan=cell(ns,1); y=cell(ns,1); x=cell(ns,1); imp_input = 0;
for i=1:ns
    if i==1
        x_0=zeros(nf(rr(i)),1);
        t(i)=scale*randi([1 10],1,1);
        dtt{i}=0:dt:t(i);
    else
        t(i)=t(i-1)+scale*randi([1 10],1,1);
        dtt{i}=t(i-1):dt:t(i);
        if rr(i)==rr(i-1)
            x_0=Map{rr(i),rr(i-1)}*x{i-1}(end,:)';
        else
            imp_input = imp_input + norm(u_ex{rr(i-1)}(t(i-1)))^2;
            x_0=Map{rr(i),rr(i-1)}*x{i-1}(end,:)'+Bimp{rr(i),rr(i-1)}*u_ex{rr(i-1)}(t(i-1));
        end
    end
    odefun=@(t,y) (J{rr(i)}*y)+Bdiff{rr(i)}*u(t);

    options = odeset('RelTol',max([tol*1e-2,1e-12]),'AbsTol',max([tol*1e-2,1e-12]),'Mass',EE{rr(i)});
    [tspan{i},x{i}]=ode15s(odefun,dtt{i},x_0,options);
    if i<ns
        y{i}=[Cdiff{rr(i)}*(x{i}');Cimp{rr(i+1),rr(i)}*(x{i}')];
    else
        y{i}=[Cdiff{rr(i)}*(x{i}');Cimp{rr(i),rr(i)}*(x{i}')];
    end
end
%% Evaluate the ROM
yred=cell(ns,nr); Xplus = cell(nr,ns); Xminus = cell(nr,ns);
xred = cell(nr,ns); MAX_EIG_c = zeros(nm,nr); MAX_EIG_o = zeros(nm,nr);
time_ROM = zeros(nr,1);
for jj=1:nr
    r=rv(jj);
    % Reduced operators
    Jred=cell(nm,1); Bred=cell(nm,1); Cred=cell(nm,1); Ered=cell(nm,1);
    U=cell(nm,1); H=cell(nm,1); Ss=cell(nm,1); Vs=cell(nm,1); VV=cell(nm,1); W=cell(nm,1);
    r=r*ones(nm,1);
    for i=1:nm
        
        H{i}=Z_m{i}'*S_m{i}; [U{i},Ss{i},Vs{i}]=svd(H{i}); Ss{i}=diag(Ss{i})+tol;
        % Projection Matrices VV and W
        if (r(i)>size(Vs{i},2))||(r(i)>size(U{i},2))
           r(i)=min(size(Vs{i},2),size(U{i},2));
        end
        VV{i}=Z_m{i}*U{i}(:,1:r(i))*diag((Ss{i}(1:r(i))).^(-0.5));
        W{i}=S_m{i}*(Vs{i}(:,1:r(i)))*diag((Ss{i}(1:r(i))).^(-0.5)); 
        Jred{i}=(W{i}'*((J{i}*VV{i})));
        G = eig(Jred{i}); G=max(G); 
        if G>0
            disp(G)
        end
        Ered{i}=(W{i}'*((EE{i}*VV{i})));
        Bred{i}=(W{i}'*Bdiff{i});
        Cred{i}=Cdiff{i}*VV{i};
        for j = 1:nm
            BimpRed{i,j} = W{i}'*Bimp{i,j};
            CimpRed{j,i} = Cimp{j,i}*VV{i};
        end

        % Store conditions
        Lyap_c = Jred{i}*diag((Ss{i}(1:r(i))))+diag((Ss{i}(1:r(i))))*(Jred{i}') + Bred{i}*Bred{i}';
        Lyap_o = (Jred{i}')*diag((Ss{i}(1:r(i))))+diag((Ss{i}(1:r(i))))*Jred{i} + Cred{i}'*Cred{i};
       
        MAX_EIG_c(i,jj) = max(eig(Lyap_c)); MAX_EIG_o(i,jj) = max(eig(Lyap_o));
        
     
    end
    % Reduced solution size r
    tic
    for i=1:ns

        if i==1
            x_0=zeros(r(rr(i)),1);
        else
            if rr(i)==rr(i-1)
                x_0=(W{rr(i)}'*Map{rr(i),rr(i-1)}*VV{rr(i-1)})*(xred{jj,i-1}(end,:)');
            else
                x_0=(W{rr(i)}'*Map{rr(i),rr(i-1)}*VV{rr(i-1)})*(xred{jj,i-1}(end,:)')+BimpRed{rr(i),rr(i-1)}*u_ex{rr(i-1)}(t(i-1));
            end
            Xminus{jj,i-1} = xred{jj,i-1}(end,:)';
            Xplus{jj,i-1} = x_0;
        end
        odefun=@(t,y) Jred{rr(i)}*y+Bred{rr(i)}*u(t);
        options = odeset('RelTol',max([tol*1e-2,1e-12]),'AbsTol',max([tol*1e-2,1e-12]),'Mass',Ered{rr(i)});
        [~,xred{jj,i}]=ode15s(odefun,dtt{i},x_0,options);
        if i<ns
            yred{i,jj}=[Cred{rr(i)}*(xred{jj,i}');CimpRed{rr(i+1),rr(i)}*(xred{jj,i}')];
        else
            yred{i,jj}=[Cred{rr(i)}*(xred{jj,i}');CimpRed{rr(i),rr(i)}*(xred{jj,i}')];
        end
    end
    time_ROM(jj)=toc;
    Xminus{jj,ns} = xred{jj,i}(end,:)'; 
end
time_ROM(2)=time_ROM(20);
Gc_GLE = zeros(nr,ns); Go_GLE = zeros(nr,ns);
for jj=2:nr
    for i = 1:(ns-1)
        if (jj<=numel(Ss{rr(i)}))&&((jj<=numel(Ss{rr(i+1)})))
    
        Gc_GLE(jj-1,i)  = ((Ss{rr(i)}(jj))^2)*norm((Xminus{jj,i}+[Xminus{jj-1,i};0])'*(diag((Ss{rr(i)}(1:jj)).^(-0.5))))^2 ...
                        - ((Ss{rr(i+1)}(jj))^2)*norm((Xplus{jj,i}+[Xplus{jj-1,i};0])'*(diag((Ss{rr(i+1)}(1:jj)).^(-0.5))))^2;
        
        Go_GLE(jj-1,i) =  norm((Xminus{jj,i}-[Xminus{jj-1,i};0])'*(diag((Ss{rr(i)}(1:jj)).^(0.5))))^2 ...
                       -  norm((Xplus{jj,i}-[Xplus{jj-1,i};0])'*(diag((Ss{rr(i+1)}(1:jj)).^(0.5))))^2;

        end
    end
    if (jj<=numel(Ss{rr(i)}))&&((jj<=numel(Ss{rr(i+1)})))
        Go_GLE(jj,ns) =  norm((Xminus{jj,ns}-[Xminus{jj-1,ns};0])'*(diag((Ss{rr(ns)}(1:jj)).^(0.5))))^2;
        Gc_GLE(jj,ns) =  ((Ss{rr(ns)}(jj))^2)*norm((Xminus{jj,ns}+[Xminus{jj-1,ns};0])'*(diag((Ss{rr(ns)}(1:jj)).^(-0.5))))^2;
    end
end
%% 4) Figure 3b
rP=20;
Out_to_Plot =1;
figure
for i=1:ns

    if i==1
        plot(tspan{i},(y{i}(Out_to_Plot,:)'),'-b','LineWidth',LW)
        hold on
        plot(tspan{i},(yred{i,rP}(Out_to_Plot,:)'),'--r','LineWidth',LW)

    end
    plot(tspan{i},(y{i}(Out_to_Plot,:)'),'-b','LineWidth',LW)
    hold on

    plot(tspan{i},(yred{i,rP}(Out_to_Plot,:)'),'--r','LineWidth',LW)

end

xlabel('t','Interpreter','Latex')
ylabel('$\mathbf{y}(t)$','Interpreter','Latex')
lgd=legend(['FOM $n=' num2str(n) '$'],['ROM $r=' num2str(rv(rP)) '$']);
set(lgd, 'Interpreter','Latex','Location','best');

set(gca,'Fontname',FN,'Fontsize',FS);
set(gcf, 'Color', 'w');
% cleanfigure
% matlab2tikz('CMSS_Plot_3b.tex','showInfo', false);
%% Evaluate the norm of X_o, X_c for standard and modified GLEs approximate solutions
norm_Xo_m=zeros(nr,1); norm_Xo=zeros(nr,1); norm_Xo_r=zeros(nr,1); norm_Xo_m_r=zeros(nr,1); 
for jj=2:nr
    for j=1:ns

        if (jj<=numel(Ss{rr(j)}))
            norm_Xo(jj)=norm_Xo(jj)+MAX_EIG_o(rr(j),jj)*(((norm(xred{jj,j}(2:(end-1),:)-[xred{jj-1,j}(2:(end-1),:), zeros(numel(xred{jj-1,j}(2:(end-1),1)),1)],"fro").^2))...
                +0.5*(norm(xred{jj,j}(1,:)-[xred{jj-1,j}(1,:),0]).^2)...
                +0.5*(norm(xred{jj,j}(end,:)-[xred{jj-1,j}(end,:), 0]).^2));

             norm_Xo_r(jj)=norm_Xo_r(jj)+MAX_EIG_c(rr(j),jj)*(((norm(xred{jj,j}(2:(end-1),:)+[xred{jj-1,j}(2:(end-1),:), zeros(numel(xred{jj-1,j}(2:(end-1),1)),1)],"fro").^2))...
                +0.5*(norm(xred{jj,j}(1,:)+[xred{jj-1,j}(1,:),0]).^2)...
                +0.5*(norm(xred{jj,j}(end,:)+[xred{jj-1,j}(end,:), 0]).^2));

        end
    end

    if norm_Xo(jj)>0
        norm_Xo(jj)=sqrt(dt)*sqrt(norm_Xo(jj));
    else
        norm_Xo(jj) = 0;
    end
    if norm_Xo_m(jj)>0
        norm_Xo_m(jj)=sqrt(dt)*sqrt(norm_Xo_m(jj));
    else
        norm_Xo_m(jj) = 0;
    end
    if norm_Xo_r(jj)>0
        norm_Xo_r(jj)=sqrt(dt)*sqrt(norm_Xo_r(jj));
    else
        norm_Xo_r(jj) = 0;
    end

end
%% Computation of the jumps error estimate component
Gc_GLE = sum(Gc_GLE,2); Go_GLE = sum(Go_GLE,2); 
Jump_Err = Gc_GLE + Go_GLE;
Jump_Err(Jump_Err>0) = 0;
Jump_Err = cumsum((-Jump_Err).^(0.5),'reverse');
%% Computation of the LMI  error estimate component
LMI_ERR_COMP = norm_Xo_r+norm_Xo; 
LMI_ERR_COMP(LMI_ERR_COMP<1e-12) = 0; 
LMI_ERR_COMP = cumsum(LMI_ERR_COMP,'reverse');
%% Computation of the SV error component
Dim = zeros(nm,1); Ss_pad = cell(nm,1); Ss_pad_m = cell(nm,1); 
for i=1:nm
  Dim(i) = numel(Ss{i}); 
end
Dmax = max([max(Dim),rv(end)]); SS = zeros(Dmax,nm); 
for i=1:nm
   Ss_pad{i} = [Ss{i}; zeros( Dmax - length(Ss{i}),1)];
   SS(:,i) = Ss_pad{i};
end

SS = max(SS, [], 2); 
ind = find(SS<tol+0.01*tol); 
SS(ind) = 0; 
SS(ind(1)) = LMI_Phat_Z(2,1);
SS = 2*cumsum(SS,'reverse'); 
%% Computing L_2 errors using trapezoidal quadrature rule
normU=0; 
for jj=1:nr
final_error_L2=0; normy=0; final_error_L2_m=0;

    for j=1:ns
        final_error_L2=final_error_L2+((norm(y{j}(:,2:(end-1))-yred{j,jj}(:,2:(end-1)),'fro').^2))...
            +0.5*(norm(y{j}(:,1)-yred{j,jj}(:,1)).^2)...
            +0.5*(norm(y{j}(:,end)-yred{j,jj}(:,end)).^2);


        normy=normy+(sum(norm(y{j}(:,2:(end-1)),'fro').^2))...
            +0.5*(norm(y{j}(:,1)).^2)...
            +0.5*(norm(y{j}(:,end)).^2);

        if jj==1
            normU=normU+sum(norm(u(dtt{j}(2:(end-1))'),'fro')^2)+0.5*norm(u(dtt{j}(1)'))^2+...
                0.5*norm(u(dtt{j}(end)'))^2;

            normU=sqrt(dt)*sqrt(normU+imp_input);
        end


    end
    final_error_L2=sqrt(dt)*sqrt(final_error_L2);
    normy = sqrt(dt)*sqrt(normy);

    final_error_L2_r(jj,3)=final_error_L2/(normU);
end
%% Data for figure 4a 
Err_Estimate(:,3) = SS(rv)+(Jump_Err+LMI_ERR_COMP)/normU;
%% Evaluate ROM without impulsive components
yred=cell(ns,nr); Xplus = cell(nr,ns); Xminus = cell(nr,ns);
xred = cell(nr,ns); MAX_EIG_c = zeros(nm,nr); MAX_EIG_o = zeros(nm,nr);
time_ROM = zeros(nr,1);
for jj=1:nr
    r=rv(jj);
    % Reduced operators
    Jred=cell(nm,1); Bred=cell(nm,1); Cred=cell(nm,1); Ered=cell(nm,1);
    U=cell(nm,1); H=cell(nm,1); Ss=cell(nm,1); Vs=cell(nm,1); VV=cell(nm,1); W=cell(nm,1);
    r=r*ones(nm,1);
    for i=1:nm
        
        H{i}=Z2_m{i}'*S2_m{i}; [U{i},Ss{i},Vs{i}]=svd(H{i}); Ss{i}=diag(Ss{i})+tol;
        % Projection Matrices VV and W
        if (r(i)>size(Vs{i},2))||(r(i)>size(U{i},2))
           r(i)=min(size(Vs{i},2),size(U{i},2));
        end
        VV{i}=Z2_m{i}*U{i}(:,1:r(i))*diag((Ss{i}(1:r(i))).^(-0.5));
        W{i}=S2_m{i}*(Vs{i}(:,1:r(i)))*diag((Ss{i}(1:r(i))).^(-0.5)); 
        Jred{i}=(W{i}'*((J{i}*VV{i})));
        G = eig(Jred{i}); G=max(G); 
        if G>0
            disp(G)
        end
        Ered{i}=(W{i}'*((EE{i}*VV{i})));
        Bred{i}=(W{i}'*Bdiff{i});
        Cred{i}=Cdiff{i}*VV{i};
        for j = 1:nm
            BimpRed{i,j} = W{i}'*Bimp{i,j};
            CimpRed{j,i} = Cimp{j,i}*VV{i};
        end

        % Store conditions
        Lyap_c = Jred{i}*diag((Ss{i}(1:r(i))))+diag((Ss{i}(1:r(i))))*(Jred{i}') + Bred{i}*Bred{i}';
        Lyap_o = (Jred{i}')*diag((Ss{i}(1:r(i))))+diag((Ss{i}(1:r(i))))*Jred{i} + Cred{i}'*Cred{i};
       
        MAX_EIG_c(i,jj) = max(eig(Lyap_c)); MAX_EIG_o(i,jj) = max(eig(Lyap_o));
    end
    % Reduced solution size r
    tic
    for i=1:ns

        if i==1
            x_0=zeros(r(rr(i)),1);
        else
            if rr(i)==rr(i-1)
                x_0=(W{rr(i)}'*Map{rr(i),rr(i-1)}*VV{rr(i-1)})*(xred{jj,i-1}(end,:)');
            else
                x_0=(W{rr(i)}'*Map{rr(i),rr(i-1)}*VV{rr(i-1)})*(xred{jj,i-1}(end,:)')+BimpRed{rr(i),rr(i-1)}*u_ex{rr(i-1)}(t(i-1));
            end
            Xminus{jj,i-1} = xred{jj,i-1}(end,:)';
            Xplus{jj,i-1} = x_0;
        end
        odefun=@(t,y) Jred{rr(i)}*y+Bred{rr(i)}*u(t);
        options = odeset('RelTol',max([tol*1e-2,1e-12]),'AbsTol',max([tol*1e-2,1e-12]),'Mass',Ered{rr(i)});
        [~,xred{jj,i}]=ode15s(odefun,dtt{i},x_0,options);
        if i<ns
            yred{i,jj}=[Cred{rr(i)}*(xred{jj,i}');CimpRed{rr(i+1),rr(i)}*(xred{jj,i}')];
        else
            yred{i,jj}=[Cred{rr(i)}*(xred{jj,i}');CimpRed{rr(i),rr(i)}*(xred{jj,i}')];
        end
    end
    time_ROM(jj)=toc;
    Xminus{jj,ns} = xred{jj,i}(end,:)'; 
end
time_ROM(2)=time_ROM(20);
Gc_GLE = zeros(nr,ns); Go_GLE = zeros(nr,ns);
for jj=2:nr
    for i = 1:(ns-1)
        if (jj<=numel(Ss{rr(i)}))&&((jj<=numel(Ss{rr(i+1)})))
    
        Gc_GLE(jj-1,i)  = ((Ss{rr(i)}(jj))^2)*norm((Xminus{jj,i}+[Xminus{jj-1,i};0])'*(diag((Ss{rr(i)}(1:jj)).^(-0.5))))^2 ...
                        - ((Ss{rr(i+1)}(jj))^2)*norm((Xplus{jj,i}+[Xplus{jj-1,i};0])'*(diag((Ss{rr(i+1)}(1:jj)).^(-0.5))))^2;
        
        Go_GLE(jj-1,i) =  norm((Xminus{jj,i}-[Xminus{jj-1,i};0])'*(diag((Ss{rr(i)}(1:jj)).^(0.5))))^2 ...
                       -  norm((Xplus{jj,i}-[Xplus{jj-1,i};0])'*(diag((Ss{rr(i+1)}(1:jj)).^(0.5))))^2;

        end
    end
    if (jj<=numel(Ss{rr(i)}))&&((jj<=numel(Ss{rr(i+1)})))
        Go_GLE(jj,ns) =  norm((Xminus{jj,ns}-[Xminus{jj-1,ns};0])'*(diag((Ss{rr(ns)}(1:jj)).^(0.5))))^2;
        Gc_GLE(jj,ns) =  ((Ss{rr(ns)}(jj))^2)*norm((Xminus{jj,ns}+[Xminus{jj-1,ns};0])'*(diag((Ss{rr(ns)}(1:jj)).^(-0.5))))^2;
    end
end
%% Evaluate the norm of X_o, X_c for standard and modified GLE approximate solutions
norm_Xo_m=zeros(nr,1); norm_Xo=zeros(nr,1); norm_Xo_r=zeros(nr,1); norm_Xo_m_r=zeros(nr,1); 
for jj=2:nr
    for j=1:ns

        if (jj<=numel(Ss{rr(j)}))
            norm_Xo(jj)=norm_Xo(jj)+MAX_EIG_o(rr(j),jj)*(((norm(xred{jj,j}(2:(end-1),:)-[xred{jj-1,j}(2:(end-1),:), zeros(numel(xred{jj-1,j}(2:(end-1),1)),1)],"fro").^2))...
                +0.5*(norm(xred{jj,j}(1,:)-[xred{jj-1,j}(1,:),0]).^2)...
                +0.5*(norm(xred{jj,j}(end,:)-[xred{jj-1,j}(end,:), 0]).^2));

             norm_Xo_r(jj)=norm_Xo_r(jj)+MAX_EIG_c(rr(j),jj)*(((norm(xred{jj,j}(2:(end-1),:)+[xred{jj-1,j}(2:(end-1),:), zeros(numel(xred{jj-1,j}(2:(end-1),1)),1)],"fro").^2))...
                +0.5*(norm(xred{jj,j}(1,:)+[xred{jj-1,j}(1,:),0]).^2)...
                +0.5*(norm(xred{jj,j}(end,:)+[xred{jj-1,j}(end,:), 0]).^2));

        end
    end

    if norm_Xo(jj)>0
        norm_Xo(jj)=sqrt(dt)*sqrt(norm_Xo(jj));
    else
        norm_Xo(jj) = 0;
    end
    if norm_Xo_m(jj)>0
        norm_Xo_m(jj)=sqrt(dt)*sqrt(norm_Xo_m(jj));
    else
        norm_Xo_m(jj) = 0;
    end
    if norm_Xo_r(jj)>0
        norm_Xo_r(jj)=sqrt(dt)*sqrt(norm_Xo_r(jj));
    else
        norm_Xo_r(jj) = 0;
    end

end
%% Computation of the jumps error estimate component
Gc_GLE = sum(Gc_GLE,2); Go_GLE = sum(Go_GLE,2); 
Jump_Err = Gc_GLE + Go_GLE;
Jump_Err(Jump_Err>0) = 0;
Jump_Err = cumsum((-Jump_Err).^(0.5),'reverse');
%% Computation of the LMI  error estimate component
LMI_ERR_COMP = norm_Xo_r+norm_Xo; 
LMI_ERR_COMP(LMI_ERR_COMP<1e-12) = 0; 
LMI_ERR_COMP = cumsum(LMI_ERR_COMP,'reverse');
%% Computation of the SV error component
Dim = zeros(nm,1); Ss_pad = cell(nm,1); Ss_pad_m = cell(nm,1); 
for i=1:nm
  Dim(i) = numel(Ss{i}); 
end
Dmax = max([max(Dim),rv(end)]); SS = zeros(Dmax,nm); 
for i=1:nm
   Ss_pad{i} = [Ss{i}; zeros( Dmax - length(Ss{i}),1)];
   SS(:,i) = Ss_pad{i};
end

SS = max(SS, [], 2); 
ind = find(SS<tol+0.01*tol); 
SS(ind) = 0; 
SS(ind(1)) = LMI_Phat_Z(2,1);
SS = 2*cumsum(SS,'reverse'); 
%% Computing L_2 errors using trapezoidal quadrature rule
si=size(B{1},2);
normU=0; 
for jj=1:nr
final_error_L2=0; normy=0; final_error_L2_m=0;

    for j=1:ns
        final_error_L2=final_error_L2+((norm(y{j}(:,2:(end-1))-yred{j,jj}(:,2:(end-1)),'fro').^2))...
            +0.5*(norm(y{j}(:,1)-yred{j,jj}(:,1)).^2)...
            +0.5*(norm(y{j}(:,end)-yred{j,jj}(:,end)).^2);


        normy=normy+(sum(norm(y{j}(:,2:(end-1)),'fro').^2))...
            +0.5*(norm(y{j}(:,1)).^2)...
            +0.5*(norm(y{j}(:,end)).^2);

        if jj==1
            normU=normU+sum(norm(u(dtt{j}(2:(end-1))'),'fro')^2)+0.5*norm(u(dtt{j}(1)'))^2+...
                0.5*norm(u(dtt{j}(end)'))^2;

            normU=sqrt(dt)*sqrt(normU+imp_input);
        end


    end
    final_error_L2=sqrt(dt)*sqrt(final_error_L2);
    normy = sqrt(dt)*sqrt(normy);

    final_error_L2_r(jj,4)=final_error_L2/(normU);
end
%% Saving data
Err_Estimate(:,4) = SS(rv)+(Jump_Err+LMI_ERR_COMP)/normU;
%% 5) Figure 4a
figure
semilogy(rv,final_error_L2_r(:,1),'-ob','LineWidth',LW)
hold on
semilogy(rv,final_error_L2_r(:,3),'--+b','LineWidth',LW)
semilogy(rv,Err_Estimate(:,1),'-or','LineWidth',LW)
semilogy(rv,Err_Estimate(:,3),'--+r','LineWidth',LW)

xlabel('$r$','Interpreter','Latex')
lgd=legend('$\varepsilon(r,\mathbf{u}_1)$','$\varepsilon(r,\mathbf{u}_2)$','$\tau_{rel}(r,\mathbf{u}_1)$','$\tau_{rel}(r,\mathbf{u}_2)$');
set(lgd, 'Interpreter','Latex','Location','best');
set(gca,'Fontname',FN,'Fontsize',1.2*FS);
set(gcf, 'Color', 'w');
%% 5) Figure 4b
figure
semilogy(rv,final_error_L2_r(:,1),'-ob','LineWidth',LW)
hold on
semilogy(rv,final_error_L2_r(:,2),'--+b','LineWidth',LW)
semilogy(rv,final_error_L2_r(:,3),'-or','LineWidth',LW)
semilogy(rv,final_error_L2_r(:,4),'--+r','LineWidth',LW)

xlabel('$r$','Interpreter','Latex')
lgd=legend('$\varepsilon(r,\mathbf{u}_1)$','$\varepsilon(r,\mathbf{u}_1)$ (no imp.)','$\varepsilon(r,\mathbf{u}_1)$','$\varepsilon(r,\mathbf{u}_1)$ (no imp.)');
set(lgd, 'Interpreter','Latex','Location','best');
set(gca,'Fontname',FN,'Fontsize',1.2*FS);
set(gcf, 'Color', 'w');




