%% Code based on the work:
%[1] M. Manucci and B. Unger, Balancing-based model reduction for switched descriptor systems
% ArXiv e-print 2404.10511, 2024.
%% --------------------------------------------------------------
clc
clearvars
close all
addpath(genpath('./Functions/'))
FS = 15;       % Fontsize
FN = 'times';  % Fontname
LW = 2;      % Linewidth
MS = 7.8;    % Markersize
%% Declare User Input
flag=input('Digit 1 for the mass-spring problem, digit 2 for the Stokes problem\n');
%% Script for Large Scale DAEs
if flag==1
    g=500;  %input('Size of the problem g, the dimension will be n=2g+1: \n');
    nm=5;    %input('Number of Modes, minimum 2: \n');
    tol=1e-10; %input('Accuracy required for the ROM: \n');
    n=2*g+1; nf=n-3;
    DAE_Index=3;
    dimKerE=1;
    mas=ones(g,1);
    k1=1.5*ones(g-1,1); k2=2*ones(g,1);
    d1=0.7*ones(g-1,1); d2=0.7*ones(g,1);
    %% Definition of the First System
    [E1, A1, B1, C1, M1, D1, K1, G1] = msd_ind3(g, mas, k1, k2, d1, d2);
    A{1}=A1; E{1}=E1; B{1}=B1; C{1}=C1; I=speye(n);
    si=size(B1,2); so=size(C1,1);
    seed = 123; rng(seed); %fixing the seed for generation of random numbers (in order to make experiments reproducible)
    for i=2:nm

        D2=D1+rand*(0.35)*speye(g);
        AA=sparse(n,n);
        AA(1:g,g+1:2*g)=speye(g);
        AA(g+1:2*g,1:g)=K1;
        AA(g+1:2*g,g+1:2*g)=D2;
        AA(2*g+1:end,1:g)=G1;
        AA(g+1:2*g,2*g+1:end)=-G1';
        A{i}=AA;
        % Change Algebraic contriants
        A{i}(end,i+1)=0.5;
        %% Test probelm without Input-Dependent Jump
        E{i}=E1+rand*E1; E{i}=sparse(E{i}); dimKerE(i)=1;
        B{i}=B1+I(:,i)+I(:,g+i);
        C{i}=C1;
        %% Test Problem for Input Dependent Jump (only two mode)
        if nm>2
            fprintf('To run Input Dependent Jump test problem nm must be equal 2 \n');
        else
            E{i}=E1; E{i}(end-i+1,(end-i+1):end)=0; E{i}=sparse(E{i}); dimKerE(i)=i;
            q=find(A{i}(end-i+1,:)); v1=A{i}(end-i+1,:); v1=circshift(v1,-1); A{i}(end-i,:)=0; A{i}(end-i,:)=v1;
            B{i}=B1; C{i}=C1;
        end
    end
end
%% Stokes
if flag==2
    nx=input('Discretization points in the first dimension\n');
    ny=input('Discretization points in the second dimension\n');
    si=3; %input('Number of input\n');
    so=3; %input('Number of output\n');
    nm=5; %input('Number of Modes, minimum 2: \n');
    tol=1e-10; %input('Accuracy required for the ROM: \n');
    DAE_Index=2;
    n = (nx-1)*ny+(ny-1)*nx+nx*ny-1; I=speye(n);
    dimKerE(1)=nx*ny-1; opts=[];
    [E1, A1, B1, C1, nf] = stokes_ind2(si, so, nx, ny,opts);
    A{1}=A1; E{1}=E1; B{1}=B1; C{1}=C1;

    AA=sparse(n,n);
    AA(1:((nx-1)*ny+(ny-1)*nx),1:((nx-1)*ny+(ny-1)*nx))=A{1}(1:((nx-1)*ny+(ny-1)*nx),1:((nx-1)*ny+(ny-1)*nx));
    mu=linspace(-0.35,0.35,nm);
    seed = 123; rng(seed);
    for i=2:nm
        dimKerE(i)=nx*ny-1;
        A{i}=A1+mu(i)*AA;
        alpha=1;
        E{i}=E1+alpha*rand*E1; E{i}=sparse(E{i});
        B{i}=B1+[I(:,i),I(:,i+1),I(:,i+(nx-1)*ny+(ny-1))];
        C{i}=C1+[I(:,i),I(:,i+1),I(:,i+(nx-1)*ny+(ny-1))]';
    end
end
%% Determing the Quasi-Waiestrass Form for each mode
% NOTE: V and W are formed explicitly for the two test problems here considered, in general
% this should (and could) be avoided.
V=cell(nm,1); W=cell(nm,1);
T=cell(nm,1); S=cell(nm,1);
DS=cell(nm,1); IS=cell(nm,1);
n1=cell(nm,1); time_Wong=zeros(nm,1); I=speye(n); Iimp=cell(nm,1);
for i=1:nm
    tic
    [V{i},W{i}] = Wong_Sequance(A{i},I(:,(n-dimKerE(i)+1):end),DAE_Index);
    time_Wong(i)=toc;
    T{i}=[V{i},W{i}];
    S{i}=([E{i}*V{i},A{i}*W{i}]); S{i}(abs(S{i})<1e-12)=0; S{i}=sparse(S{i});

    n1{i}=size(V{i},2);
    Q=sparse(n,n); Q(1:n1{i},1:n1{i})=speye(n1{i});

    DS{i}=T{i}*Q;     DS{i}(abs(DS{i})<1e-12)=0; DS{i}=sparse(DS{i});
    IS{i}=T{i}*(I-Q); IS{i}(abs(IS{i})<1e-12)=0; IS{i}=sparse(IS{i});
    Iimp{i}=(I-Q);
end
AV_WONG_TIME=mean(time_Wong);

%% Defining Input-Output matrices of Differential-Impulsive components
Bd=cell(nm,1); Bi=cell(nm,1); ni=cell(nm,1);
Cd=cell(nm,1); Ci=cell(nm,1); no=cell(nm,1);
for i=1:nm
    Q=sparse(n,n);              Q(1:n1{i},1:n1{i})=speye(n1{i});
    Bd{i}= DS{i}*(S{i}\B{i});   Cd{i}=(T{i}'\((C{i}*T{i}*Q)'))';
    Bi{i}=IS{i}*(S{i}\B{i});    Ci{i}=C{i}-Cd{i};
    ni{i}=size(B{i},2);         no{i}=size(C{i},1);
    for j=1:(DAE_Index-1)
        Bi{i}(:,(ni{i}*j+1):(ni{i}*(j+1)))=IS{i}*(S{i}\(E{i}*Bi{i}(:,(ni{i}*(j-1)+1):(ni{i}*(j)))));
    end
end
BB=[]; CC=[]; TT=T{1};
I_nJ1_nJj=cell(nm,1); I_nJj=cell(nm,1); BBimp=cell(nm,1);
for i=1:nm
    I_nJ1_nJj{i}=I(1:n1{1},1:n1{i}); I_nJ1_nJj{i}((n1{i}+1):end,:)=0;
    I_nJj{i}=I(1:n1{i},:); I_nJj{i}(:,(n1{i}+1):end)=0;
    BB=[BB,I_nJ1_nJj{i}*I_nJj{i}*(S{i}\B{i})];
    CC=[CC;(T{i}'\(C{i}*T{i}*I_nJj{i}'*I_nJj{i})')'*T{1}*(I_nJj{1}')];

    BBimp{i}=Iimp{i}*(S{i}\B{i}); Bimp_n=BBimp{i};
    for ii=1:(DAE_Index-1)
        Bimp_n=Iimp{i}*(S{i}\((E{i}*T{i})*Iimp{i}*Bimp_n));
        BBimp{i}=[BBimp{i}, Bimp_n];
    end
    BBimp{i}= I_nJj{1}*(T{1}\(T{i}*BBimp{i}));
end
BB_t=BB;
for i=1:nm

    if normest(BBimp{i})>1e-12
        BB_t=[BB_t,BBimp{i}];
    end

end
%% Computing the Gramians
% Approximate Solution of GLE for P
tic
[Z,rho_P]=Solve_LS_GLE(A,S,T,BB,n1,tol,1); % Use BB_t instead of BB is you want to test for input dependent jumps
time_P=toc;
% Approximate Solution of GLE for Q
tic
[S_q,rho_Q]=Solve_LS_GLE(A,S,T,CC,n1,tol,2);
time_Q=toc;

SZ=svd(Z); SS_q=svd(S_q); 

%% Computing the constant CP
deltaZ=zeros(numel(SZ),1);
deltaZ(1)=SZ(1)*(SZ(1)^2-SZ(2)^2);
for i=2:(numel(SZ)-1)
    deltaZ(i)=max(SZ(i)*(SZ(i)^2-SZ(i+1)^2),SZ(i)*(SZ(i-1)^2-SZ(i)^2));
    if deltaZ(i)>1/sqrt(tol)
        deltaZ(i)=1/sqrt(tol);
    end
end
deltaZ(numel(SZ))=max(SZ(end)/(SZ(end)^2),SZ(end)/(SZ(end-1)^2-SZ(end)^2));
if deltaZ(end)>1/sqrt(tol)
        deltaZ(end)=1/sqrt(tol);
end
CP=sum(deltaZ);

%% Computing the constant CQ

deltaS=zeros(numel(SS_q),1);
deltaS(1)=SS_q(1)/(SS_q(1)^2-SS_q(2)^2);
for i=2:(numel(SS_q)-1)
    deltaS(i)=max(SS_q(i)/(SS_q(i)^2-SS_q(i+1)^2),SS_q(i)/(SS_q(i-1)^2-SS_q(i)^2));
    if deltaS(i)>1/sqrt(tol)
        deltaS(i)=1/sqrt(tol);
    end
end
deltaS(numel(SS_q))=max(SS_q(end)/(SS_q(end)^2),SS_q(end)/(SS_q(end-1)^2-SS_q(end)^2));
if deltaS(end)>1/sqrt(tol)
   deltaS(end)=1/sqrt(tol);
end
CQ=sum(deltaS);

%% ---------------------------------

% Write Gramians in the system general coordinates and full dimension
Z_l = T{1}*(I_nJj{1}'*Z);
S_l = (T{1}'\(I_nJj{1}'*S_q));

H=Z_l'*S_l; [U,Sss,Vs]=svd(H);
% Rescale the singluar values
Ss=diag(Sss)*sqrt(rho_P)*sqrt(rho_Q);

NH=size(Ss,1); normQin=norm(pinv(S_q)); normPin=norm(pinv(Z));
Fun_Err=@(i) 2*sum(Ss(i+1:end))+2*(NH-i)*((normPin+CP)*norm(S_q)*tol+(normQin+CQ)*norm(Z)*tol...
              +norm(Z)*SS_q(end)+norm(S_q)*SZ(end));
decay=zeros(NH,1); decay_Sv=zeros(NH,1);
for i=0:(NH-1)
    decay(i+1)=Fun_Err(i);
    decay_Sv(i+1)=2*sum(Ss(i+1:end));
end
%% Plot decay of Sum of Singular Values and Modified Sum of Singular Values
figure
semilogy(decay,'--b','LineWidth',LW)
hold on
semilogy(decay_Sv,'-r','LineWidth',LW)
xlabel('$r$','Interpreter','Latex')

set(gca,'Fontname',FN,'Fontsize',FS);
set(gcf, 'Color', 'w');
%% Control on LMI (only when dimension of the problem is <1000)
if n<1000
    lambda=zeros(nm,1);
    for i=1:nm
        Zj=I_nJj{i}*(T{i}\(T{1}*I_nJj{1}'*Z)); %Write Gramians in terms of the coordinates of the ith mode
        J=I_nJj{i}*(S{i}\(A{i}*T{i}))*(I_nJj{i}'); % Form J_i explicitly
        % Need to rescale the Gramian by rho_P factor
        Zj=sqrt(rho_P)*Zj;
        M=J*(Zj*Zj')+(Zj*Zj')*J'+((I_nJ1_nJj{i}*I_nJj{i}*(S{i}\B{i}))*(I_nJ1_nJj{i}*I_nJj{i}*(S{i}\B{i}))');
        lambda(i)=max(real(eigs(M,n)));
        if lambda(i)>0
            fprintf('LMI is not satisfied for mode %d, largest eigenvalue is %d \n',i,lambda(i))
        end
    end
end
%% Function to compute diffential and impulsive subspaces
function [V,W] = Wong_Sequance(A, K, DAE_Index)
%% Space W
n=size(A,1); I=speye(n);
Kn_old=K; W=K; normA=sqrt(A(:)'*A(:));
for i=1:DAE_Index-1
    Kn=A*Kn_old/normA; Kn_oldTKn_old=W'*W; Kn_old_new=W;
    for j=1:size(Kn,2)
        alpha_W=(Kn_oldTKn_old)\(Kn_old_new'*Kn(:,j));
        normRes_W(j)=alpha_W'*Kn_old_new'*Kn_old_new*alpha_W-2*Kn(:,j)'*Kn_old_new*alpha_W+Kn(:,j)'*Kn(:,j);
        if abs(normRes_W(j))>1e-10
            Kn_oldTKn_old=[Kn_oldTKn_old,Kn_old_new'*Kn(:,j);Kn(:,j)'*Kn_old_new,Kn(:,j)'*Kn(:,j)];
            Kn_old_new=[Kn_old_new,Kn(:,j)];
        end
    end
    zero_columns = find(abs(normRes_W)<1e-11);
    Kn_old(:,zero_columns)=[];
    Kn_old=A*Kn_old/normA;  W=[W,Kn_old/normest(Kn_old)];
end
%% Determine V as the complementar space to W
V=[];
Kort=[I(:,(1:(n-size(K,2))))];
Kort2=Kort;
KK=[Kort,K];
WW=W(:,(size(K,2)+1):end);
alpha=KK'*WW; g=[]; flag=0;
while flag==0
    res=-alpha(:,1)'*alpha(:,1)+WW(:,1)'*WW(:,1);
    if abs(res)<1e-14
        q=find(alpha(:,1));  j=1;
        v=KK(:,q(1)); Kort(:,q(1))=[];
        KK(:,q(1))=[];
        if (size(Kort,2)+size(W,2))==n
            break;
        end
        % Orthgonalize v with respect to the sparse matrix
        for jj=1:3
            v=v-KK*(KK'*v);
            v=v/normest(v);
        end
        KK=[KK,v]; WW=WW(:,2:end);
        alpha=KK'*WW;
    end
end
V=Kort;
end