%% Code based on the work:
%[1] M. Manucci and B. Unger, Solving Generalized Lyapunov Equations with guarantees:
% application to the Reduction of Linear Switched Systems
% ArXiv e-print ..., 2024.
%% --------------------------------------------------------------
function [S,rho] = Solve_LSS_GLE_2(A,E,B,opts_GLE,flag)
if isfield(opts_GLE,'rescale')
    rescale = opts_GLE.rescale;
    gamma= 1;
else
    rescale = 1;
end
if isfield(opts_GLE,'RelativeErr')
    flag_RE = opts_GLE.RelativeErr;
    optsKS.RE=flag_RE;
else
    flag_RE = 0;
    optsKS.RE = 0;
end
tol=opts_GLE.tol;
sigma_min = opts_GLE.sigma_min;
tol_L=opts_GLE.tol2;
tol_tSV=1e-1*tol; %Truncation tolerance
opts.maxit=1e3; opts.tol=1e-4; %SVD convergence options
n=size(A{1},1);
% Coefficent for Exit tolerance
c_SI=0.9; c_KS=0.1;
%% Solve for Reachability Gramian
if flag==1
    %% Computing Constant beta for the rescaling
    normD=0;
    for i=1:numel(A)
        Dfun3=@(v,flag) Dfun2(A{i},A{1},E{i},E{1},v,flag);
        normD=normD+svds(Dfun3,[n,n],1,'largest',opts)^2;
    end
    rho=(2*normD/(2*sigma_min));
    %% --------------------------------------------------------------
    nit=200; % Max number of iterations for KS
    %If norm of known term is larger than one then rescale with the
    %constant rho.
    tol_L = tol*sigma_min/(1+2*gamma);
    [Vs,Y] = solve_LSS_KS(A{1},E{1},B,tol_L,n,optsKS); 
    [U,Sig,~]=svd(Y,'econ');
    [s,~]=find(((diag(Sig)))<tol_L);
    if isempty(s)
        s=size(Y,2)+1;
    end
    s=s(1);
    Z_l=Vs*U(:,1:(s-1))*sqrt(Sig(1:(s-1),1:(s-1))); 
    maxit=50; % Maximum number of iterations for GLEs solver
    tol2=tol_L;
    %% Main Loop for GLE solver
    for j=1:maxit
        B_n=[];
        for i=2:numel(A)
            DP_l=Dfun(A{i}, A{1}, E{i}, E{1}, Z_l);
            B_n=[B_n, DP_l];
        end
        % Rescale the known term
        if rescale==1
            B_n=[B,(1/sqrt(rho))*B_n];
        else
            B_n=[B,B_n];
        end

        [Vs,Y] = solve_LSS_KS(A{1},E{1},B_n,tol2,nit,optsKS);
        Sig_old=Sig; U_old=U;
        [U,Sig,~]=svd(Y,'econ'); 
        [s,~]=find(((diag(Sig)))<tol_tSV); % Extract the relevant information from Y
        if isempty(s)
            [s,~]=find(((diag(Sig)./Sig(1,1)))<tol_tSV); % Extract the relevant information from Y
            if isempty(s)
                s=size(Y,2);
            end
        end
        s=s(1);
        Z_l_new=Vs*U(:,1:s)*sqrt(Sig(1:s,1:s));

        if size(Z_l_new,2)>size(Z_l,2)
            diff=size(Z_l_new,2)-size(Z_l,2);
            Sig_old=[Sig_old(1:size(Z_l,2),1:size(Z_l,2)),zeros(size(Z_l,2),diff); zeros(diff,size(Z_l_new,2))];
            U_old=[U_old(1:size(Z_l,2),1:size(Z_l,2)),zeros(size(Z_l,2),diff); zeros(diff,size(Z_l_new,2))];
            
            Z_l=[Z_l,zeros(size(Z_l,1),diff)];
        else
            Z_l=Z_l(:,1:size(Z_l_new,2));
        end
        
          GLE_ET=norm((abs(Z_l_new)-abs(Z_l))*sqrt(Sig(1:s,1:s)),'fro')...
              +norm((abs(Z_l)-(abs(Z_l_new)))*sqrt(Sig_old(1:s,1:s)),'fro');
          
      
        display(GLE_ET)
        % Exit tolerance for the method
        if flag_RE==1
            if GLE_ET<(c_SI*tol*Sig(1,1))/gamma

                fprintf('Solution of generalized Lypunov for Reachability Gramian converged \n')
                S=Z_l_new;
                break
            end
        else
            if GLE_ET<(c_SI*tol)/gamma

                fprintf('Solution of generalized Lypunov for Reachability Gramian converged \n')
                S=Z_l_new;
                break
            end
        end
        tol2=max(min((c_KS/c_SI)*GLE_ET,tol2),c_KS*2*tol*sigma_min/(1+2*gamma));
        Z_l=Z_l_new;
    end
    if j>=maxit
        
        fprintf('Solution of Generalized Lypunov for Reachability Gramian did not converge, Error Estimate is: %d \n',GLE_ET)
        S=Z_l_new;
    end
end
%% Solve for Observability Gramian
if flag==2
    %% Computing constant beta for the rescaling
    normD=0;
    for i=1:numel(A)
        Dfun4=@(v,flag) Dfun5(A{i},A{1},E{i},E{1},v,flag);
        normD=normD+svds(Dfun4,[n,n],1,'largest',opts)^2;
    end
    rho=(2*normD/(2*sigma_min));
    %% --------------------------------------------------------------
    nit=200; % Max number of iterations for KS
    [Vs,Y] = solve_LSS_KS_t(A{1},E{1},B',tol_L,n,optsKS); 
    [U,Sig,~]=svd(Y,'econ'); 
    [s,~]=find(((diag(Sig)))<tol_tSV); % Extract the relevant information from Y
    if isempty(s)
        s=size(Y,2);
    end
    s=s(1);
    S_l=Vs*U(:,1:s)*sqrt(Sig(1:s,1:s)); 
    maxit=100; 
    tol2=tol;
    for j=1:maxit
        B_n=[];
        for i=2:numel(A)
            DP_l=Dfun_t(A{i}, A{1}, E{i}, E{1}, S_l);
            B_n=[B_n, DP_l];
        end
        if rescale==1
            B_n=[B',sqrt(1/rho)*B_n];
        else
            B_n=[B',B_n];
        end

        [Vs,Y] = solve_LSS_KS_t(A{1},E{1},B_n,tol2,nit,optsKS);
        Sig_old=Sig; U_old=U;
        [U,Sig,~]=svd(Y,'econ'); [s,~]=find(((diag(Sig)))<tol_tSV);
     
        if isempty(s)
            s=size(Y,2);
        end
        s=s(1);
        S_l_new=Vs*U(:,1:s)*sqrt(Sig(1:s,1:s));

        if size(S_l_new,2)>size(S_l,2)
            diff=size(S_l_new,2)-size(S_l,2);
            Sig_old=[Sig_old(1:size(S_l,2),1:size(S_l,2)),zeros(size(S_l,2),diff); zeros(diff,size(S_l_new,2))];
            U_old=[U_old(1:size(S_l,2),1:size(S_l,2)),zeros(size(S_l,2),diff); zeros(diff,size(S_l_new,2))];
            

            S_l=[S_l,zeros(size(S_l,1),diff)];
        else
            S_l=S_l(:,1:size(S_l_new,2));
        end

        % Error Estimate
        GLE_ET=norm((abs(S_l_new)-abs(S_l))*(U(1:s,1:s))*sqrt(Sig(1:s,1:s)),'fro')...
              +norm((abs(S_l)-(abs(S_l_new)))*U_old(1:s,1:s)*sqrt(Sig_old(1:s,1:s)),'fro');
     
        display(GLE_ET)
        %% Exit criteria
        if flag_RE==0
            if GLE_ET<c_SI*tol/gamma

                fprintf('Solution of generalized Lypunov for Observability Gramian converged \n')
                S=S_l_new;
                break
            end
        else
            if GLE_ET<c_SI*Sig(1,1)*tol/gamma

                fprintf('Solution of generalized Lypunov for Observability Gramian converged \n')
                S=S_l_new;
                break
            end
        end
        tol2=max(min((c_KS/c_SI)*GLE_ET,tol2),c_KS*2*tol*sigma_min/(1+2*gamma));
        S_l=S_l_new;
    end
    if j>=maxit
        
        fprintf('Solution of Generalized Lypunov for Observability Gramian did not converge, error estimate is: %d \n',GLE_ET)
        S=S_l_new;
    end
end
end

%% -----------------------------------------------------------------------------------------------------------------------

function [Dv] = Dfun(A, A1, E, E1, v)

    Jiv = E\(A*v); J1v= E1\(A1*v);

    Dv=Jiv-J1v;

end

function [Dvt] = Dfun_t(A, A1, E, E1, v)
    
    Jivt= A'*(E'\v); J1vt= A1'*(E1'\v);
  
    Dvt=Jivt-J1vt;

end


function [v]=Dfun2(A,A1,E,E1,v,flag)
     switch flag
         case 'notransp'
             v=  Dfun(A, A1, E, E1, v);
         case 'transp'
             v=Dfun_t(A, A1, E, E1, v);
     end
end

function [v]=Dfun5(A,A1,E,E1,v,flag)
     switch flag
         case 'notransp'
             v=    Dfun_t(A, A1, E, E1, v);
         case 'transp'
             v=      Dfun(A, A1, E, E1, v);
     end
end

