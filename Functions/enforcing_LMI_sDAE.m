function [Z_hat,LMI_Phat] = enforcing_LMI_sDAE(A,E,Z,B,rho_P)
%% Function to enforce LMI
n = size(A{1},1); nm = numel(A); AA = A{1}; optsKS = [];
Z_hat=Z;
RES = AA*(Z_hat*Z_hat')+(Z_hat*Z_hat')*(AA')+B*B';
for j=1:nm
    RES = RES + (1/(rho_P))*(A{j}-AA)*(Z_hat*Z_hat')*((A{j}-AA))';
end
[U,EIG] = eig(RES);

EIG=diag(EIG); ind=find(EIG>eps);
[Vs,Y] = solve_LSS_KS(AA,E,U(:,ind)*(diag(EIG(ind))).^(1/2),1e-15,1000,optsKS);
[D,DD] = svd(Y); DD=diag(DD); ind2=find((DD)<1e-15);
if isempty(ind2)
    ind2=size(Y,1);
end
Z_hat = [Z_hat, Vs*(D(:,1:ind2(1))*diag(DD(1:ind2(1)).^(0.5)))];
LMI_Phat(1,1)=LMI_control(Z_hat,AA,B);
LMI_Phat(1,2)=LMI_control(Z_hat,AA,sparse(n,1));
if LMI_Phat(1,1)>0
    fprintf('LMI is not satisfied for mode %d, and the associated approximated modified Gramian, largest eigenvalue is %d \n',1,LMI_Phat(1,1))
end
end