%% Code based on the work:
%[1] M. Manucci and B. Unger, ...
% ArXiv e-print ..., 2024.
%% --------------------------------------------------------------
function [V,Y] = solve_LSS_KS_t(A,E,B,tol,k,optsKS)
if isfield(optsKS,'RE')
    flagRE = optsKS.RE;
else
    flagRE = 0;
end
if flagRE==0
% Iteration 0
Vk=full(B); tolGS=1e-14; V=B(:,1)/norm(B(:,1));
for ii=2:size(Vk,2)
    
        Vk(:,ii)=Vk(:,ii) - V*(V'*Vk(:,ii)); 
        Vk(:,ii)=Vk(:,ii) - V*(V'*Vk(:,ii));
        Vk(:,ii)=Vk(:,ii) - V*(V'*Vk(:,ii)); 
    
    if norm(Vk(:,ii))>tolGS
        Vk(:,ii)=Vk(:,ii)/norm(Vk(:,ii));
        V=[V,Vk(:,ii)];
    end
end
V0=V; V=V0; n=size(V,2); AV=A'*(E'\V); VAV=V'*AV; BV=V'*B;
for i=1:k
 
    Vk=A'*(E'\V0);

    for ii=1:size(Vk,2)

        Vk(:,ii)=Vk(:,ii) - V*(V'*Vk(:,ii)); 
        Vk(:,ii)=Vk(:,ii) - V*(V'*Vk(:,ii));
        Vk(:,ii)=Vk(:,ii) - V*(V'*Vk(:,ii)); 

        if norm(Vk(:,ii))>tolGS
            Vk(:,ii)=Vk(:,ii)/norm(Vk(:,ii));
            V=[V,Vk(:,ii)];
        end

    end
    n_new=size(V,2);
    if i>1
       Vj=V(:,n+1:n_new);
       Res=2*Vj'*(AV*Y);
       KS_RES=norm(Res,"fro");
       if KS_RES<tol
           V=V(:,1:n);
           break
       end
    end
 
    V0=V(:,n+1:n_new);  AVnew=A'*(E'\V0);
    
    VAV=[       VAV,        V(:,1:n)'*AVnew ;...
        V(:,n+1:n_new)'*AV ,V(:,n+1:n_new)'*AVnew];

    AV=[AV,AVnew];
    BV=[BV;V0'*B]; 
    
    Y=lyap(VAV,BV*BV');
    n=n_new;
end
else
   % Iteration 0
Vk=full(B); tolGS=1e-14; V=B(:,1)/norm(B(:,1));
for ii=2:size(Vk,2)
    
        Vk(:,ii)=Vk(:,ii) - V*(V'*Vk(:,ii)); 
        Vk(:,ii)=Vk(:,ii) - V*(V'*Vk(:,ii));
        Vk(:,ii)=Vk(:,ii) - V*(V'*Vk(:,ii)); 
    
    if norm(Vk(:,ii))>tolGS
        Vk(:,ii)=Vk(:,ii)/norm(Vk(:,ii));
        V=[V,Vk(:,ii)];
    end
end
V0=V; V=V0; n=size(V,2); AV=A'*(E'\V); VAV=V'*AV; BV=V'*B;
for i=1:k
 
    Vk=A'*(E'\V0);

    for ii=1:size(Vk,2)

        Vk(:,ii)=Vk(:,ii) - V*(V'*Vk(:,ii)); 
        Vk(:,ii)=Vk(:,ii) - V*(V'*Vk(:,ii));
        Vk(:,ii)=Vk(:,ii) - V*(V'*Vk(:,ii)); 

        if norm(Vk(:,ii))>tolGS
            Vk(:,ii)=Vk(:,ii)/norm(Vk(:,ii));
            V=[V,Vk(:,ii)];
        end

    end
    n_new=size(V,2);
    if i>1
       Vj=V(:,n+1:n_new);
       Res=2*Vj'*(AV*Y);
       KS_RES=norm(Res,"fro");
       if KS_RES<tol*norm(Y)
           V=V(:,1:n);
           break
       end
    end
 
    V0=V(:,n+1:n_new);  AVnew=A'*(E'\V0);
    
    VAV=[       VAV,        V(:,1:n)'*AVnew ;...
        V(:,n+1:n_new)'*AV ,V(:,n+1:n_new)'*AVnew];

    AV=[AV,AVnew];
    BV=[BV;V0'*B]; 
    
    Y=lyap(VAV,BV*BV');
    n=n_new;
end 
end
end