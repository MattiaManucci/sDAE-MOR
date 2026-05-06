function [V,W] = QWF(A, E, K, DAE_Index)
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
for i=1:DAE_Index
    V=A\(E*V);
end
end