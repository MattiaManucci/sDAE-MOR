function [lambda] = LMI_control(Z,A,B)
M=A*(Z*Z')+(Z*Z')*A'+B*B';
lambda=max(real(eig(M)));
end