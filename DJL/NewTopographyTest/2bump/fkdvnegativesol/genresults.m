% Just saving some fkdv sol, truncated from L=200 to L=50
% note l=8 here ...

load('L200l8.mat')

B_cont_neg=v(2^12+1-2^10:2^12+2^10,:);
gamma_cont_neg=gamma;
lambda_cont_neg=lambda;
save('gamma_cont_neg.mat','gamma_cont_neg');
save('lambda_cont_neg.mat','lambda_cont_neg');
save('B_cont_neg.mat','B_cont_neg');