clear;clc
tic
for ii=1:10000
% Simple chebtau for -u_xx+lambda*u=f with dirichlet BCs u(1)=u(-1)=0
global N
N=2^8+1;
k=(0:N-1)';
x=cos(pi*k/(N-1));

ue=exp(x.^2)-exp(1);
lambda=1;
RHS=-exp(1)-exp(x.^2).*(1+4*x.^2);

f_hat=real(fct(RHS));

A=zeros((N+1)/2,3);
B=zeros((N+1)/2-1,3);

g=zeros((N+1)/2,1);
h=zeros((N+1)/2-1,1);

%--------------------------------------------------------------------------
% ODD POINTS (even k, since matlab starts at 1)
%--------------------------------------------------------------------------

% lower diag vector
A1=zeros((N+1)/2,1);
% main diag vector
A2=zeros((N+1)/2,1);
% upper diag vector
A3=zeros((N+1)/2,1);

% k odd
k_odd=(3:2:N-2)';  

A2(2:end-1)=1+lambda./(2*((k_odd-1).^2-1));
A1(2:end-1)=-cc(k_odd-3)*lambda./(4*(k_odd-1).*(k_odd-2));
A3(2:end-1)=-lambda./(4*(k_odd-1).*k_odd);
    
% RHS
g(2:end-1)=-cc(k_odd-3)./(4*(k_odd-1).*(k_odd-2)).*f_hat(k_odd-2)+ ...
        1./(2*((k_odd-1).^2-1)).*f_hat(k_odd)-1./(4*(k_odd-1).*k_odd).*f_hat(k_odd+2);
    
% k=0 BC
g(1)=0;
A2(1)=1;
A3(1)=1;

% k=N BC
A2((N+1)/2)=1;
A1((N+1)/2)=-lambda/(4*(N-1)*(N-2));
g((N+1)/2)=-1/(4*(N-1)*(N-2))*f_hat(N-2);

%--------------------------------------------------------------------------
% EVEN POINTS (odd k, since matlab starts at 1)
%--------------------------------------------------------------------------

% lower diag
B1=zeros((N-1)/2,1);
% main diag
B2=zeros((N-1)/2,1);
% upper diag
B3=zeros((N-1)/2,1);

% k even
k_even=(4:2:N-3)';
    
B2(2:end-1)=1+lambda./(2*((k_even-1).^2-1));
B1(2:end-1)=-lambda./(4*(k_even-1).*(k_even-2));
B3(2:end-1)=-lambda/(4*(k_even-1).*k_even);

% RHS
h(2:end-1)=-1./(4*(k_even-1).*(k_even-2)).*f_hat(k_even-2)+ ...
        1./(2*((k_even-1).^2-1)).*f_hat(k_even)-1./(4*(k_even-1).*k_even).*f_hat(k_even+2);
    
% k=1 BC
h(1)=0;
B2(1)=1;
B3(1)=1;

% k=N
B2((N-1)/2)=1;
B1((N-1)/2)=-lambda/(4*(N-1)*(N-2));
g((N-1)/2)=-1/(4*(N-1)*(N-2))*f_hat(N-2);

% k=N-1
B2((N-1)/2)=1;
B1((N-1)/2)=-lambda/(4*(N-2)*(N-3));
h((N-1)/2)=-1/(4*(N-2)*(N-3))*f_hat(N-3);


% index row (main,lower,upper)
row=[(1:(N+1)/2)';(2:(N+1)/2)';(1:(N+1)/2-1)';ones((N+1)/2-2,1)];

% index columns
column=[(1:(N+1)/2)';(1:(N+1)/2-1)';(2:(N+1)/2)';(3:(N+1)/2)'];

% matrix entry
entry=[A2;A1(2:end);A3(1:end-1);ones((N+1)/2-2,1)];

% Sparse matrix
AA=sparse(row,column,entry,(N+1)/2,(N+1)/2);

% Generate matrix
% index row (main,lower,upper)
row=[(1:(N-1)/2)';(2:(N-1)/2)';(1:(N-1)/2-1)';ones((N-1)/2-2,1)];

% index columns
column=[(1:(N-1)/2)';(1:(N-1)/2-1)';(2:(N-1)/2)';(3:(N-1)/2)'];

% matrix entry
entry=[B2;B1(2:end);B3(1:end-1);ones((N-1)/2-2,1)];

% Sparse matrix
BB=sparse(row,column,entry,(N-1)/2,(N-1)/2);

% Solve ...
% ODD
u1_hat=AA\g;
% EVEN
u2_hat=BB\h;

% Combine:
u_hat=zeros(N,1);
u_hat(1:2:end)=u1_hat;
u_hat(2:2:end)=u2_hat;


% Transform back
u=ifct(u_hat);

% Check residual ...
r=RHS+ifct(chebdiff(fct(u),2))-lambda*u;
% disp(rms(r))

end
toc
% Plots!
% plot(x,u-ue);title('error plot')

function c=cc(k)

    c=2*(k==0)+(k>=1);

end
