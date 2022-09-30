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

% cc=2*(k==0)+(k>=1);
% k=2:N-2

%--------------------------------------------------------------------------
% ODD POINTS
%--------------------------------------------------------------------------

% for i=3:2:N-2

index_odd=(3:2:N-2)';
    
    A(2:end-1,2)=1+lambda./(2*((index_odd-1).^2-1));
    A(1:end-2,1)=-cc(index_odd-3)*lambda./(4*(index_odd-1).*(index_odd-2));
    A(3:end,3)=-lambda./(4*(index_odd-1).*index_odd);
    
    g(2:end-1)=-cc(index_odd-3)./(4*(index_odd-1).*(index_odd-2)).*f_hat(index_odd-2)+ ...
        1./(2*((index_odd-1).^2-1)).*f_hat(index_odd)-1./(4*(index_odd-1).*index_odd).*f_hat(index_odd+2);
    
    AA=spdiags(A,-1:1,(N+1)/2,(N+1)/2);

% k=0 BC1
g(1)=0;
AA(1,:)=1;
% end

% k=N
AA((N+1)/2,(N+1)/2)=1;
AA((N+1)/2,(N+1)/2-1)=-lambda/(4*(N-1)*(N-2));
g((N+1)/2)=-1/(4*(N-1)*(N-2))*f_hat(N-2);

%--------------------------------------------------------------------------
% EVEN POINTS
%--------------------------------------------------------------------------

index_even=(4:2:N-3)';
    
B(2:end-1,2)=1+lambda./(2*((index_even-1).^2-1));
B(1:end-2,1)=-lambda./(4*(index_even-1).*(index_even-2));
B(3:end,3)=-lambda/(4*(index_even-1).*index_even);
    
h(2:end-1)=-1./(4*(index_even-1).*(index_even-2)).*f_hat(index_even-2)+ ...
        1./(2*((index_even-1).^2-1)).*f_hat(index_even)-1./(4*(index_even-1).*index_even).*f_hat(index_even+2);
    
BB=spdiags(B,-1:1,(N+1)/2-1,(N+1)/2-1);
% k=1 BC2
BB(1,:)=1;
h(1)=0;

% k=N-1
BB((N-1)/2,(N-1)/2)=1;
BB((N-1)/2,(N-1)/2-1)=-lambda/(4*(N-2)*(N-3));
h((N-1)/2)=-1/(4*(N-2)*(N-3))*f_hat(N-3);

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
