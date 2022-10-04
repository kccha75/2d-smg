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

A=zeros((N+1)/2,(N+1)/2);
B=zeros((N+1)/2-1,(N+1)/2-1);
g=zeros((N+1)/2,1);
h=zeros((N+1)/2-1,1);

% k=0 BC1
A(1,:)=1;
g(1)=0;

% k=2:N-2

% ODD POINTS

for i=3:2:N-2
    
    A((i+1)/2,(i+1)/2)=1+lambda*beta(i-1)/(2*((i-1)^2-1));
    A((i+1)/2,(i+1)/2-1)=-cc(i-3)*lambda/(4*(i-1)*(i-2));
    A((i+1)/2,(i+1)/2+1)=-lambda*beta(i+1)/(4*(i-1)*i);
    
    g((i+1)/2)=-cc(i-3)/(4*(i-1)*(i-2))*f_hat(i-2)+ ...
        beta(i-1)/(2*((i-1)^2-1))*f_hat(i)-beta(i+1)/(4*(i-1)*i)*f_hat(i+2);
    
end

% k=N
A((N+1)/2,(N+1)/2)=1;
A((N+1)/2,(N+1)/2-1)=-cc(N-3)*lambda/(4*(N-1)*(N-2));
g((N+1)/2)=-cc(N-3)/(4*(N-1)*(N-2))*f_hat(N-2);

% EVEN POINTS

% k=1 BC2
B(1,:)=1;
h(1)=0;

for i=4:2:N-3
    
    B(i/2,i/2)=1+lambda*beta(i-1)/(2*((i-1)^2-1));
    B(i/2,i/2-1)=-cc(i-3)*lambda/(4*(i-1)*(i-2));
    B(i/2,i/2+1)=-lambda*beta(i+1)/(4*(i-1)*i);
    
    h(i/2)=-cc(i-3)/(4*(i-1)*(i-2))*f_hat(i-2)+ ...
        beta(i-1)/(2*((i-1)^2-1))*f_hat(i)-beta(i+1)/(4*(i-1)*i)*f_hat(i+2);
    
end

% k=N-1
B((N-1)/2,(N-1)/2)=1;
B((N-1)/2,(N-1)/2-1)=-cc(N-4)*lambda/(4*(N-2)*(N-3));
h((N-1)/2)=-cc(N-4)/(4*(N-2)*(N-3))*f_hat(N-3);

% Solve ...
% ODD
u1_hat=A\g;
% EVEN
u2_hat=B\h;

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
    if k==0
        c=2;
    elseif k>=1
        c=1;
    end

end

function b=beta(k)
global N    
    if k<=N-2 && k>=0
        b=1;
    elseif k>N-2
        b=0;
    end

end