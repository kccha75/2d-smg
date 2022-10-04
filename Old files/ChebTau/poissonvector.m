% Simple poisson with lots of loops (bad bad) but works :)
% fixed the mess up with x and y ...

% tilde is FFT
% hat is FCT

clear;%clc
tic
for ii=1:10
% Simple chebtau for -au_xx-bu_yy=f with
% x-fourier z-cheb with dirichlet bcs

Nx=2^8;
Ny=2^8+1;

kx=[0:Nx/2-1 -Nx/2 -Nx/2+1:-1]';
ky=(0:Ny-1)';

x=2*pi*(-Nx/2:Nx/2-1)'/Nx;
y=cos(pi*ky/(Ny-1));

[X,Y]=ndgrid(x,y);
[KX,KY]=ndgrid(kx,ky);

ue=(cosh(Y)-cosh(1)).*exp(sin(X));
RHS=-exp(sin(X)).*(cosh(Y).*(1+cos(X).^2-sin(X))+cosh(1)*(-cos(X).^2+sin(X)));

a=1;
b=1;
lambda=a*kx.^2/b;

g=zeros(Nx,(Ny+1)/2);
h=zeros(Nx,(Ny+1)/2-1);

f_tilde=fft(RHS)/b;

u_tilde_hat=zeros(size(X));

f_tilde_hat=fct(transpose(f_tilde));

f_tilde_hat=transpose(f_tilde_hat);

%--------------------------------------------------------------------------
% ODD POINTS (even k, since matlab starts at 1)
%--------------------------------------------------------------------------

    % lower diag vector
    A1=zeros(Nx,(Ny+1)/2);
    % main diag vector
    A2=zeros(Nx,(Ny+1)/2);
    % upper diag vector
    A3=zeros(Nx,(Ny+1)/2);

    % k odd
    k_odd=(3:2:Ny-2);  

    A2(:,2:end-1)=1+lambda./(2*((k_odd-1).^2-1));
    A1(:,2:end-1)=-cc(k_odd-3).*lambda./(4*(k_odd-1).*(k_odd-2));
    A3(:,2:end-1)=-lambda./(4*(k_odd-1).*k_odd);
    
    % RHS
    g(:,2:end-1)=-cc(k_odd-3)./(4*(k_odd-1).*(k_odd-2)).*f_tilde_hat(:,k_odd-2)+ ...
        1./(2*((k_odd-1).^2-1)).*f_tilde_hat(:,k_odd)-1./(4*(k_odd-1).*k_odd).*f_tilde_hat(:,k_odd+2);
    
    % k=0 BC
    g(:,1)=0;
    A2(:,1)=1;
    A3(:,1)=1;

    % k=N BC
    A2(:,(Ny+1)/2)=1;
    A1(:,(Ny+1)/2)=-lambda./(4*(Ny-1)*(Ny-2));
    g(:,(Ny+1)/2)=-1/(4*(Ny-1)*(Ny-2))*f_tilde_hat(:,Ny-2);

    %--------------------------------------------------------------------------
    % EVEN POINTS (odd k, since matlab starts at 1)
    %--------------------------------------------------------------------------

    % lower diag
    B1=zeros(Nx,(Ny-1)/2);
    % main diag
    B2=zeros(Nx,(Ny-1)/2);
    % upper diag
    B3=zeros(Nx,(Ny-1)/2);

    % k even
    k_even=(4:2:Ny-3);
    
    B2(:,2:end-1)=1+lambda./(2*((k_even-1).^2-1));
    B1(:,2:end-1)=-lambda./(4*(k_even-1).*(k_even-2));
    B3(:,2:end-1)=-lambda./(4*(k_even-1).*k_even);
    
    % RHS
    h(:,2:end-1)=-1./(4*(k_even-1).*(k_even-2)).*f_tilde_hat(:,k_even-2)+ ...
        1./(2*((k_even-1).^2-1)).*f_tilde_hat(:,k_even)-1./(4*(k_even-1).*k_even).*f_tilde_hat(:,k_even+2);
    
    % k=1 BC
    h(:,1)=0;
    B2(:,1)=1;
    B3(:,1)=1;

    % k=N
    B2(:,(Ny-1)/2)=1;
    B1(:,(Ny-1)/2)=-lambda/(4*(Ny-1)*(Ny-2));
    g(:,(Ny-1)/2)=-1/(4*(Ny-1)*(Ny-2))*f_tilde_hat(:,Ny-2);

    % k=N-1
    B2(:,(Ny-1)/2)=1;
    B1(:,(Ny-1)/2)=-lambda/(4*(Ny-2)*(Ny-3));
    h(:,(Ny-1)/2)=-1/(4*(Ny-2)*(Ny-3))*f_tilde_hat(:,Ny-3);

    %--------------------------------------------------------------------------
    % Sort and generate sparse matrix
    %--------------------------------------------------------------------------

    % index row (main,lower,upper)
    rowA=[(1:(Ny+1)/2)';(2:(Ny+1)/2)';(1:(Ny+1)/2-1)';ones((Ny+1)/2-2,1)];
    
    % index columns
    columnA=[(1:(Ny+1)/2)';(1:(Ny+1)/2-1)';(2:(Ny+1)/2)';(3:(Ny+1)/2)'];
    
    % matrix entry
    entryA=[A2 A1(:,2:end) A3(:,1:end-1) ones(Nx,(Ny+1)/2-2)];
    
    % index row (main,lower,upper)
    rowB=[(1:(Ny-1)/2)';(2:(Ny-1)/2)';(1:(Ny-1)/2-1)';ones((Ny-1)/2-2,1)];

    % index columns
    columnB=[(1:(Ny-1)/2)';(1:(Ny-1)/2-1)';(2:(Ny-1)/2)';(3:(Ny-1)/2)'];

    % matrix entry
    entryB=[B2 B1(:,2:end) B3(:,1:end-1) ones(Nx,(Ny-1)/2-2)];


    
%--------------------------------------------------------------------------
% Solve ...
%--------------------------------------------------------------------------

for j=1:Nx

    % ODD
    u1_tilde_hat=sparse(rowA,columnA,entryA(j,:),(Ny+1)/2,(Ny+1)/2)\transpose(g(j,:));
    % EVEN
    u2_tilde_hat=sparse(rowB,columnB,entryB(j,:),(Ny-1)/2,(Ny-1)/2)\transpose(h(j,:));

    % Combine:
    u_tilde_hat(j,1:2:end)=u1_tilde_hat;
    u_tilde_hat(j,2:2:end)=u2_tilde_hat;

end

% Inverse transforms
u=real(ifft(transpose(ifct(transpose(u_tilde_hat)))));

% Residual
uxx=real(ifft(-kx.^2.*fft(u)));

uyy=ifct(chebdiff(fct(transpose(u)),2));
uyy=transpose(uyy);

r=RHS+a*uxx+b*uyy;
% disp(rms(r(:)))
end
toc

function c=cc(k)

    c=2*(k==0)+(k>=1);

end