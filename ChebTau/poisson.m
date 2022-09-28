clear;clc

% Simple chebtau for -au_xx-bu_yy=f with
% x-fourier z-cheb with dirichlet bcs
global Ny
Ny=2^6+1;
Nx=2^6;
ky=(0:Ny-1)';
kx=[0:Nx/2-1 -Nx/2 -Nx/2+1:-1]';
x=cos(pi*ky/(Ny-1));
y=2*pi*(-Nx/2:Nx/2-1)'/Nx;

[X,Y]=ndgrid(x,y);
[KX,KY]=ndgrid(kx,ky);

ue=(cosh(X)-cosh(1)).*exp(sin(Y));
RHS=-exp(sin(Y)).*(cosh(X).*(1+cos(Y).^2-sin(Y))+cosh(1)*(-cos(Y).^2+sin(Y)));

a=1;
b=1;
lambda=a*kx.^2/b;

f_tilde=fft(transpose(RHS))/b;
f_tilde=transpose(f_tilde);

% f_hat_tilde=fct(f_tilde);

u=zeros(size(X));

for j=1:Nx
    
    f_hat_tilde=fct2(f_tilde(:,j));
    
    A=zeros((Ny+1)/2,(Ny+1)/2);
    B=zeros((Ny+1)/2-1,(Ny+1)/2-1);
    
    g=zeros((Ny+1)/2,1);
    h=zeros((Ny+1)/2-1,1);
    
    % k=0 BC1
    A(1,:)=1;
    g(1)=0;

    % k=2:N-2

    % ODD POINTS

    for i=3:2:Ny-2
    
        A((i+1)/2,(i+1)/2)=1+lambda(j)*beta(i-1)/(2*((i-1)^2-1));
        A((i+1)/2,(i+1)/2-1)=-cc(i-3)*lambda(j)/(4*(i-1)*(i-2));
        A((i+1)/2,(i+1)/2+1)=-lambda(j)*beta(i+1)/(4*(i-1)*i);
    
        g((i+1)/2)=-cc(i-3)/(4*(i-1)*(i-2))*f_hat_tilde(i-2)+ ...
        beta(i-1)/(2*((i-1)^2-1))*f_hat_tilde(i)-beta(i+1)/(4*(i-1)*i)*f_hat_tilde(i+2);
    
    end

    % k=N
    A((Ny+1)/2,(Ny+1)/2)=1;
    A((Ny+1)/2,(Ny+1)/2-1)=-cc(Ny-3)*lambda(j)/(4*(Ny-1)*(Ny-2));
    g((Ny+1)/2)=-cc(Ny-3)/(4*(Ny-1)*(Ny-2))*f_hat_tilde(Ny-2);

    % EVEN POINTS

    % k=1 BC2
    B(1,:)=1;
    h(1)=0;

    for i=4:2:Ny-3
    
        B(i/2,i/2)=1+lambda(j)*beta(i-1)/(2*((i-1)^2-1));
        B(i/2,i/2-1)=-cc(i-3)*lambda(j)/(4*(i-1)*(i-2));
        B(i/2,i/2+1)=-lambda(j)*beta(i+1)/(4*(i-1)*i);
    
        h(i/2)=-cc(i-3)/(4*(i-1)*(i-2))*f_hat_tilde(i-2)+ ...
        beta(i-1)/(2*((i-1)^2-1))*f_hat_tilde(i)-beta(i+1)/(4*(i-1)*i)*f_hat_tilde(i+2);
    
    end

    % k=N-1
    B((Ny-1)/2,(Ny-1)/2)=1;
    B((Ny-1)/2,(Ny-1)/2-1)=-cc(Ny-4)*lambda(j)/(4*(Ny-2)*(Ny-3));
    h((Ny-1)/2)=-cc(Ny-4)/(4*(Ny-2)*(Ny-3))*f_hat_tilde(Ny-3);

    % Solve ...
    % ODD
    u1_hat=A\g;
    % EVEN
    u2_hat=B\h;

    % Combine:
    u_hat=zeros(Ny,1);
    u_hat(1:2:end)=u1_hat;
    u_hat(2:2:end)=u2_hat;

    % Transform back
    u(:,j)=ifct2(u_hat);

end

u=real(ifft(transpose(u)));
u=u';

% Residual
uxx=real(ifft(-kx.^2.*fft(u')));
uxx=uxx';

uyy=ifct(chebdiff(fct(u),2));

r=RHS+a*uxx+b*uyy;
disp(rms(r(:)))

function c=cc(k)
    if k==0
        c=2;
    elseif k>=1
        c=1;
    end

end

function b=beta(k)
global Ny    
    if k<=Ny-2 && k>=0
        b=1;
    elseif k>Ny-2
        b=0;
    end

end