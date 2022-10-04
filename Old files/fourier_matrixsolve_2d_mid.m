% Solves -(a*u_x)_x-(b*u_y)_y+c*u=f using Fourier matrix inversion using
% midpoints
%
% Inputs:
% pde.a
% pde.b
% pde.c
% pde.f
% N
% k
% 
% Ouputs:
% v - solution

function v=fourier_matrixsolve_2d_mid(~,pde,domain,option)

Nx=domain.N(1);
Ny=domain.N(2);
kx=domain.k{1};
ky=domain.k{2};

% Check if coefficients constant and find midpoint of a and b (at i+1)
if length(pde.a)==1
    
    a_mid=pde.a*ones(Nx,Ny);
    
else

    a_mid=real(ifft(fft(pde.a).*exp(pi*1i*kx/Nx)));

end

if length(pde.b)==1
    
    b_mid=pde.b*ones(Nx,Ny);
    
else
    
    b_mid=real(ifft(fft(pde.b').*exp(pi*1i*ky/Ny)));
    b_mid=b_mid';
    
end

% 1D Fourier Dx diff matrix at i+1/2
D1xf=ifft(1i.*kx.*fft(eye(Nx,Nx)).*exp(pi*1i*kx/Nx));
% 2D Fourier Dx diff matrix at i+1/2
D2xf=kron(speye(Ny),D1xf);

% 1D Fourier Dy diff matrix at i+1/2
D1yf=ifft(1i.*ky.*fft(eye(Ny,Ny)).*exp(pi*1i*ky/Ny));
% 2D Fourier Dy diff matrix at i+1/2
D2yf=kron(D1yf,speye(Nx));

% 1D Fourier Dx diff matrix at i-1/2
D1xb=ifft(1i.*kx.*fft(eye(Nx,Nx)).*exp(-pi*1i*kx/Nx));
% 2D Fourier Dx diff matrix at i-1/2
D2xb=kron(speye(Ny),D1xb);

% 1D Fourier Dy diff matrix at i-1/2
D1yb=ifft(1i.*ky.*fft(eye(Ny,Ny)).*exp(-pi*1i*ky/Ny));
% 2D Fourier Dy diff matrix at i-1/2
D2yb=kron(D1yb,speye(Nx));

% Solve + reshape
A=-D2xb*(a_mid(:).*D2xf)-D2yb*(b_mid(:).*D2yf)+pde.c(:).*speye(Nx*Ny);

% Check if Poisson type problem, then solve for mean 0 solution
if max(abs(pde.c(:)))<1e-12

    A(1,:)=1;
    pde.f(1)=0;
    
end

v=A\pde.f(:);
v=real(reshape(v,Nx,Ny));

end