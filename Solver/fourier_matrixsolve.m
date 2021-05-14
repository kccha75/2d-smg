% Solves -(a*u_x)_x-(b*u_y)_y+c*u=f using Fourier matrix inversion
%
% Inputs:
% pde.a
% pde.b
% pde.c
% pde.f
% scheme.kx
% scheme.ky
% 
% Ouputs:
% v - solution

function v=fourier_matrixsolve(v,pde,scheme,option)

[Nx,Ny]=size(v);

% 1D Fourier Dx diff matrix
D1x=ifft(1i.*scheme.kx.*fft(eye(Nx,Nx)));
% 2D Fourier Dx diff matrix
D2x=kron(speye(Ny),D1x);

% 1D Fourier Dy diff matrix
D1y=ifft(1i.*scheme.ky.*fft(eye(Ny,Ny)));
% 2D Fourier Dy diff matrix
D2y=kron(D1y,speye(Nx));

% Operator diff matrix
A=-D2x*(pde.a(:).*D2x)-D2y*(pde.b(:).*D2y)+pde.c(:).*eye(Nx*Ny);

% Check if Poisson type problem, then solve for mean 0 solution
if max(abs(pde.c(:)))<1e-12

    A(1,:)=1;
    pde.f(1)=0;
    
end

% Solve + reshape
v=A\pde.f(:);
v=real(reshape(v,Nx,Ny));

end