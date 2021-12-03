% Solves au_xx+bu_yy+c*u=f with 2nd-order FD scheme (cheb-fourier mesh)
% WITH NEUMANN
%
% Uses matrix inversion (backslash)
% Assumes periodic boundary conditions
%
% Inputs:
% v - not used except for size purposes
% pde.a
% pde.b
% pde.f
% domain.dx
% domain.N
% 
% Ouputs:
% v - solution
%

function v=cheb_fourier_matrixsolve_neumann(~,pde,domain,~)

Nx=domain.N(1);
Ny=domain.N(2);
kx=domain.k{1};
ky=domain.k{2};
a=pde.a;
b=pde.b;
c=pde.c;

% Check if coefficients constant
if length(pde.a)==1
    a=pde.a*ones(domain.N);
end
if length(pde.b)==1
    b=pde.b*ones(domain.N);
end
if length(pde.c)==1
    c=pde.c*ones(domain.N);
end

% -------------------------------------------------------------------------
% Cheb au_xx
% -------------------------------------------------------------------------

% 1D Cheb Dxx diff matrix
Dxx=ifct(chebdiff(fct(eye(Nx,Nx)),2));

% BCs
% Right Dirichlet
Dxx(1,:)=0;
Dxx(1,:)=sum(fct(eye(Nx,Nx)).*kx.^2);

% Left Neumann
Dxx(end,:)=0;Dxx(end,end)=1;

a(1,:)=1;a(end,:)=1;
b(1,:)=0;b(end,:)=0;
c(1,:)=0;c(end,:)=0;

% 2D Matrix
D2xx=kron(speye(Ny),Dxx);

% -------------------------------------------------------------------------
% Fourier bu_yy
% -------------------------------------------------------------------------

% 1D Fourier Dyy diff matrix
Dyy=real(ifft(-ky.^2.*fft(eye(Ny,Ny))));

% 2D Matrix
D2yy=kron(Dyy,speye(Nx));

% -------------------------------------------------------------------------
% Spectral matrix
A=a(:).*D2xx+b(:).*D2yy+c(:).*speye(Nx*Ny);

% -------------------------------------------------------------------------
% Boundary conditions
% -------------------------------------------------------------------------

% Left

% A(1:Nx)=

% Right

% A((Ny-1)*Nx+1:Nx*Ny)=

% % Top
% A(1:Nx:(Ny-1)*Nx+1,:)=sparse(Ny,Nx*Ny); % Zero specific rows to BCs
% index=sub2ind(size(A),1:Nx:(Ny-1)*Nx+1,1:Nx:(Ny-1)*Nx+1); % Get index
% A(index)=1; % Set main diag element of rows to 1
% 
% % Bottom
% A(Nx:Nx:Nx*Ny,:)=sparse(Ny,Nx*Ny);
% index=sub2ind(size(A),Nx:Nx:Nx*Ny,Nx:Nx:Nx*Ny);
% A(index)=1;

% -------------------------------------------------------------------------
% Check if Poisson type problem, then solve for mean 0 solution
if max(abs(pde.c(:)))<1e-12

    A(1,:)=1;
    pde.f(1)=0;
    
end

v=A\pde.f(:);
v=reshape(v,Nx,Ny);

end