% Solves au_xx+bu_yy+c*u=f with 2nd-order FD scheme (cheb-fourier mesh)
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

function v=cheb_fourier_FD_2d(~,pde,domain,~)

Nx=domain.N(1);
Ny=domain.N(2);
dx=domain.dx{1};
dy=domain.dx{2};
a=pde.a;
b=pde.b;

% Check if coefficients constant
if length(pde.a)==1
    a=pde.a*ones(domain.N);
end
if length(pde.b)==1
    b=pde.b*ones(domain.N);
end

% -------------------------------------------------------------------------
% Cheb au_xx FD
% -------------------------------------------------------------------------

% Step size
dx1=dx(1:end-1);dx2=dx(2:end);

% Preset diagonal array
B=zeros(Nx,3);

% Main diag u_i
B(2:Nx-1,2)=-(1./dx1+1./dx2).*2./(dx1+dx2);

% Off diag u_i+1
B(3:Nx,3)=1./dx2.*2./(dx1+dx2);

% Off diag u_i-1
B(1:Nx-2,1)=1./dx1.*2./(dx1+dx2);

% Diagonal position array
d=[-1 0 1];

% Generate sparse matrix
Dxx=spdiags(B,d,Nx,Nx);

% Dirichlet BCs
Dxx(1,1)=1;Dxx(end,end)=1;

% 2D Matrix
D2xx=kron(speye(Ny),Dxx);

% -------------------------------------------------------------------------
% Fourier bu_yy FD
% -------------------------------------------------------------------------

% Preset diagonal array
B=zeros(Ny,3);

B(1:Ny,2)=(-2)/dy^2;
B(1:Ny-1,1)=(1)/(dy^2);
B(2:Ny,3)=(1)/(dy^2);

% Diagonal position array
d=[-1 0 1];

% Generate sparse matrix
Dyy=spdiags(B,d,Ny,Ny);

% First row (cyclic)
Dyy(1,Ny)=(1)/(dy^2);

% Last row (cyclic)
Dyy(Ny,1)=(1)/(dy^2);

% 2D Matrix
D2yy=kron(Dyy,speye(Nx));

% Solve + reshape
A=a(:).*D2xx+b(:).*D2yy+pde.c(:).*speye(Nx*Ny);

% Check if Poisson type problem, then solve for mean 0 solution
if max(abs(pde.c(:)))<1e-12

    A(1,:)=1;
    pde.f(1)=0;
    
end

v=A\pde.f(:);
v=reshape(v,Nx,Ny);

end