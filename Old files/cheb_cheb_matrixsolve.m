% Solves au_xx+bu_yy+c*u=f with 2nd-order FD scheme (cheb-cheb mesh)
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

function v=cheb_cheb_matrixsolve(~,pde,domain,~)

Nx=domain.N(1);
Ny=domain.N(2);

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
% Cheb au_xx
% -------------------------------------------------------------------------

% 1D Cheb Dxx diff matrix
Dxx=ifct(chebdiff(fct(eye(Nx,Nx)),2));

% Dirichlet BCs
Dxx(1,:)=0;Dxx(end,:)=0;Dxx(1,1)=1;Dxx(end,end)=1;

% 2D Matrix
D2xx=kron(speye(Ny),Dxx);

% -------------------------------------------------------------------------
% Fourier bu_yy
% -------------------------------------------------------------------------

% 1D Cheb Dxx diff matrix
Dyy=ifct(chebdiff(fct(eye(Ny,Ny)),2));

% Dirichlet BCs
Dyy(1,:)=0;Dyy(end,:)=0;Dyy(1,1)=1;Dyy(end,end)=1;

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