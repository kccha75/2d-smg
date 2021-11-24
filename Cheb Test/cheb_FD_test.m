function v=cheb_FD_test(~,pde,domain,~)

Nx=domain.N(1);
dx=domain.dx{1};


% Check if coefficients constant
if length(pde.a)==1
    a=pde.a*ones(domain.N+1);
end
if length(pde.b)==1
    b=pde.b*ones(domain.N+1);
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

A=pde.a(:).*Dxx+speye(Nx);

% Dirichlet BCs
A(1,1)=1;A(end,end)=1;

v=A\pde.f(:);
end