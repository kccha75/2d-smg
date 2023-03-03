% Function initialises PDE parameters
%
% Inputs:
%
% domain.dim - dimension of problem
% domain.X - ndgrid of x,y
% fKdV.L
% fKdV.gamma
% fKdV.delta
%
% Outputs:
%
% pde - structure of pde array coefficients
% domain.BC - boundary conditions (4x2 cell)
%

function [fKdV,pde,domain]=fKdVpdeinitialise(fKdV,domain)

dim=domain.dim;
x=domain.x;
X=domain.X;

L=fKdV.L;
gamma=fKdV.gamma;
delta=fKdV.delta;
topography=fKdV.topography;

% -------------------------------------------------------------------------
% Set up PDE
% -------------------------------------------------------------------------
% PDE Parameters

a=@(X) (2*pi/L)^2;
b=@(X) 0;
c=@(X) -delta;

% RHS
f=@(X) -gamma*topography(L/(2*pi)*X);

% -------------------------------------------------------------------------
% Set up BCs
% -------------------------------------------------------------------------
% Boundary conditions for each discretisation (Fourier not used)

% Boundary conditions for each discretisation (if fourier not used)
% x(1) a11*u+b11*u'=rhs11 
alpha1{1}=@(y) 1;
beta1{1}=@(y) 0;
BCRHS1{1}=@(y) 0*y;

% x(end) a21*u+b21*u'= rhs21
alpha2{1}=@(y) 1;
beta2{1}=@(y) 0; 
BCRHS2{1}=@(y) 0*y;

% y(1) a12*u+b12*u'=rhs12 
alpha1{2}=@(x) 1;
beta1{2}=@(x) 0;
BCRHS1{2}=@(x) 0*x;

% y(end) a22*u+b22*u'= rhs22
alpha2{2}=@(x) 1;
beta2{2}=@(x) 0;
BCRHS2{2}=@(x) 0*x;

BC=cell(4,dim);
BCRHS=cell(2,dim);

y=fliplr(x);

for i=1:dim
    
    BC{1,i}=alpha1{i}(x{i});
    BC{2,i}=beta1{i}(x{i});
    BC{3,i}=alpha2{i}(x{i});
    BC{4,i}=beta2{i}(x{i});

    BCRHS{1,i}=BCRHS1{i}(y{i});
    BCRHS{2,i}=BCRHS2{i}(y{i});
    
end

% -------------------------------------------------------------------------
% Sort into structures
% -------------------------------------------------------------------------
domain.BC = BC;

a=a(X);
b=b(X);
c=c(X);
f=f(X);

pde.a = a;
pde.b = b;
pde.c = c;
pde.f = f;

% -------------------------------------------------------------------------
% Apply boundary conditions nonlinear, DIRICHLET
% -------------------------------------------------------------------------
index=cell(2,domain.dim); % Left and right, for each dimension

Nx=domain.N(1);
Ny=domain.N(2);

index{1,1}=1:Nx:Nx*(Ny-1)+1; % Top boundary of matrix (x(1))
index{2,1}=Nx:Nx:Nx*Ny; % Bottom boundary of matrix (x(end))

index{1,2}=1:Nx; % Left boundary of matrix (y(1))
index{2,2}=Nx*(Ny-1)+1:Nx*Ny; % Right boundary of matrix (y(end))

for i=1:domain.dim
    
    if domain.discretisation(i)~=1 % If not Fourier, set BCs
        fKdV.v(index{1,i})=BCRHS{1,i}; % x(1) boundary
        fKdV.v(index{2,i})=BCRHS{2,i}; % x(end) boundary
    end
    
end

end