% Function initialises PDE parameters
%
% Inputs:
%
% domain.dim - dimension of problem
% domain.X - ndgrid of x,y
% DJL.mu
% DJL.L
% DJL.u - wave speed
% DJL.KAI - x domain length
%
% Outputs:
%
% pde - structure of pde array coefficients
% domain.BC - boundary conditions (4x2 cell)
%

function [pde,domain]=DJLpdeinitialisecords(DJL,domain)

dim=domain.dim;
x=domain.x;
X=domain.X;

mu=DJL.mu;
pde.u=DJL.u;
KAI=DJL.KAI;


% -------------------------------------------------------------------------
% Set up PDE
% -------------------------------------------------------------------------
% PDE Parameters

Lx=KAI/mu/pi; % ?? domain to [-pi pi]
Ly=1/2; % [0 1] domain to [-1 1]

a=@(X,Y) 1/Lx^2;
b=@(X,Y) 1/Ly^2;
c=@(X,Y) 0*X;

% RHS
f=@(X,Y) 0*X;

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

a=a(X{1},X{2});
b=b(X{1},X{2});
c=c(X{1},X{2});
f=f(X{1},X{2});

pde.a = a;
pde.b = b;
pde.c = c;
pde.f = f;

% -------------------------------------------------------------------------
% Apply boundary conditions
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
        pde.f(index{1,i})=BCRHS{1,i}; % x(1) boundary
        pde.f(index{2,i})=BCRHS{2,i}; % x(end) boundary
    end
    
end

end