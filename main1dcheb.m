clear;close all;%clc
%
% Simple main for Cheb Fourier SMG, uses homogenous Dirichlet / Neumann BCs
% 
% -------------------------------------------------------------------------
% Solve PDE au_xx + bu_x + cu = f using Fourier Cheb Spectral Multigrid
% -------------------------------------------------------------------------
% INPUT PARAMETERS
% -------------------------------------------------------------------------

% Dimension of problem
dim=1;

% Discretisation flag for each dimension
% 1 - Fourier
% 2 - Cheb
discretisation=2;

% Boundary conditions for each discretisation (if fourier not used)
% x(1) a1*u+b1*u'=rhs1
alpha1{1}=@(y) 0;
beta1{1}=@(y) 1;
BCRHS1{1}=@(y) 2*exp(1);

% x(end) a21*u+b21*u'= rhs21
alpha2{1}=@(y) 1;
beta2{1}=@(y) 0; 
BCRHS2{1}=@(y) -1/exp(1);

% Grid size
finestgrid = 6;
coarsestgrid = 3;

% PDE Parameters
a=@(X) 2+sin(X);
b=@(X) 1+X.^2;
c=@(X) exp(X);

% RHS
f=@(X) exp(X).*(5+X.*(3+exp(X)+X+X.^2)+(2+X).*sin(X));

% Exact solution
ue=@(X) X.*exp(X);

% Initial guess
v0=@(X) rand(size(X));

% -------------------------------------------------------------------------
% Multigrid Options here
% -------------------------------------------------------------------------

% Number of V-cycles if option is chosen, otherwise number of v-cycles done
% after FMG
option.num_vcycles=1;

% Solver / solution tolerance
option.tol=1e-12;

% Relaxations on the up and down cycle during Multigrid
option.Nd=1;
option.Nu=1;

% Multigrid solver options:'V-cycle' or 'FMG'
option.solver='FMG';

% Multigrid scheme: 'Correction' or 'FAS'
option.mgscheme='FAS';

% Operator, coarse grid solver, Relaxation
option.operator=@Lu;
option.coarsegridsolver=@specmatrixsolve;
option.relaxation=@MRR;

% Restriction for pde coefficients
option.restriction=@(vf) cheb_restrict(vf);

% Restriction for residual and RHS
option.restriction_residual=@(vf) cheb_restrict_residual(vf);

% Prolongation
option.prolongation=@(vc) cheb_prolong(vc);

% Preconditioner
option.preconditioner=@FDmatrixsolve;

% Number of preconditioned relaxations
option.prenumit=1;

% -------------------------------------------------------------------------
% Set up parameters
% -------------------------------------------------------------------------
N=ones(1,max(dim,2)); 
x=cell(1,dim);
k=cell(1,dim);
dx=cell(1,dim);

for i=1:length(discretisation)
    
    switch discretisation(i)

        % Fourier discretisation
        case 1
            N(i) = 2^finestgrid;
            k{i} = [0:N(i)/2-1 -N(i)/2 -N(i)/2+1:-1]';
            x{i} = 2*pi*(-N(i)/2:N(i)/2-1)'/N(i);
            dx{i} = x{i}(2)-x{i}(1);
            
       % Chebyshev discretisation
        case 2
            N(i) = 2^finestgrid+1;
            k{i} = (0:N(i)-1)';
            x{i} = cos(pi*k{i}/(N(i)-1));
            dx{i} = x{i}(2:end)-x{i}(1:end-1); % is negative but ok :)
            
    end
    
end

X=x{1};

a=a(X);
b=b(X);
c=c(X);
f=f(X);

ue=ue(X);
v0=v0(X);

% -------------------------------------------------------------------------
% Set up BCs
% -------------------------------------------------------------------------
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
domain.dim = dim;
domain.discretisation = discretisation;
domain.BC = BC;

domain.N = N;
domain.k = k;
domain.dx = dx;

pde.a = a;
pde.b = b;
pde.c = c;
pde.f = f;

option.finestgrid=finestgrid;
option.coarsestgrid=coarsestgrid;
option.grids=finestgrid-coarsestgrid+1;

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

% -------------------------------------------------------------------------
% SOLVE HERE
% -------------------------------------------------------------------------

tic
[v,r]=mg(v0,pde,domain,option);
% option.numit=30;
% [v,r]=MRR(v0,pde,domain,option);
% [v,r]=bicgstab(v0,pde,domain,option);
toc
disp(rms(r(:)))


plot(X,v);xlabel('x');ylabel('v');title('Numerical solution of ODE')
figure;plot(X,abs(v-ue));xlabel('x');ylabel('Error');title('Error compared to exact solution')

