clear;close all;%clc
%
% -------------------------------------------------------------------------
% Solve PDE au_xxxx + bu_xx + cu + du^2= f using Fourier Cheb Spectral Multigrid
% Fifth-order KdV equation (see Yang)
% -------------------------------------------------------------------------
% INPUT PARAMETERS
% -------------------------------------------------------------------------

% Dimension of problem
dim=1;

% Discretisation flag for each dimension
% 1 - Fourier
% 2 - Cheb
discretisation=1;

% Boundary conditions for each discretisation (if fourier not used)
% x(1) a1*u+b1*u'=rhs1
alpha1{1}=@(y) 0;
beta1{1}=@(y) 1;
BCRHS1{1}=@(y) 0;

% x(end) a21*u+b21*u'= rhs21
alpha2{1}=@(y) 1;
beta2{1}=@(y) 0; 
BCRHS2{1}=@(y) 0;

% Grid size
finestgrid = 10;
coarsestgrid = 7;

% PDE Parameters
mu=-1.2;
a=@(X) 1/15^4;
b=@(X) 2/15^2;
c=@(X) -mu;

% RHS
f=@(X) 0*X;

% Initial guess
v0=@(X) 0.25*sech(0.3*X*15).*cos(X*15);

% -------------------------------------------------------------------------
% Multigrid Options here
% -------------------------------------------------------------------------

% Number of V-cycles if option is chosen, otherwise number of v-cycles done
% after FMG or number of decent iterations
option.numit=100;

% Solver / solution tolerance
option.tol=1e-10;

% Relaxations on the up and down cycle during Multigrid
option.Nd=1;
option.Nu=1;

% Multigrid solver options:'V-cycle' or 'FMG'
option.solver='FMG';

% Multigrid scheme: 'Correction' or 'FAS'
option.mgscheme='FAS';

% Operator, coarse grid solver, Relaxation
option.operator=@Lu_kdv;
option.coarsegridsolver=@cg;
option.relaxation=@MRR;

% Restriction for pde coefficients
option.restriction=@(vf) fourier_restrict_filtered(vf);

% Restriction for residual and RHS
option.restriction_residual=@(vf) fourier_restrict_filtered(vf);

% Prolongation
option.prolongation=@(vc) fourier_prolong_filtered(vc);

% Preconditioner
option.preconditioner=@yang_kdv_pre;

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

% -------------------------------------------------------------------------
% NEWTON HERE
% -------------------------------------------------------------------------

% New b(x) function in Newton
cnew=c+6*v0;
v=v0;

% Error guess (keep at 0)
e0=zeros(Nx,Ny);

tic
for i=1:20
    
    pde.c=c;
    % Initial RHS of linear equation
    pde.f=f-(option.operator(v,pde,domain)+3*v.^2);
    
    r=rms(rms(pde.f));
    fprintf('Residual Newton = %d\n',r)
    if r<=1e-10
        fprintf('Converged after %d Newton Iterations \n',i-1)
        break
    end
    
    % Solve linear equation
    pde.c=cnew;

%     option.tol=1e-1*r;
    [e,r]=cg(e0,pde,domain,option);
%     [e,r]=mg(v,pde,domain,option);

    % Update correction
    v=v+real(e);
    
    cnew=c+6*v;
    
end

if i==20
    
    fprintf('Did not converge to required tolerance after %d Newton Iterations\n',i)
    
end 
toc

plot(v);

