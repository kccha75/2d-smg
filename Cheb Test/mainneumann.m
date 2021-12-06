clear;close all;%clc
% -------------------------------------------------------------------------
% Solve PDE au_xx + bu_yy + cu = f using Fourier Cheb Spectral Multigrid
% -------------------------------------------------------------------------
% INPUT PARAMETERS
% -------------------------------------------------------------------------

% Dimension of problem
dim=2;

% Discretisation flag for each dimension
% 1 - Fourier
% 2 - Cheb
discretisation=[2 1];

% Boundary flags (row) (x(1) x(end)) for each dimension
% 1 - Dirichlet
% 2 - Neumann
% 3 - Fourier
BCflag=[2 1 3 3];

% Boundary condition values (column vector) (Non-fourier only)
BC=[0 0];

finestgrid = 5;
coarsestgrid = 3;

% PDE Parameters
a=@(X,Y) 1;
b=@(X,Y) 1;
c=@(X,Y) 1;

% RHS
f=@(X,Y) 1/4*exp(sin(Y)).*(cosh(1/2*(-1+X)).*(5+4*cos(Y).^2-4*sin(Y))-4*cosh(1)*(1+cos(Y).^2-sin(Y)));

% Exact solution
ue=@(X,Y) (cosh(1/2*(X-1))-cosh(1)).*exp(sin(Y));

% Initial guess
v0=@(X,Y) 0*X;

% -------------------------------------------------------------------------
% Multigrid Options here
% -------------------------------------------------------------------------

% Number of V-cycles if option is chosen, otherwise number of v-cycles done
% after FMG
option.num_vcycles=5;

% Solver / solution tolerance
option.tol=1e-12;

% Relaxations on the up and down cycle during Multigrid
option.Nd=1;
option.Nu=1;

% Multigrid solver options:'V-cycle' or 'FMG'
option.solver='V-cycle';

% Multigrid scheme: 'Correction' or 'FAS'
option.mgscheme='Correction';

% Operator, coarse grid solver, Relaxation
option.operator=@Lu_2d;
option.coarsegridsolver=@cheb_fourier_matrixsolve;
option.relaxation=@MRR;

% Restriction
x_restrict=@cheb_restrict_dirichlet;
y_restrict=@fourier_restrict_filtered;
option.restriction=@(vf) restrict_2d(vf,x_restrict,y_restrict);

% Prolongation
x_prolong=@cheb_prolong;
y_prolong=@fourier_prolong_filtered;
option.prolongation=@(vc) prolong_2d(vc,x_prolong,y_prolong);

% Preconditioner
option.preconditioner=@cheb_fourier_FD_neumann;
% Number of preconditioned relaxations
option.prenumit=1;

% -------------------------------------------------------------------------
% Set up parameters
% -------------------------------------------------------------------------
N=zeros(1,dim);
x=cell(1,dim);
k=x;
dx=x;

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
            dx{i} = x{i}(2:end)-x{i}(1:end-1);
            
    end
    
end

[X,Y] = ndgrid(x{1},x{2});

a=a(X,Y);
b=b(X,Y);
c=c(X,Y);
f=f(X,Y);

ue=ue(X,Y);
v0=v0(X,Y);

% -------------------------------------------------------------------------
% Sort into structures
% -------------------------------------------------------------------------
domain.dim = dim;
domain.discretisation = discretisation;
domain.BCflag = BCflag;
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

index=cell(length(domain.BCflag),1);

Nx=domain.N(1);
Ny=domain.N(2);

index{1}=1:Nx:Nx*(Ny-1)+1; % Top boundary of matrix (x(1))
index{2}=Nx:Nx:Nx*Ny; % Bottom boundary of matrix (x(end))
index{3}=1:Nx; % Left boundary of matrix (y(1))
index{4}=Nx*(Ny-1)+1:Nx*Ny; % Right boundary of matrix (y(end))

for i=1:length(domain.BCflag)
    
    if domain.BCflag(i)~=3 % If not Fourier, set BCs
        pde.f(index{i})=BC(:,i);
    end
    
end

% -------------------------------------------------------------------------
% SOLVE HERE
% -------------------------------------------------------------------------

tic
% [v,r]=mg(v0,pde,domain,option);
option.numit=30;
[v,r]=MRR(v0,pde,domain,option);
% [v,r]=bicgstab(v0,pde,domain,option);
toc
disp(rms(r(:)))
