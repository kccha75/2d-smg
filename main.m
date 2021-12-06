clear;close all;%clc
% -------------------------------------------------------------------------
% Solve PDE au_xx + bu_yy + cu = f using Fourier Cheb Spectral Multigrid
% -------------------------------------------------------------------------
% INPUT PARAMETERS
% -------------------------------------------------------------------------

% Discretisation flag for each direction
% 1 - Fourier
% 2 - Cheb
discretisation=[2 1];

finestgrid = 4;
coarsestgrid = 3;

% PDE Parameters
a=@(X,Y) 1;
b=@(X,Y) 1;
c=@(X,Y) 1;

% RHS
f=@(X,Y) exp(sin(Y)).*(cosh(X).*(2+cos(Y).^2-sin(Y))+cosh(1).*(-1-cos(Y).^2+sin(Y)));

% Exact solution
ue=@(X,Y) (cosh(X)-cosh(1)).*exp(sin(Y));

% Initial guess
v0=@(X,Y) 0*X;

% -------------------------------------------------------------------------
% Multigrid Options here
% -------------------------------------------------------------------------

% Number of V-cycles if option is chosen, otherwise number of v-cycles done
% after FMG
option.num_vcycles=0;

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
option.operator=@cheb_fourier_Lu_2d;
option.coarsegridsolver=@cheb_fourier_matrixsolve;
option.relaxation=@MRR;

% Restriction
x_restrict=@cheb_restrict_dirichlet;
y_restrict=@fourier_restrict_filtered;
option.restriction=@(vf) restrict_2d(vf,x_restrict,y_restrict);

% Prolongation
x_prolong=@cheb_prolong_dirichlet;
y_prolong=@fourier_prolong_filtered;
option.prolongation=@(vc) prolong_2d(vc,x_prolong,y_prolong);

% Preconditioner
option.preconditioner=@cheb_fourier_FD_2d;
% Number of preconditioned relaxations
option.prenumit=1;

% -------------------------------------------------------------------------
% Set up parameters
% -------------------------------------------------------------------------
d=length(discretisation); % Dimension of problem
N=zeros(1,d);
x=cell(size(discretisation));
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

option.discretisation = discretisation;
% -------------------------------------------------------------------------
% SOLVE HERE
% -------------------------------------------------------------------------
pde.f(1,:)=0;pde.f(end,:)=0;

tic
[v,r]=mg(v0,pde,domain,option);
% option.numit=10;
% [v,r]=MRR(v0,pde,domain,option);
% [v,r]=bicgstab(v0,pde,domain,option);
toc
