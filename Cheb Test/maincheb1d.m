clear;close all;%clc
% -------------------------------------------------------------------------
% Solve PDE u_xx + au_x + bu = f using Cheb Spectral Multigrid at
% Cheb collocation points 
% -------------------------------------------------------------------------
% INPUT PARAMETERS
% -------------------------------------------------------------------------
L(1) = 30*pi;
finestgrid = 3;
coarsestgrid = 3;

% PDE Parameters
a=@(X) 0;
b=@(X) exp(X);

f=@(X) exp(X).*(X-1).*(X+1)+2;

% Exact solution
% ue=@(x) sin(x).^2;

% Initial guess
v0=@(X) 0*X;

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
option.solver='V-cycle';

% Multigrid scheme: 'Correction' or 'FAS'
option.mgscheme='Correction';

% Operator, coarse grid solver, Relaxation, Restriction, Prolongation options
option.operator=@cheb_Lu_1d;
option.coarsegridsolver=@bicgstab;
option.relaxation=@MRR;
option.restriction=@cheb_restrict;
option.prolongation=@cheb_prolong;

% Preconditioner
option.preconditioner=@fd;
% Number of precondition relaxations
option.prenumit=0;

% -------------------------------------------------------------------------
% Set up parameters
% -------------------------------------------------------------------------
N(1) = 2^finestgrid;
N(2) = 1;

% Spectral Wave numbers
k{1} = (0:N(1)-1)';

x{1} = cos(pi*k{1}/(N(1)-1));
X=ndgrid(x{1});

a=a(X);
b=b(X);
f=f(X);

% ue=ue(x);
v0=v0(X);

% -------------------------------------------------------------------------
% Sort into structures
% -------------------------------------------------------------------------
% Assuming constant dx
dx = x{1}(2:end)-x{1}(1:end-1);

% Sort into structures
domain.L = L;
domain.N = N;
domain.k = k;
domain.dx = dx;

pde.a = a;
pde.b = b;
pde.f = f;

option.finestgrid=finestgrid;
option.coarsestgrid=coarsestgrid;
option.grids=finestgrid-coarsestgrid+1;

% -------------------------------------------------------------------------
% SOLVE HERE
% -------------------------------------------------------------------------
tic

% [v,r]=mg(v0,pde,domain,option);
[v,r]=bicgstab(v0,pde,domain,option);

toc