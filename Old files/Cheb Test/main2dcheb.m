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

L(1) = 2;
L(2) = 2*pi;
finestgrid = 6;
coarsestgrid = 3;

% PDE Parameters
a=@(X,Y) 1;
b=@(X,Y) 1;
c=@(X,Y) 1;

% RHS
% f=@(X,Y) exp(sin(Y)).*(cosh(X)+exp(X+Y).*(-cosh(1)+cosh(X))+(cosh(1)-cosh(X)).*(-cos(Y).^2+sin(Y)));
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
option.num_vcycles=2;

% Solver / solution tolerance
option.tol=1e-12;

% Relaxations on the up and down cycle during Multigrid
option.Nd=1;
option.Nu=1;

% Multigrid solver options:'V-cycle' or 'FMG'
option.solver='V-cycle';

% Multigrid scheme: 'Correction' or 'FAS'
option.mgscheme='FAS';

% Operator, coarse grid solver, Relaxation, Restriction, Prolongation options
option.operator=@cheb_fourier_Lu_2d;
option.coarsegridsolver=@cheb_fourier_matrixsolve;
option.relaxation=@MRR;

% Restriction
x_restrict=@cheb_restrict;
y_restrict=@fourier_restrict_filtered;
option.restriction=@(vf) restrict_2d(vf,x_restrict,y_restrict);

% Prolongation
x_prolong=@cheb_prolong;
y_prolong=@fourier_prolong_filtered;
option.prolongation=@(vc) prolong_2d(vc,x_prolong,y_prolong);

% Preconditioner
option.preconditioner=@cheb_fourier_FD_2d;
% Number of precondition relaxations
option.prenumit=1;
% -------------------------------------------------------------------------
% Set up parameters
% -------------------------------------------------------------------------
N(1) = 2^finestgrid+1; % cheb +1 points!s
N(2) = 2^finestgrid;

% Spectral Wave numbers
k{1} = (0:N(1)-1)'; % Cheb
k{2} = 2*pi/L(2)*[0:N(2)/2-1 -N(2)/2 -N(2)/2+1:-1]'; % Fourier

x{1} = cos(pi*k{1}/(N(1)-1)); % Cheb
x{2} = L(2)*(-N(2)/2:N(2)/2-1)'/N(2); % Fourier

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
% dx
dx{1} = x{1}(2:end)-x{1}(1:end-1);
dx{2} = x{2}(2)-x{2}(1);

% Sort into structures
domain.L = L;
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
% option.numit=10;


[pde,domain]=setstructures_test(v0,pde,domain,option);

tic
[v,r]=mg2(v0,pde,domain,option);
% [v,r]=MRR(v0,pde,domain,option);
% [v,r]=bicgstab(v0,pde,domain,option);
toc

% surf(pde.f-option.operator(ue,pde,domain))

% v=cheb_fourier_FD_2d(1,pde,domain,1);
% r=pde.f-option.operator(v,pde,domain);
% disp(rms(r(:)))
% figure;surf(v-ue)