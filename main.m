clear;close all;%clc
% -------------------------------------------------------------------------
% Solve PDE -(au_x)_x-(bu_y)_y+cu = f using Fourier Spectral Multigrid at
% Fourier collocation points
%
% -------------------------------------------------------------------------
% INPUT PARAMETERS
% -------------------------------------------------------------------------
L(1) = 2 * pi;
L(2) = 2 * pi;
finestgrid = 7;
coarsestgrid = 3;

% PDE Parameters

epsilon=0.0;
a=@(X,Y) 1;%+epsilon*exp(cos(X+Y));
b=@(X,Y) 1;%+epsilon*exp(cos(X+Y));
c=@(X,Y) 0*X;

f=@(X,Y) -pi*(-1*(1+exp(cos(X+Y))*epsilon).*cos(Y).*cos(pi/4+pi*cos(Y)).*sin(pi/4+pi*cos(X))+ ...
    exp(cos(X+Y))*epsilon.*cos(pi/4+pi*cos(Y)).*sin(Y).*sin(X+Y).*sin(pi/4+pi*cos(X))-...
    ((1+exp(cos(X+Y))*epsilon).*cos(X).*cos(pi/4+pi*cos(X))-exp(cos(X+Y))*epsilon.*cos(pi/4+pi*cos(X)).*sin(X).*sin(X+Y)+...
    pi*(1+exp(cos(X+Y))*epsilon).*(sin(X).^2+sin(Y).^2).*sin(pi/4+pi*cos(X))).*sin(pi/4+pi*cos(Y)));

% Exact solution
ue=@(X,Y) sin(1*pi*cos(X)+pi/4).*sin(1*pi*cos(Y)+pi/4)-1/2*besselj(0,pi)^2;

% Initial guess
Nx=2^finestgrid;Ny=2^finestgrid;
v0=@(X,Y) ue(X,Y)+2*(rand(Nx,Ny)-0.5);

% -------------------------------------------------------------------------
% Multigrid Options here
% -------------------------------------------------------------------------

% Number of V-cycles if option is chosen, otherwise number of v-cycles done
% after FMG
option.num_vcycles=10;

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
option.operator=@fourier_Lu_2d_mid;
option.coarsegridsolver=@fourier_matrixsolve_2d;
option.relaxation=@MRR;
option.restriction=@fourier_restrict_2d_mid;
option.prolongation=@fourier_prolong_2d_mid;

% Preconditioner
option.preconditioner=@fourier_RBGSLineRelax_2d;
% Number of precondition relaxations
option.prenumit=0;

% -------------------------------------------------------------------------
% Set up parameters
% -------------------------------------------------------------------------
N(1) = 2^finestgrid;
N(2) = 2^finestgrid;

% Spectral Wave numbers
k(:,1) = 2*pi/L(1)*[0:N(1)/2-1 -N(1)/2 -N(1)/2+1:-1]';
k(:,2) = 2*pi/L(2)*[0:N(2)/2-1 -N(2)/2 -N(2)/2+1:-1]';
[KX,KY] = ndgrid(k(:,1),k(:,2));

x(:,1) = L(1)*(-N(1)/2:N(1)/2-1)'/N(1);
x(:,2) = L(2)*(-N(2)/2:N(2)/2-1)'/N(2);
[X,Y] = ndgrid(x(:,1),x(:,2));

a=a(X,Y);
b=b(X,Y);
c=c(X,Y);
f=f(X,Y);

ue=ue(X,Y);
v0=v0(X,Y);

% -------------------------------------------------------------------------
% Sort into structures
% -------------------------------------------------------------------------
% Assuming constant dx
dx(1) = x(2,1)-x(1,1);
dx(2) = x(2,2)-x(1,2);

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

% -------------------------------------------------------------------------
% Solve
% -------------------------------------------------------------------------
% MG here

tic
[v,r]=mg(v0,pde,domain,option);

% Check if Poisson type problem, then scale for mean 0 solution
if max(abs(pde.c(:)))<1e-12
    v=v-1/(Nx(1)*N(2))*sum(sum(v));
end

toc

% surf(real(fft2(ue-v)))

% option.numit=3;
% [v,r]=MRR(v0,pde,scheme,option);
% surf(real(fft2(v-ue)))

% Pcg comparison
% option.preconditioner=@fourier_yang_pre;
% option.prenumit=0;
% tic
% [vv,rr]=bicgstab1(v0,pde,scheme,option);
% vv=vv-1/(Nx*Ny)*sum(sum(vv));
% toc
% 
% % MRR comparison
% option.numit=20;
% tic
% [vv,rr]=MRR(v0,pde,scheme,option);
% vv=vv-1/(Nx*Ny)*sum(sum(vv));
% toc
% disp(rms(rms(rr)))