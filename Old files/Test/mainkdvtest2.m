clear;close all;%clc
% -------------------------------------------------------------------------
% Solve KP using Newton's method in the form -au_xx+cu+3du^2+b*intint(u_yy)dxdx=f
% Fourier collocation points
% Kind of works??? does converge but add any Y component and it doesn't
% like it ...
% -------------------------------------------------------------------------
% INPUT PARAMETERS
% -------------------------------------------------------------------------
L(1) = 100;
L(2) = 100;
finestgrid = 9;
coarsestgrid = 6;

% PDE Parameters

a=@(X,Y) 1;%+epsilon*exp(cos(X+Y));
b=@(X,Y) -1/2;%+epsilon*exp(cos(X+Y));
c=@(X,Y) 0;

% RHS function
f=@(X,Y) -8*(sech(X).^2);

% Initial guess
v0=@(X,Y) 2.5*sech(X).^2;

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
option.solver='FMG';

% Multigrid scheme: 'Correction' or 'FAS'
option.mgscheme='Correction';

% Operator, coarse grid solver, Relaxation, Restriction, Prolongation options
option.operator=@fourier_KPu_2d;
option.coarsegridsolver=@cg;
option.relaxation=@MRR;
option.restriction=@fourier_restrict_2d_filtered;
option.prolongation=@fourier_prolong_2d_filtered;

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
option.grids=finestgrid-coarsestgrid+1;

% -------------------------------------------------------------------------
% NEWTON HERE
% -------------------------------------------------------------------------

% New b(x) function in Newton
d=-3;
cnew=c+2*d*v0;
% cnew=c+3*alias_2d_test(v0,v0);
v=v0;

% Error guess (keep at 0)
e0=zeros(N);

tic
for i=1:20
    
    pde.c=c;
    % Initial RHS of linear equation
    pde.f=f-(option.operator(v,pde,domain)+d.*v.^2);
   
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
%     e=fourier_matrixsolve_mid(e0,pde,domain,option);
%     [e,r]=mg(e0,pde,domain,option);

    % Update correction
    v=v+e;
    
    cnew=c+2*d.*v;
%     cnew=c+3*alias_2d_test(v,v);
end

if i==20
    
    fprintf('Did not converge to required tolerance after %d Newton Iterations\n',i)
    
end 
toc