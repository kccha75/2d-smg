clear;close all
% -------------------------------------------------------------------------
% Solve Yang NLS using Newton's method
% Note: there is a problem if coarsestgrid is 6 or lower 
% 
% -------------------------------------------------------------------------
% INPUT PARAMETERS
% -------------------------------------------------------------------------
L(1) = 10 * pi;
L(2) = 10 * pi;
finestgrid = 8;
coarsestgrid = 7;

% PDE parameters

a=@(X,Y) 1;
b=@(X,Y) 1;

V0=6;mu=-4.11;
% c(x) function
c=@(X,Y) V0*(sin(X).^2+sin(Y).^2)+mu;

% RHS function
f=@(X,Y) 0*X;

% Initial guess
v0=@(X,Y) 0.49*sech(2*sqrt(X.^2+Y.^2));

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
% Multigrid Options here
% -------------------------------------------------------------------------

% Number of V-cycles if option is chosen, otherwise number of v-cycles done
% after FMG
option.num_vcycles=100;

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
option.operator=@fourier_Lu_2d_mid;
option.coarsegridsolver=@cg;
option.relaxation=@MRR;
option.restriction=@fourier_restrict_2d_mid;
option.prolongation=@fourier_prolong_2d_mid;

% Preconditioner
option.preconditioner=@fourier_yang_pre_2d;
% Number of precondition relaxations
option.prenumit=1;

% -------------------------------------------------------------------------
% Sort into sctructures
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
% NEWTON HERE
% -------------------------------------------------------------------------

% New b(x) function in Newton
cnew=c-3*v0.^2;
v=v0;

% Error guess (keep at 0)
e0=zeros(N);

tic
for i=1:20
    
    pde.c=c;
    % Initial RHS of linear equation
    pde.f=f-(option.operator(v,pde,domain)-v.^3);
    
    r=rms(rms(pde.f));
    fprintf('Residual Newton = %d\n',r)
    if r<=1e-10
        fprintf('Converged after %d Newton Iterations \n',i-1)
        break
    end
    
    % Solve linear equation
    pde.c=cnew;

    option.tol=1e-1*r;
%     [e,r]=cg(e0,pde,scheme,option);
    [e,r]=mg(v,pde,domain,option);

    % Update correction
    v=v+e;
    
    cnew=c-3*v.^2;
    
end

if i==20
    
    fprintf('Did not converge to required tolerance after %d Newton Iterations\n',i)
    
end 
toc