clear;close all;%clc
% -------------------------------------------------------------------------
% Solve PDE (-u_xx + au_x + bu + cu^2)_xx = f_xx using Fourier Spectral Multigrid at
% Fourier collocation points and Newton iterations
% Used on fKdV equation (Ee and Clarke form)
% -------------------------------------------------------------------------
% INPUT PARAMETERS
% -------------------------------------------------------------------------
L(1) = 50;
finestgrid = 9;
coarsestgrid = 5;

% PDE Parameters

a=@(x) 0*x;
b=@(x) 0;
c=@(x) -3;
gamma=8;
f=@(x) -16*(-2+cosh(2*x)).*sech(x).^4;

% Exact solution
% ue=@(x) sin(x).^2;

% Initial guess
v0=@(x) 1.999999999*sech(x).^2;

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
option.operator=@fourier_Ku_1d;
option.coarsegridsolver=@cg;
option.relaxation=@MRR;
option.restriction=@fourier_restrict_filtered;
option.prolongation=@fourier_prolong_filtered;

% Preconditioner
option.preconditioner=@fourier_FD_1d;
% Number of precondition relaxations
option.prenumit=0;

% -------------------------------------------------------------------------
% Set up parameters
% -------------------------------------------------------------------------
N(1) = 2^finestgrid;
N(2) = 1;

% Spectral Wave numbers
k(:,1) = 2*pi/L(1)*[0:N(1)/2-1 -N(1)/2 -N(1)/2+1:-1]';

x(:,1) = L(1)*(-N(1)/2:N(1)/2-1)'/N(1);

a=a(x);
b=b(x);
c=c(x);
f=f(x);

% ue=ue(x);
v0=v0(x);

% -------------------------------------------------------------------------
% Sort into structures
% -------------------------------------------------------------------------
% Assuming constant dx
dx(1) = x(2,1)-x(1,1);

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
c=pde.c;
bnew=b+2*c.*v0;
v=v0;

% Error guess (keep at 0)
e0=zeros(N);

tic
for i=1:20
    
    pde.b=b;
    % Initial RHS of linear equation
    pde.c=c;
    pde.f=f-(option.operator(v,pde,domain));
   
    r=rms(rms(pde.f));
    fprintf('Residual Newton = %d\n',r)
    if r<=1e-10
        fprintf('Converged after %d Newton Iterations \n',i-1)
        break
    end
    
    % Solve linear equation
    pde.b=bnew;
    pde.c=0*x;

%     option.tol=1e-2*r;
    [e,r]=cg(e0,pde,domain,option);
%     e=fourier_matrixsolve_1d(e0,pde,domain,option);
%     [e,r]=mg(e0,pde,domain,option);

    % Update correction
    v=v+e;
    
    bnew=b+2*c.*v;

end

if i==20
    
    fprintf('Did not converge to required tolerance after %d Newton Iterations\n',i)
    
end 
toc