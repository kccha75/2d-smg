clear;close all;%clc
% -------------------------------------------------------------------------
% Solve PDE -u_xxxx + au_xx + bu + cu^2 = f using Fourier Spectral Multigrid at
% Fourier collocation points and Newton iterations
% Used on 5th order KdV equation
% -------------------------------------------------------------------------
% INPUT PARAMETERS
% -------------------------------------------------------------------------
L(1) = 30*pi;
finestgrid = 9;
coarsestgrid = 6;

% PDE Parameters
mu=-1.2;
a=@(X) -2;
b=@(X) mu;
c=@(X) -3;

f=@(X) 0*X;

% Exact solution
% ue=@(x) sin(x).^2;

% Initial guess
v0=@(X) 0.25*sech(0.3*X).*cos(X);

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
option.mgscheme='Correction';

% Operator, coarse grid solver, Relaxation, Restriction, Prolongation options
option.operator=@fourier_Lu_1d_kdv5;
option.coarsegridsolver=@cg;
option.relaxation=@MRR;
option.restriction=@fourier_restrict_filtered;
option.prolongation=@fourier_prolong_filtered;

% Preconditioner
option.preconditioner=@fourier_yang_pre_1d_kdv5;
% Number of precondition relaxations
option.prenumit=1;

% -------------------------------------------------------------------------
% Set up parameters
% -------------------------------------------------------------------------
N(1) = 2^finestgrid;
N(2) = 1;

% Spectral Wave numbers
k{1} = 2*pi/L(1)*[0:N(1)/2-1 -N(1)/2 -N(1)/2+1:-1]';

x{1} = L(1)*(-N(1)/2:N(1)/2-1)'/N(1);
X=ndgrid(x{1});

a=a(X);
b=b(X);
c=c(X);
f=f(X);

% ue=ue(x);
v0=v0(X);

% -------------------------------------------------------------------------
% Sort into structures
% -------------------------------------------------------------------------
% Assuming constant dx
dx(1) = x{1}(2)-x{1}(1);

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

v=v0;

% Error guess (keep at 0)
e0=zeros(N);

tic
for i=1:20
    
    % Calculate Jacobian for linear equation
    jacobian=jacobian_kdv5_1d(v,pde,domain);
   
    r=rms(rms(jacobian.f));
    fprintf('Residual Newton = %d\n',r)
    if r<=1e-10
        fprintf('Converged after %d Newton Iterations \n',i-1)
        break
    end
    
    % Solve linear equation
%     option.tol=1e-2*r;
%     [e,r]=cg(e0,pde,domain,option);
    [e,r]=mg(e0,jacobian,domain,option);

    % Update correction
    v=v+e;

end

if i==20
    
    fprintf('Did not converge to required tolerance after %d Newton Iterations\n',i)
    
end 
toc