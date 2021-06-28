clear;close all
% -------------------------------------------------------------------------
% Solve Yang KP using Newton's method
% in the form (-u_xxxx+au_xx+bu+cu^2)_xx-u_yy=f

% -------------------------------------------------------------------------
% INPUT PARAMETERS
% -------------------------------------------------------------------------
L(1) = 20;
L(2) = 20;
finestgrid = 8;
coarsestgrid = 5;

% PDE parameters
B=4/3;lambda=-1;
a=@(X,Y) 1/2*(B-1/3);
b=@(X,Y) lambda;
c=@(X,Y) -3/4;
d=@(X,Y) -1/2;

% RHS function
f=@(X,Y) 0*X;

% Initial guess
v0=@(X,Y) 16*lambda*(3-((-2*lambda/(B-1/3))^(1/2)*X).^2+(abs(2*lambda)/(B-1/3)^(1/2)*Y).^2)./(3+((-2*lambda/(B-1/3))^(1/2)*X).^2+(abs(2*lambda)/(B-1/3)^(1/2)*Y).^2).^2;
% v0=@(X,Y) -4*sech(X.^2+Y.^2);
% -------------------------------------------------------------------------
% Set up parameters
% -------------------------------------------------------------------------
N(1) = 2^finestgrid;
N(2) = 2^finestgrid;

% Spectral Wave numbers
k{1} = 2*pi/L(1)*[0:N(1)/2-1 -N(1)/2 -N(1)/2+1:-1]';
k{2} = 2*pi/L(2)*[0:N(2)/2-1 -N(2)/2 -N(2)/2+1:-1]';
[KX,KY] = ndgrid(k{1},k{2});

x{1} = L(1)*(-N(1)/2:N(1)/2-1)'/N(1);
x{2} = L(2)*(-N(2)/2:N(2)/2-1)'/N(2);
[X,Y] = ndgrid(x{1},x{2});

a=a(X,Y);
b=b(X,Y);
c=c(X,Y);
d=d(X,Y);
f=f(X,Y);

v0=v0(X,Y);

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
option.operator=@fourier_Ku_2d;
option.coarsegridsolver=@bicgstab;
option.relaxation=@MRR;
option.restriction=@fourier_restrict_2d_filtered;
option.prolongation=@fourier_prolong_2d_filtered;

% Preconditioner
option.preconditioner=@fourier_Ku_pre_2d;
% Number of precondition relaxations
option.prenumit=1;

% -------------------------------------------------------------------------
% Sort into sctructures
% -------------------------------------------------------------------------
% Assuming constant dx
dx(1) = x{1}(2)-x{1}(1);
dx(2) = x{2}(2)-x{2}(1);

% Sort into structures
domain.L = L;
domain.N = N;
domain.k = k;
domain.dx = dx;

pde.a = a;
pde.b = b;
pde.c = c;
pde.d = d;
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
    jacobian=jacobian_Ku_2d(v,pde,domain);
   
    % Check nonlinear residual
    r=rms(rms(jacobian.f));
    fprintf('Residual Newton = %d\n',r)
    if r<=1e-10
        fprintf('Converged after %d Newton Iterations \n',i-1)
        break
    end
    
    % Solve linear equation
%     option.tol=1e-4*r;
    [e,r]=bicgstab(e0,jacobian,domain,option);
%     [e,r]=mg(e0,jacobian,domain,option);

    % Update correction
    v=v+e;
    
end

if i==20
    
    fprintf('Did not converge to required tolerance after %d Newton Iterations\n',i)
    
end 
toc
