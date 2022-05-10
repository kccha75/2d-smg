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
discretisation=[1 2];

% Boundary conditions for each discretisation

% x(1) a1*u+b1*u'=0 x(end) a2*u+b2*u'=0
alpha1{1}=@(X,Y) 1;
beta1{1}=@(X,Y) 0;
alpha2{1}=@(X,Y) 1;
beta2{1}=@(X,Y) 0; 

% Fourier ... not needed
alpha1{2}=@(X,Y) 1;
beta1{2}=@(X,Y) 0;
alpha2{2}=@(X,Y) 1;
beta2{2}=@(X,Y) 0;

% Boundary condition values (column vector) (Non-fourier only)
% BC=[0 0]; assume they are 0 for now ...

finestgrid = 9;
coarsestgrid = 5;

% PDE Parameters
a=@(X,Y) 1;
b=@(X,Y) 1;
c=@(X,Y) 0;

% RHS
f=@(X,Y) 0*X;

% Exact solution
ue=@(X,Y) (cosh(1/2*(Y-1))-cosh(1)).*exp(sin(X));

% Initial guess
v0=@(X,Y) 0*X;

% -------------------------------------------------------------------------
% Multigrid Options here
% -------------------------------------------------------------------------

% Number of V-cycles if option is chosen, otherwise number of v-cycles done
% after FMG
option.num_vcycles=12;

% Solver / solution tolerance
option.tol=1e-12;

% Relaxations on the up and down cycle during Multigrid
option.Nd=1;
option.Nu=1;

% Multigrid solver options:'V-cycle' or 'FMG'
option.solver='FMG';

% Multigrid scheme: 'Correction' or 'FAS'
option.mgscheme='Correction';

% Operator, coarse grid solver, Relaxation
option.operator=@Lu_2d;
option.coarsegridsolver=@specmatrixsolve_2d;
option.relaxation=@MRR;

% Restriction
option.restriction=@(vf) restrict_2d(vf,@fourier_restrict_filtered,@cheb_restrict);

option.restriction_residual=@(vf) restrict_2d(vf,@fourier_restrict_filtered,@cheb_restrict_residual);

% Prolongation
option.prolongation=@(vc) prolong_2d(vc,@fourier_prolong_filtered,@cheb_prolong);

% Preconditioner
option.preconditioner=@FDmatrixsolve_2d;
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
            N(i) = 2^(finestgrid-2)+1;
            k{i} = (0:N(i)-1)';
            x{i} = cos(pi*k{i}/(N(i)-1));
            dx{i} = x{i}(1:end-1)-x{i}(2:end); % due to x(1)=1, x(end)=-1
            
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
% Set up BCs
% -------------------------------------------------------------------------
BC=cell(4,dim);

for i=1:dim
    
    BC{1,i}=alpha1{i}(X,Y);
    BC{2,i}=beta1{i}(X,Y);
    BC{3,i}=alpha2{i}(X,Y);
    BC{4,i}=beta2{i}(X,Y);
    
end

% -------------------------------------------------------------------------
% Sort into structures
% -------------------------------------------------------------------------
domain.dim = dim;
domain.discretisation = discretisation;
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
index=cell(2,domain.dim); % Left and right, for each dimension

Nx=domain.N(1);
Ny=domain.N(2);

index{1,1}=(1:Nx:Nx*(Ny-1)+1)'; % Top boundary of matrix (x(1))
index{2,1}=(Nx:Nx:Nx*Ny)'; % Bottom boundary of matrix (x(end))

index{1,2}=(1:Nx)'; % Left boundary of matrix (y(1))
index{2,2}=(Nx*(Ny-1)+1:Nx*Ny)'; % Right boundary of matrix (y(end))

for i=1:domain.dim
    
    if domain.discretisation(i)~=1 % If not Fourier, set BCs
        pde.f(index{1,i})=0; % Assume 0 for now ...
        pde.f(index{2,i})=0;
    end
    
end

% -------------------------------------------------------------------------
% Conformal mapping test...
% -------------------------------------------------------------------------

h = @(x) .1*(1-tanh((x-.5)/.1).^2)+.1*(1-tanh((x+.5)/.1).^2); % Bump function

H0=1; % initial height

loops=10;

u=x{1};
x_old=u;
kinv=[0;1./(1i*k{1}(2:N(1)))];

for i=1:loops
    
    y_bc=h(x_old); % assume initially x=u 
    L=H0-trapI(y_bc,dx{1})/(2*pi); % New L value (see boundary integral) extra terms due to Fourier periodic
   
    v=L/2*(x{2}+1); % To set the new domain to be [0 L] (needs fixing)
    
    % Solve laplace's equation, poisson's after transformation
    % need to change coefficients a,b
    % define new RHS, change to homogeneous BCs
    
    [U,V]=ndgrid(u,v);
    Y=y_bc+(H0-y_bc).*V/L;

 
    pde.b=1;
    pde.f=-Lu_2d(Y,pde,domain);

    pde.b=1/(L/2)^2;
    % BCs
    pde.f(index{1,2})=0;
    pde.f(index{2,2})=0;
    
    % Solve here
    [Y,r]=mg(v0,pde,domain,option);
    
    % Transform back to original coordinates
    y=y_bc.*(1-V/L)+H0*V/L+Y;

    % find dy/dv on domain
    dy=2/L*ifct(chebdiff(fct(y'),1));
    dy=dy';
    
    % find new x integrate wrt u

    % look at x on bottom boundary
    xx=dy(:,1);
    x_new=mean(xx)*(x{1}+pi)-pi+real(ifft(kinv.*fft(xx)));
    

    
    % compare old x with new x, break if tol met, loop otherwise
    if norm(x_new - x_old) < 1e-8
        fprintf('Linear residual is %d\n',rms(r(:)))
        fprintf('x diff is %d\n',norm(x_new - x_old))
        fprintf('Reached tolerance after %d iterations!\n',i)
        break;
    end
    
    fprintf('Linear residual is %d\n',rms(r(:)))
    fprintf('x diff is %d\n',norm(x_new - x_old))
    
    x_old=x_new;
    
end
