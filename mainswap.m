clear;close all;%clc
%
% Not so simple main for Cheb Fourier SMG, uses nonhomogenous BCs with
% swapped directions for testing
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

% Boundary conditions for each discretisation (if fourier not used)
% x(1) a11*u+b11*u'=rhs11 
alpha1{1}=@(y) 123123;
beta1{1}=@(y) 123124;
BCRHS1{1}=@(y) 1125125;

% x(end) a21*u+b21*u'= rhs21
alpha2{1}=@(y) 11231;
beta2{1}=@(y) 124124; 
BCRHS2{1}=@(y) 235235;

% y(1) a12*u+b12*u'=rhs12 
alpha1{2}=@(x) 0;
beta1{2}=@(x) 1;
BCRHS1{2}=@(x) 1+1/2*exp(sin(x))*sinh(1);

% y(end) a22*u+b22*u'= rhs22
alpha2{2}=@(x) 1;
beta2{2}=@(x) 0; 
BCRHS2{2}=@(x) -1+exp(sin(x)).*(1-cosh(1));

% Grid size
finestgrid = 6;
coarsestgrid = 3;

% PDE Parameters
a=@(X,Y) 1;
b=@(X,Y) 1;
c=@(X,Y) 1;

% RHS
f=@(X,Y) Y+exp(sin(X)).*cosh(1).*(-1-cos(X).^2+sin(X))-1/4*exp(sin(X)).*cosh((1+Y)/2).*(-5-4*cos(X).^2+4*sin(X));

% Exact solution
ue=@(X,Y) (cosh(1/2*(-Y-1))-cosh(1)).*exp(sin(X))+Y;

% Initial guess
v0=@(X,Y) rand(size(X));

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
option.solver='FMG';

% Multigrid scheme: 'Correction' or 'FAS'
option.mgscheme='Correction';

% Operator, coarse grid solver, Relaxation
option.operator=@Lu_2d;
option.coarsegridsolver=@specmatrixsolve_2d;
option.relaxation=@MRR;

% Restriction for pde coefficients
option.restriction=@(vf) restrict_2d(vf,@fourier_restrict_filtered,@cheb_restrict);

% Restriction for residual and RHS
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
k=cell(1,dim);
dx=cell(1,dim);

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
BCRHS=cell(2,dim);

y=fliplr(x);

for i=1:dim
    
    BC{1,i}=alpha1{i}(x{i});
    BC{2,i}=beta1{i}(x{i});
    BC{3,i}=alpha2{i}(x{i});
    BC{4,i}=beta2{i}(x{i});

    BCRHS{1,i}=BCRHS1{i}(y{i});
    BCRHS{2,i}=BCRHS2{i}(y{i});
    
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

index{1,1}=1:Nx:Nx*(Ny-1)+1; % Top boundary of matrix (x(1))
index{2,1}=Nx:Nx:Nx*Ny; % Bottom boundary of matrix (x(end))

index{1,2}=1:Nx; % Left boundary of matrix (y(1))
index{2,2}=Nx*(Ny-1)+1:Nx*Ny; % Right boundary of matrix (y(end))

for i=1:domain.dim
    
    if domain.discretisation(i)~=1 % If not Fourier, set BCs
        pde.f(index{1,i})=BCRHS{1,i}; % x(1) boundary
        pde.f(index{2,i})=BCRHS{2,i}; % x(end) boundary
    end
    
end


% -------------------------------------------------------------------------
% SOLVE HERE
% -------------------------------------------------------------------------

tic
[v,r]=mg(v0,pde,domain,option);
% option.numit=30;
% [v,r]=MRR(v0,pde,domain,option);
% [v,r]=bicgstab(v0,pde,domain,option);
toc
disp(rms(r(:)))
tic
option.numit=5;
[vv,rr]=MRR(v0,pde,domain,option);
disp(rms(r(:)))
toc

surf(X,Y,v);xlabel('x');ylabel('y');title('Numerical solution of Poissons equation')
figure;contour(X,Y,v);xlabel('x');ylabel('y')
figure;surf(X,Y,abs(v-ue));xlabel('x');ylabel('y');title('Error compared to exact solution')

