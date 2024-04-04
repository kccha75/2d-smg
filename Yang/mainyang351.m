clear;close all;%clc

% -------------------------------------------------------------------------
% Yang paper example 3.4 first part
% coarse grid cannot go below 6, but higher is slightly? better
% -------------------------------------------------------------------------

% x domain [-Dx pi, Dx pi]
Dx=5;

% y domain [-Dy pi, Dy pi]
Dy=5;

% Dimension of problem
dim=2;

% Discretisation flag for each dimension
% 1 - Fourier
% 2 - Cheb
discretisation=[1 1];

% Boundary conditions for each discretisation (if fourier not used)
% x(1) a11*u+b11*u'=rhs11 
alpha1{1}=@(y) 1;
beta1{1}=@(y) 0;
BCRHS1{1}=@(y) 0;

% x(end) a21*u+b21*u'= rhs21
alpha2{1}=@(y) 1;
beta2{1}=@(y) 0; 
BCRHS2{1}=@(y) 0;

% y(1) a12*u+b12*u'=rhs12 
alpha1{2}=@(x) -1;
beta1{2}=@(x) -1;
BCRHS1{2}=@(x) -23121;

% y(end) a22*u+b22*u'= rhs22
alpha2{2}=@(x) -1;
beta2{2}=@(x) -1;
BCRHS2{2}=@(x) -2311;

% Grid size
finestgrid = 9;
coarsestgrid = 7;

% PDE Parameters
aa=@(X,Y) 1/Dx^2;
bb=@(X,Y) 1/Dy^2;

V0=6;mu=5;
cc=@(X,Y) -V0*(sin(X*Dx).^2+sin(Y*Dy).^2)+mu;

% RHS
ff=@(X,Y) 0*X;

% Initial guess
vv0=@(X,Y) 1.15*sech(2*sqrt((X*Dx).^2+(Y*Dy).^2));

% -------------------------------------------------------------------------
% Multigrid Options here
% -------------------------------------------------------------------------

% Number of V-cycles if option is chosen, otherwise number of v-cycles done
% after FMG or number of decent iterations
option.numit=1;

% Solver / solution tolerance
option.tol=1e-10;

% Relaxations on the up and down cycle during Multigrid
option.Nd=1;
option.Nu=1;

% Multigrid solver options:'V-cycle' or 'FMG'
option.solver='FMG';

% Multigrid scheme: 'Correction' or 'FAS'
option.mgscheme='Correction';

% Operator, coarse grid solver, Relaxation
option.operator=@Lu;
option.coarsegridsolver=@cg;
option.relaxation=@MRR;

% Restriction for pde coefficients
option.restriction=@(vf) restrict_2d(vf,@fourier_restrict_filtered,@fourier_restrict_filtered);

% Restriction for residual and RHS
option.restriction_residual=@(vf) restrict_2d(vf,@fourier_restrict_filtered,@fourier_restrict_filtered);

% Prolongation
option.prolongation=@(vc) prolong_2d(vc,@fourier_prolong_filtered,@fourier_prolong_filtered);

% Preconditioner
option.preconditioner=@yang_NLS_pre;

% Number of preconditioned relaxations
option.prenumit=1;

% -------------------------------------------------------------------------
% Set up parameters
% -------------------------------------------------------------------------
m=3;
t_mg=zeros(1,m);
t_cg=zeros(1,m);
mg_tol=zeros(1,m);

for jj=1:m % loop grid sizes

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
            dx{i} = x{i}(2:end)-x{i}(1:end-1); % is negative but ok :)
            
    end
    
end

[X,Y] = ndgrid(x{1},x{2});

a=aa(X,Y);
b=bb(X,Y);
c=cc(X,Y);
f=ff(X,Y);

v0=vv0(X,Y);

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
C=c;
% -------------------------------------------------------------------------
% NEWTON HERE
% -------------------------------------------------------------------------

% New b(x) function in Newton
cnew=C-3*v0.^2;
v=v0;

% Error guess (keep at 0)
e0=zeros(Nx,Ny);

% -------------------------------------------------------------------------
% SMG
% -------------------------------------------------------------------------
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

    option.tol=1e-10;
    [e,r]=mg(v,pde,domain,option);
    mg_tol(jj)=rms(r(:));

    % Update correction
    v=v+real(e);
    
    cnew=c-3*v.^2;
    
end

if i==20
    
    fprintf('Did not converge to required tolerance after %d Newton Iterations\n',i)
    
end 
t_mg(jj)=toc;

% -------------------------------------------------------------------------
% CG
% -------------------------------------------------------------------------

% New b(x) function in Newton
cnew=C-3*v0.^2;
v=v0;

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

    option.tol=max(1e-10,mg_tol(jj));
    [e,r]=cg(e0,pde,domain,option);

    % Update correction
    v=v+real(e);
    
    cnew=c-3*v.^2;
    
end

if i==20
    
    fprintf('Did not converge to required tolerance after %d Newton Iterations\n',i)
    
end 
t_cg(jj)=toc;

finestgrid = finestgrid+1;
end

% -------------------------------------------------------------------------
% Plot
% -------------------------------------------------------------------------
figure('Position',[300 300 600 300]); fsz=15; lw=2;

subplot(1,2,1)
surf(X*Dx,Y*Dy,v,'LineStyle','none')
xlabel('$x$','interpreter','latex','fontsize',fsz)
ylabel('$y$','interpreter','latex','fontsize',fsz)
zlabel('$u$','interpreter','latex','fontsize',fsz)

subplot(1,2,2)
M=linspace(1,m,m)+coarsestgrid;
semilogy(M,t_cg,'-x',M,t_mg,'-o')
xlabel('$2^N$','interpreter','latex','fontsize',fsz)
ylabel('$t$','interpreter','latex','fontsize',fsz)