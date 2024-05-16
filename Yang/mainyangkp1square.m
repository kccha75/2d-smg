clear;close all;%clc

% -------------------------------------------------------------------------
% Yang paper example 5.1 solves +ve case with various coarse grids
% Uses (N-1)*N grid size (grid anisotrophy)
% Very very rough code ... is not clean
% +ve initial condition can use only k=2, finestgrid=10 coarsest grid=9
% -------------------------------------------------------------------------

% Dimension of problem
dim=2;

% Domain size
Lx=120*pi;
Ly=60*pi;

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
finestgrid = 10;
coarsestgrid = 9;

% PDE Parameters
aa=@(X,Y) 1;
bb=@(X,Y) 2;
mu=-1.2;
cc=@(X,Y) -mu;
dd=@(X,Y) 1;

% RHS
ff=@(X,Y) 0*X;

% Initial guess
vv0=@(X,Y) 0.43*sech(0.3*sqrt((X).^2+(Y).^2)).*cos(X);

% -------------------------------------------------------------------------
% Multigrid Options here
% -------------------------------------------------------------------------

% Number of V-cycles if option is chosen, otherwise number of v-cycles done
% after FMG or number of decent iterations
option.numit=0;

% Solver / solution tolerance
% option.tol=1e-5;

% Relaxations on the up and down cycle during Multigrid
option.Nd=1;
option.Nu=1;

% Multigrid solver options:'V-cycle' or 'FMG'
option.solver='FMG';

% Multigrid scheme: 'Correction' or 'FAS'
option.mgscheme='Correction';

% Operator, coarse grid solver, Relaxation
option.operator=@Lu_kp;
option.coarsegridsolver=@bicgstabmg;
option.relaxation=@MRR;

% Restriction for pde coefficients
option.restriction=@(vf) restrict_2d(vf,@fourier_restrict_filtered,@fourier_restrict_filtered);

% Restriction for residual and RHS
option.restriction_residual=@(vf) restrict_2d(vf,@fourier_restrict_filtered,@fourier_restrict_filtered);

% Prolongation
option.prolongation=@(vc) prolong_2d(vc,@fourier_prolong_filtered,@fourier_prolong_filtered);

% Preconditioner
option.preconditioner=@yang_kp_pre;

% Number of preconditioned relaxations
option.prenumit=1;

% -------------------------------------------------------------------------
% Set up parameters
% -------------------------------------------------------------------------
m=4;
t_mg=zeros(1,m);
t_cg=zeros(1,m);
mg_tol=zeros(1,m);
finestgrid=finestgrid+m-1;
for jj=m:-1:1 % loop grid sizes
if finestgrid==14
    option.tol=1e-5; % fine grid N=14 OR finestgrid+m-1
elseif finestgrid==13
    option.tol=1e-6; % N=13
elseif finestgrid==12
    option.tol=1e-8; % N=12
elseif finestgrid==11
    option.tol=1e-9; % N=11
elseif finestgrid==10
    option.tol=1e-10; % N=10
end
tol=option.tol;
N=zeros(1,dim);
x=cell(1,dim);
k=cell(1,dim);
dx=cell(1,dim);

N(1)=2^(finestgrid);
k{1}=2*pi/Lx*[0:N(1)/2-1 -N(1)/2 -N(1)/2+1:-1]';
x{1}=Lx*(-N(1)/2:N(1)/2-1)'/N(1);

N(2)=2^(finestgrid-1);
k{2}=2*pi/Ly*[0:N(2)/2-1 -N(2)/2 -N(2)/2+1:-1]';
x{2}=Ly*(-N(2)/2:N(2)/2-1)'/N(2);

[X,Y] = ndgrid(x{1},x{2});

a=aa(X,Y);
b=bb(X,Y);
c=cc(X,Y);
d=dd(X,Y);
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
pde.d = d;
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
cnew=C+6*v0;
v=v0;

% Error guess (keep at 0)
e0=zeros(Nx,Ny);

% -------------------------------------------------------------------------
% SMG
% -------------------------------------------------------------------------
for kk=2:-1:1
    option.coarsestgrid=coarsestgrid-kk+1;
    option.grids=option.finestgrid-option.coarsestgrid+1;
    cnew=C+6*v0;
    v=v0;
tic
for i=1:20
    
    pde.c=c;
    % Initial RHS of linear equation
    pde.f=f-(option.operator(v,pde,domain)+3*real(ifft(-k{1}.^2.*fft(v.^2))));
    
    r=rms(rms(pde.f));
    fprintf('Residual Newton = %d\n',r)
    if r<=option.tol
        fprintf('Converged after %d Newton Iterations \n',i-1)
        break
    end
    
    % Solve linear equation
    pde.c=cnew;

    option.reltol=tol/r;
    [e,r]=mg(v,pde,domain,option);
    mg_tol(jj)=rms(r(:));

    % Update correction
    v=v+real(e);

    % Mean 0 solution
    v=fft2(v);
    v(1,1)=0;
    v=real(ifft2(v));
    
    cnew=c+6*v;
    
end

if i==20
    
    fprintf('Did not converge to required tolerance after %d Newton Iterations\n',i)
    
end 
t_mg(jj,kk)=toc;
end

% -------------------------------------------------------------------------
% CG
% -------------------------------------------------------------------------

% New b(x) function in Newton
cnew=C+6*v0;
v=v0;

tic
for i=1:20
    
    pde.c=c;
    % Initial RHS of linear equation
    pde.f=f-(option.operator(v,pde,domain)+3*real(ifft(-k{1}.^2.*fft(v.^2))));
    
    r=rms(rms(pde.f));
    fprintf('Residual Newton = %d\n',r)
    if r<=option.tol
        fprintf('Converged after %d Newton Iterations \n',i-1)
        break
    end
    
    % Solve linear equation
    pde.c=cnew;

    option.reltol=max(tol,mg_tol(jj))/r;
    [e,r]=option.coarsegridsolver(e0,pde,domain,option);

    % Update correction
    v=v+real(e);

    % Mean 0 solution
    v=fft2(v);
    v(1,1)=0;
    v=real(ifft2(v));
    
    cnew=c+6*v;    
    
end

if i==20
    
    fprintf('Did not converge to required tolerance after %d Newton Iterations\n',i)
    
end 
t_cg(jj)=toc;

finestgrid = finestgrid-1;
end

% -------------------------------------------------------------------------
% Plot
% -------------------------------------------------------------------------
figure('Position',[300 300 600 300]); fsz=15; lw=2;

subplot(1,2,1)
surf(X,Y,v,'LineStyle','none')
xlabel('$x$','interpreter','latex','fontsize',fsz)
ylabel('$y$','interpreter','latex','fontsize',fsz)
zlabel('$u$','interpreter','latex','fontsize',fsz)
axis([-25,25,-25,25])

subplot(1,2,2)
M=linspace(1,m,m)+coarsestgrid;
semilogy(M,t_cg,'-x',M,t_mg(:,2),'-o',M,t_mg(:,1),'-o')
xlabel('$N$','interpreter','latex','fontsize',fsz)
ylabel('$t$','interpreter','latex','fontsize',fsz)
legend('CG','SMG coarse grid N=8','SMG coarse grid N=9','Location','SouthEast')