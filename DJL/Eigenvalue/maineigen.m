clear;close all;%clc
% -------------------------------------------------------------------------
% Solve PDE au_xx + bu_yy + cu = f using Fourier Cheb Spectral Multigrid
% -------------------------------------------------------------------------
% INPUT PARAMETERS
% -------------------------------------------------------------------------

% Dimension of problem
dim=1;

% Discretisation flag for each dimension
% 1 - Fourier
% 2 - Cheb
discretisation=[2];

% Boundary conditions for each discretisation

% x(1) a1*u+b1*u'=0 x(end) a2*u+b2*u'=0
alpha1{1}=@(X,Y) 1;
beta1{1}=@(X,Y) 0;
alpha2{1}=@(X,Y) 1;
beta2{1}=@(X,Y) 0; 

% Boundary condition values (column vector) (Non-fourier only)
% BC=[0 0]; assume they are 0 for now ...

finestgrid = 6;
coarsestgrid = 3;

% PDE Parameters
a=@(X) 1;
b=@(X) 1;
c=@(X) 1;

% RHS
f=@(X) 0*X;

% Exact solution
% ue=@(X,Y) (cosh(1/2*(X-1))-cosh(1)).*exp(sin(Y));

% Initial guess
v0=@(X) 0*X;

% -------------------------------------------------------------------------
% Multigrid Options here
% -------------------------------------------------------------------------

% Number of V-cycles if option is chosen, otherwise number of v-cycles done
% after FMG
option.num_vcycles=5;

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
option.operator=@Lu_1d;
option.coarsegridsolver=@specmatrixsolve_1d;
option.relaxation=@MRR;

% Restriction
option.restriction=@cheb_restrict;

% Prolongation
option.prolongation=@cheb_prolong;

% Preconditioner
option.preconditioner=@cheb_FD_1d;
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
            N(i) = 2^finestgrid+1;
            k{i} = (0:N(i)-1)';
            x{i} = cos(pi*k{i}/(N(i)-1));
            dx{i} = x{i}(1:end-1)-x{i}(2:end); % due to x(1)=1, x(end)=-1
            
    end
    
end

X = ndgrid(x{1});

a=a(X);
b=b(X);
c=c(X);
f=f(X);

% ue=ue(X,Y);
v0=v0(X);

% -------------------------------------------------------------------------
% Set up BCs
% -------------------------------------------------------------------------
BC=cell(4,dim);

for i=1:dim
    
    BC{1,i}=alpha1{i}(X);
    BC{2,i}=beta1{i}(X);
    BC{3,i}=alpha2{i}(X);
    BC{4,i}=beta2{i}(X);
    
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

if domain.discretisation~=1 % If not Fourier, set BCs
	pde.f(1)=0; % Assume 0 for now ...
	pde.f(end)=0;
end


% -------------------------------------------------------------------------
% SOLVE HERE
% -------------------------------------------------------------------------

% tic
% % [v,r]=mg(v0,pde,domain,option);
% % option.numit=30;
% % [v,r]=MRR(v0,pde,domain,option);
% [v,r]=bicgstab(v0,pde,domain,option);
% toc
% disp(rms(r(:)))
% tic
% option.numit=5;
% [vv,rr]=MRR(v0,pde,domain,option);
% disp(rms(r(:)))
% toc
% 
% plot(X,v);xlabel('x');title('Numerical solution of Poissons equation')

% -------------------------------------------------------------------------
% Eigenvalue matrices ...
% -------------------------------------------------------------------------

% Eigenvalues of phi_zz+z/c^2*phi=0
% with phi(0)=0;phi(1)=0;

% A=4*ifct(chebdiff(fct(eye(N,N)),2)); % 2x since change in domain to [0,1]
% A(1,:)=0;A(1,1)=0;
% A(end,:)=0;A(end,end)=0;
% 
% X2=(X+1)/2;
% B=-diag(X2);
% B(1,1)=1;B(end,end)=1;
% 
% [V,D,flag]=eigs(A,B,65);
% 
% eigen=diag(D);
% 
% [V2,D2]=eig(A,B);
% 
% eigen2=sort(diag(D2),'descend');

% -------------------------------------------------------------------------

% Eigenvalues of phi_zz+z*lambda*phi=0
% divided both sides by z first ... to get D2/z*phi=lambda*phi
A=4*ifct(chebdiff(fct(eye(N,N)),2)); % 2x since change in domain to [0,1] and 2x for 2nd derivative

X2=(X+1)/2;

A=A./X2;

A(1,:)=0;A(1,1)=1;
A(end,:)=0;A(end,end)=1;

% [v,eigen]=eig(A,'vector');

% Find smallest 3 eigenvectors
[eigenvector,eigenvalue]=eigs(A,3,'smallestabs');
eigenvalue=diag(eigenvalue); % turn into vector

% disregarding the 2 that do not satisfy BCs
eigenvector=eigenvector(:,3);
eigenvalue=eigenvalue(3);

plot(X2,eigenvector)

% Integrals int(phi^2)
int_phi2=sum((eigenvector(1:end-1).^2+eigenvector(2:end).^2)/2.*dx{1}/2); %midpoint rule
int_phi2_spec=clenshaw_curtis(eigenvector.^2)/2; % clenshaw_curtis (divide 2 for domain)

% differentiate
deigenvector=2*ifct(chebdiff(fct(eigenvector),1)); % 2x since change in domain to [0,1]

% phi_z(0)
phiz0=deigenvector(end);

% integrals int(phi_z^2) and int(phi_z^3)
int_phi_z_2=clenshaw_curtis(deigenvector.^2)/2;
int_phi_z_3=clenshaw_curtis(deigenvector.^3)/2;

% % Test stuff ...
% f=exp(X2);
% sum((f(1:end-1).^2+f(2:end).^2)/2.*dx{1}/2) % midpoint rule
% clenshaw_curtis(f.^2)/2 % clenshaw_curtis
% exact=3.194528049465325
