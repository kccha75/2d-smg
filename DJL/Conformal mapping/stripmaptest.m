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
coarsestgrid = 6;

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

%--------------------------------------------------------------------------

%
% add:
% documentation
% calculation of Jacobian = x_u^2+y_v^2
%
NI = 10;
%
N = 512;
M = 128;
H0 = 1;
h = @(x) .1*(1-tanh((x-.5)/.1).^2)+.1*(1-tanh((x+.5)/.1).^2);
% h = @(x) .5*(1-tanh(x/.25).^2);
% h = @(x) .25*exp(-(x/.25).^2);

u = -pi+2*pi*(0:N-1)/N;
K = 1i*[0:N/2-1 -N/2:-1];
Kinv = [.0 1./K(2:N)];
xi = .5*(1-cos(pi*(0:M)/M));
L1 = [0:M];
%
xb = u;
for nit = 1:NI
	yb = h(xb);
	L = H0-1/N*sum(yb);
	v = L*xi;
	LL = -2/L*L1;
	[U,V] = ndgrid(u,v);
	F = zeros(N,M+1);
	for ix = 1:N
		F(ix,:) = yb(ix)+(H0-yb(ix)).*V(ix,:)/L;
	end
	G = zeros(N,M+1);
	%
	% initialize the RHS
	%
	F = reshape(F,[N*(M+1),1]);
	nabF = -HelmholtzOpFC(F,G,u,v,K',LL,N,M+1);
	nabF = reshape(nabF,[N,M+1]);
	nabF(:,1:M:M+1) = .0;
	%


	[Y,R] = SpecMG(1,u,v,G,nabF,@HelmholtzOpFD,@HelmholtzOpFC);
    pde.a=4;
    pde.b=1/(L/2)^2;
    pde.c=0;
    pde.f=-nabF;
%     [Y,R]=mg(v0,pde,domain,option);

%     disp(rms(rms(Y-Y1)))


	Yn = Y+reshape(F,[N,M+1]);
	Yt = Yn';
	Dy = ChebBackward(ChebDeriv(ChebForward(Yt),LL,1));
	Dx = Dy(1,:);
	Ix = mean(Dx)*(u+pi)-pi+real(ifft(Kinv.*fft(Dx)));
%    Ix = 2*pi/N*(cumsum(Dx)-Dx(1))-pi;
    fprintf('Change in x = %g \n',norm(Ix-xb));
	if norm(Ix-xb) < 1.e-8
		break;
	end
	xb = Ix;
end
%
Dx = Dy(M+1,:);
% fo = mean(Dx);
xt = mean(Dx)*(u+pi)-pi+real(ifft(Kinv.*fft(Dx)));
% xt = 2*pi/N*(cumsum(Dx)-Dx(1))-pi;
xbp = xb-u;
xtp = xt-u;
F = zeros(N,M+1);
for ix = 1:N
	F(ix,:) = xbp(ix)+(xtp(ix)-xbp(ix)).*V(ix,:)/L;
end
F = reshape(F,[N*(M+1),1]);
nabF = -HelmholtzOpFC(F,G,u,v,K',LL,N,M+1);
nabF = reshape(nabF,[N,M+1]);
nabF(:,1:M:M+1) = .0;
% X = HelmholtzOpFD(G,reshape(nabF,[N*(M+1),1]),u,v,K',LL,N,M+1);
[X,R] = SpecMG(1,u,v,G,nabF,@HelmholtzOpFD,@HelmholtzOpFC);
Xn = X+reshape(F,[N,(M+1)]);
XF = Xn+U;
YF = Yn;
%
clf
subplot(211)
contour(XF,YF,U,40); 
hold on 
contour(XF,YF,V,40)
xlabel('x')
ylabel('y')
plot(xb,h(xb))
subplot(212)
contour(U,V,XF,40); 
hold on 
contour(U,V,YF,40)
xlabel('u')
ylabel('v')
