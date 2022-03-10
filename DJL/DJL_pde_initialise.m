% -------------------------------------------------------------------------
% Set up PDE
% -------------------------------------------------------------------------
% PDE Parameters

Lx=KAI*L/mu/pi;
Ly=1/2;

a=@(X,Y) 1/Lx^2;
b=@(X,Y) 1/Ly^2;
c=@(X,Y) (Y+1)/2/u^2;

% RHS
f=@(X,Y) 0*X;

% -------------------------------------------------------------------------
% Set up BCs
% -------------------------------------------------------------------------
% Boundary conditions for each discretisation (Fourier not used)

% x(1) a1*u+b1*u'=0 x(end) a2*u+b2*u'=0
alpha1{1}=@(X,Y) 1;
beta1{1}=@(X,Y) 0;
alpha2{1}=@(X,Y) 1;
beta2{1}=@(X,Y) 0; 

% x(1) a1*u+b1*u'=0 x(end) a2*u+b2*u'=0
alpha1{2}=@(X,Y) 1;
beta1{2}=@(X,Y) 0;
alpha2{2}=@(X,Y) 1;
beta2{2}=@(X,Y) 0;

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
domain.BC = BC;

a=a(X,Y);
b=b(X,Y);
c=c(X,Y);
f=f(X,Y);

pde.a = a;
pde.b = b;
pde.c = c;
pde.f = f;

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
        pde.f(index{1,i})=0; % Assume 0 for now ...
        pde.f(index{2,i})=0;
    end
    
end
