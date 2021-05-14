clear;close all
% Solve Yang NLS using Newton's method
% (-a*u_{xxxx}+c*u_{xx}+d*u)_{xx}-b*u_{yy}
% -------------------------------------------------------------------------
% Domain parameters
% -------------------------------------------------------------------------
Lx = 120*pi;
Ly = 60*pi;
finestgrid = 10;
coarsestgrid = 9;

% -------------------------------------------------------------------------
Nx = 2^10;
Ny = 2^7;

% Spectral Wave numbers
kx = 2*pi/Lx*[0:Nx/2-1 -Nx/2 -Nx/2+1:-1]';
ky = 2*pi/Ly*[0:Ny/2-1 -Ny/2 -Ny/2+1:-1]';
[KX,KY] = ndgrid(kx,ky);

x = Lx*(-Nx/2:Nx/2-1)'/Nx;
y = Ly*(-Ny/2:Ny/2-1)'/Ny;
[X,Y] = ndgrid(x,y);

% -------------------------------------------------------------------------
% PDE parameters
% -------------------------------------------------------------------------

a=ones(Nx,Ny);
b=ones(Nx,Ny);
c=-2*ones(Nx,Ny);
mu=-1.2;
d=mu*ones(Nx,Ny);

% RHS function
f=zeros(Nx,Ny);

% Initial guess
v0=-0.43*sech(0.3*sqrt(X.^2+Y.^2)).*cos(X);

% -------------------------------------------------------------------------
% OPTIONS HERE
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
option.operator=@fourier_Lu_yang_KP;
option.coarsegridsolver=@cg;option.numit=1000;option.coarsegridsolver=@MRR;
option.relaxation=@MRR;
option.restriction=@fourier_restrict_2d_mid;
option.prolongation=@fourier_prolong_2d_mid;

% Preconditioner
option.preconditioner=@fourier_yang_KP_pre;
% Number of precondition relaxations
option.prenumit=1;

% -------------------------------------------------------------------------
% Sort into sctructures
% -------------------------------------------------------------------------
% Assuming constant dx
dx = x(2)-x(1);
dy = y(2)-y(1);

% Sort into structures
domain.Lx = Lx;
domain.Ly = Ly;
domain.Nx = Nx;
domain.Ny = Ny;
domain.x = x;
domain.y = y;

pde.a = a;
pde.b = b;
pde.c = c;
pde.d = d;
pde.f = f;

scheme.kx = kx;
scheme.ky = ky;
scheme.dx = dx;
scheme.dy = dy;

option.finestgrid=finestgrid;
option.coarsestgrid=coarsestgrid;
option.grids=finestgrid-coarsestgrid+1;

% -------------------------------------------------------------------------
% NEWTON HERE
% -------------------------------------------------------------------------

% New b(x) function in Newton
cnew=c-6*v0;
v=v0;

% Error guess (keep at 0)
e0=zeros(Nx,Ny);

tic
for i=1:20
    
    pde.c=c;
    % Initial RHS of linear equation
    pde.f=f-(option.operator(v,pde,scheme)-3*v.^2);
    
    r=rms(rms(pde.f));
    fprintf('Residual Newton = %d\n',r)
    if r<=1e-10
        fprintf('Converged after %d Newton Iterations \n',i-1)
        break
    end
    
    % Solve linear equation
    pde.c=cnew;

%     option.tol=1e-1*r;
%     [e,r]=cg(e0,pde,scheme,option);
    [e,r]=mg2(v,pde,domain,scheme,option);

    % Update correction
    v=v+e;
    
    cnew=c-3*v.^2;
    
end

if i==20
    
    fprintf('Did not converge to required tolerance after %d Newton Iterations\n',i)
    
end 
toc