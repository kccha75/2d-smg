clear;close all
% Solve Yang NLS using Newton's method
% Note: there is a problem if coarsestgrid is 6 or lower 
% 
% -------------------------------------------------------------------------
% Domain parameters
% -------------------------------------------------------------------------
Lx = 20 * pi;
Ly = 20 * pi;
finestgrid = 8;
coarsestgrid = 7;

% -------------------------------------------------------------------------
Nx = 2^finestgrid;
Ny = 2^finestgrid;

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

V0=6;mu=-4.56;
% c(x) function
c=V0*(sin(X).^2+sin(Y).^2)+mu;
c=V0*(alias_2d_test(sin(X),sin(X))+alias_2d_test(sin(Y),sin(Y)))+mu;
% RHS function
f=zeros(Nx,Ny);

% Initial guess
v0=0.56*sech(0.5*sqrt(X.^2+Y.^2)).*cos(X).*cos(Y);

% -------------------------------------------------------------------------
% OPTIONS HERE
% -------------------------------------------------------------------------

% Number of V-cycles if option is chosen, otherwise number of v-cycles done
% after FMG
option.num_vcycles=100;

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
option.operator=@fourier_Lu_mid;
option.coarsegridsolver=@cg;
option.relaxation=@MRR;
option.restriction=@fourier_restrict_2d_mid;
option.prolongation=@fourier_prolong_2d_mid;

% Preconditioner
option.preconditioner=@fourier_yang_pre;
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
cnew=c+3*v0.^2;
% cnew=c+3*alias_2d_test(v0,v0);
v=v0;

% Error guess (keep at 0)
e0=zeros(Nx,Ny);

tic
for i=1:20
    
    pde.c=c;
    % Initial RHS of linear equation
    pde.f=f-(option.operator(v,pde,scheme)+v.^3);
   
    r=rms(rms(pde.f));
    fprintf('Residual Newton = %d\n',r)
    if r<=1e-10
        fprintf('Converged after %d Newton Iterations \n',i-1)
        break
    end
    
    % Solve linear equation
    pde.c=cnew;

    option.tol=1e-1*r;
%     [e,r]=cg(e0,pde,scheme,option);
%     e=fourier_matrixsolve_mid(e0,pde,scheme,option);
    [e,r]=mg(e0,pde,domain,scheme,option);

    % Update correction
    v=v+e;
    
    cnew=c+3*v.^2;
%     cnew=c+3*alias_2d_test(v,v);
end

if i==20
    
    fprintf('Did not converge to required tolerance after %d Newton Iterations\n',i)
    
end 
toc