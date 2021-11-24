clear;close all
% -------------------------------------------------------------------------
% Solve KP using Newton's method in the form (-u_xx+au_x+bu+cu^2)_xx+du_yy=f_xx
% Note: there is a problem if coarsestgrid is 6 or lower 
% Also doesn't work!!! :)
% -------------------------------------------------------------------------
% INPUT PARAMETERS
% -------------------------------------------------------------------------
L(1) = 120;
L(2) = 120;
finestgrid = 10;
coarsestgrid = 7;

% PDE parameters

a=@(X,Y) 0;
b=@(X,Y) 0;
c=@(X,Y) -3;
d=@(X,Y) 1;

% RHS function
f=@(X,Y) -8*(-2*sech(X).^4+4*sech(X).^2.*tanh(X).^2);

% Initial guess
v0=@(X,Y) 2.01*sech(X).^2;
ue=@(X,Y) 2*sech(X).^2;

% -------------------------------------------------------------------------
% Set up parameters
% -------------------------------------------------------------------------
N(1) = 2^finestgrid;
N(2) = 2^finestgrid;

% Spectral Wave numbers
k(:,1) = 2*pi/L(1)*[0:N(1)/2-1 -N(1)/2 -N(1)/2+1:-1]';
k(:,2) = 2*pi/L(2)*[0:N(2)/2-1 -N(2)/2 -N(2)/2+1:-1]';
[KX,KY] = ndgrid(k(:,1),k(:,2));

x(:,1) = L(1)*(-N(1)/2:N(1)/2-1)'/N(1);
x(:,2) = L(2)*(-N(2)/2:N(2)/2-1)'/N(2);
[X,Y] = ndgrid(x(:,1),x(:,2));

a=a(X,Y);
b=b(X,Y);
c=c(X,Y);
d=d(X,Y);
f=f(X,Y);

v0=v0(X,Y);
ue=ue(X,Y);
% -------------------------------------------------------------------------
% Multigrid Options here
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
option.operator=@fourier_Ku_2d;
option.coarsegridsolver=@cg;
option.relaxation=@MRR;
option.restriction=@fourier_restrict_2d_filtered;
option.prolongation=@fourier_prolong_2d_filtered;

% Preconditioner
option.preconditioner=@fourier_yang_pre_2d;
% Number of precondition relaxations
option.prenumit=0;

% -------------------------------------------------------------------------
% Sort into sctructures
% -------------------------------------------------------------------------
% Assuming constant dx
dx(1) = x(2,1)-x(1,1);
dx(2) = x(2,2)-x(1,2);

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


% % -------------------------------------------------------------------------
% % NEWTON HERE
% % -------------------------------------------------------------------------


% 1d options
option1d=option;
option1d.operator=@fourier_Ku_1d;
option1d.coarsegridsolver=@cg;
option1d.relaxation=@MRR;
option1d.restriction=@fourier_restrict_filtered;
option1d.prolongation=@fourier_prolong_filtered;

% Initial guess
v=v0;

% Jacobian here (linearized equation)
bnew=b+2*c.*v0;

% Initial error guess
e0=zeros(N);e=e0;lin=e0;

tic
for i=1:20
    domain.k(1,1)=0;domain.k(1,2)=0;
    e=e0;
    pde.b=b;
    
    % Initial RHS of linear equation (nonlinear residual)
    pde.f=f-(option.operator(v,pde,domain));
    
    r=rms(rms(pde.f));
    fprintf('Residual Newton = %d\n',r)
    if r<=1e-10
        fprintf('Converged after %d Newton Iterations \n',i-1)
        break
    end
    
    % Solve Newton iteration (fixed point)
    % 1D slices
%    for jj=1:5
    pde.f=f-(option.operator(v,pde,domain))-e;
    for ii=1:N(1)
        
        pde1.a=a;
        pde1.b=bnew(:,ii);
        pde1.c=c;
        pde1.d=d;
        pde1.f=pde.f(:,ii);
        
        % Solve linear equation
        pde1.c=0;
        e(:,ii)=fourier_matrixsolve_1d_KP(e0(:,ii),pde1,domain,option1d);

        % Linear part
        lin(:,ii)=-option1d.operator(e(:,ii),pde1,domain);
        
    end

	% Nonlinear part
	pde.c=c;
	nonlin=f-option.operator(v,pde,domain);
        
	% Final e
    k(1,2)=100;
	e=1/pde.d*ifft(-1./k(:,2).^2.*fft((nonlin+lin)'));
%    end
    v=v+e';
    % Update Jacobian
    bnew=b+2.*c.*v;

end

if i==20
    
    fprintf('Did not converge to required tolerance after %d Newton Iterations\n',i)
    
end 

