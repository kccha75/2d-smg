% Function initialises DJL parameters for A_xx+A_zz+N^2(z-A)*A/u^2=0 
%
% Outputs:
%
% domain - structure of domain (see below)
% option - structure of multigrid options (see below)
% cont_option - structure of continuation options (see below)

function [domain,option,cont_option]=DJLinitialise()

% Dimension of problem
dim=2;

% Discretisation flag for each dimension
% 1 - Fourier
% 2 - Cheb
discretisation=[1 2];

finestgrid = 6;
coarsestgrid = 3;

% -------------------------------------------------------------------------
% Multigrid Options here
% -------------------------------------------------------------------------

% Number of V-cycles if option is chosen, otherwise number of v-cycles done
% after FMG
option.num_vcycles=5;

% Linear solver / Newton tolerance
option.tol=1e-12;
option.Newtontol=1e-10;
option.Newtonmaxit=20;

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
X=x;
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
            N(i) = 2^(finestgrid)+1;
            k{i} = (0:N(i)-1)';
            x{i} = cos(pi*k{i}/(N(i)-1));
            dx{i} = x{i}(1:end-1)-x{i}(2:end); % due to x(1)=1, x(end)=-1
            
    end
    
end

[X{1},X{2}] = ndgrid(x{1},x{2});

% -------------------------------------------------------------------------
% Sort into structures
% -------------------------------------------------------------------------
domain.dim = dim;
domain.discretisation = discretisation;

domain.N = N;
domain.k = k;
domain.dx = dx;

option.finestgrid=finestgrid;
option.coarsestgrid=coarsestgrid;
option.grids=finestgrid-coarsestgrid+1;

domain.x = x;
domain.X = X;

% -------------------------------------------------------------------------
% Jacobian option
% -------------------------------------------------------------------------

option.jacobian=@jacobianDJL;

% -------------------------------------------------------------------------
% Continuation Options here
% -------------------------------------------------------------------------

cont_option=option;

% Step size
cont_option.ds=0.001;
cont_option.ds_min=1e-6;
cont_option.ds_max=0.001;

% Iterations
cont_option.N_opt=4;
cont_option.Newtonmaxit=8;
cont_option.steps=200;
