% Function initialises fKdV parameters for u_xx-delta u+3u^2=-gamma sech^2 
%
% Outputs:
%
% domain - structure of domain (see below)
% option - structure of multigrid options (see below)
% cont_option - structure of continuation options (see below)

function [domain,option,cont_option]=fKdVinitialise()

% Dimension of problem
dim=1;

% Discretisation flag for each dimension
% 1 - Fourier
% 2 - Cheb
discretisation=1;

finestgrid = 9;
coarsestgrid = 3;

% -------------------------------------------------------------------------
% Multigrid Options here
% -------------------------------------------------------------------------

% Number of V-cycles if option is chosen, otherwise number of v-cycles done
% after FMG or num iterations for decent methods
option.numit=5;

% Linear solver / Newton tolerance
option.tol=1e-12;
option.Newtontol=1e-8;
option.Newtonmaxit=20;
option.Newtonsolver='mg';

% Relaxations on the up and down cycle during Multigrid
option.Nd=1;
option.Nu=1;

% Multigrid solver options:'V-cycle' or 'FMG'
option.solver='V-cycle';

% Multigrid scheme: 'Correction' or 'FAS'
option.mgscheme='Correction';

% Operator, coarse grid solver, Relaxation
option.operator=@Lu;
option.coarsegridsolver=@specmatrixsolve;
option.relaxation=@MRR;

% Restriction for pde coefficients
option.restriction=@(vf) fourier_restrict_filtered(vf);

% Restriction for residual and RHS
option.restriction_residual=@(vf) fourier_restrict_filtered(vf);

% Prolongation
option.prolongation=@(vc) fourier_prolong_filtered(vc);

% Preconditioner
option.preconditioner=@FDmatrixsolve;

% Number of preconditioned relaxations
option.prenumit=1;

% Tail tolerance
option.tailtol=1e-6;

% -------------------------------------------------------------------------
% Set up parameters
% -------------------------------------------------------------------------
N=ones(1,max(dim,2)); 
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
            N(i) = 2^(finestgrid)+1;
            k{i} = (0:N(i)-1)';
            x{i} = cos(pi*k{i}/(N(i)-1));
            dx{i} = x{i}(2:end)-x{i}(1:end-1); % is negative but ok :)
            
    end
    
end

X = x{1};

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

option.jacobian=@jacobianfKdV;

% -------------------------------------------------------------------------
% Continuation Options here
% -------------------------------------------------------------------------

% Step size
cont_option.ds=0.01;
cont_option.ds_min=1e-6;
cont_option.ds_max=0.05;

% Iterations
cont_option.N_opt=4;
cont_option.Newtonmaxit=8;
cont_option.steps=500;
