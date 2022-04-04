% Function calculates v0 initial guess for the DJL equation for N^2(e)=e
%
% Inputs:
% 
% domain.N
% domain.x
%
% DJL.epsilon - perturbation 
% DJL.L - topography length scale
% DJL.u - given wave speed
% DJL.KAI - DJL x domain size (sufficiently large)
% DJL.mode - DJL solution mode
%
% Outputs:
%
% v0 - DJL perturbation solution

function v0=DJLv0(DJL,domain)

N=domain.N;
x=domain.x;

u=DJL.u;
epsilon = DJL.epsilon;
L = DJL.L;
KAI = DJL.KAI;
mode=DJL.mode;

% -------------------------------------------------------------------------
% Eigenvalues of phi_zz+z*lambda*phi=0
% -------------------------------------------------------------------------

% divided both sides by z first to get -D2/z*phi=lambda*phi
D2z=4*ifct(chebdiff(fct(eye(N(2),N(2))),2)); % 2x since change in domain to [0,1] and 2x for 2nd derivative

% z domain [0,1]
z=(x{2}+1)/2;

A=-D2z./z;

% Boundary conditions
A(1,:)=0;A(1,1)=1;
A(end,:)=0;A(end,end)=1;

% Find smallest eigenvectors
[phi,lambda]=eigs(A,mode+2,'smallestabs');
lambda=diag(lambda); % turn into vector

% disregarding the last 2 that do not satisfy BCs
phi=phi(:,mode+2);
lambda=lambda(mode+2);

% lambda=1/C^2
C=1/sqrt(lambda);

% -------------------------------------------------------------------------
% Calculate integrals (coefficients for fKdV equation)
% -------------------------------------------------------------------------

% Integrals int(phi^2)
int_phi_2=clenshaw_curtis(phi.^2)/2; % clenshaw_curtis (divide 2 for domain)

% differentiate
dphi=2*ifct(chebdiff(fct(phi),1)); % 2x since change in domain to [0,1]

% phi_z(0)
phiz0=dphi(end);

% integrals int(phi_z^2) and int(phi_z^3)
int_phi_z_2=clenshaw_curtis(dphi.^2)/2;
int_phi_z_3=clenshaw_curtis(dphi.^3)/2;

% fKdV coefficients
r=3*C/2*int_phi_z_3/int_phi_z_2;
s=C/2*int_phi_2/int_phi_z_2;
delta=(u-C)/epsilon; % v=c+delta*epsilon
gamma=C/2*phiz0/int_phi_z_2;

% fKdV solution (after rescaling)
delta_star=delta*L^2/s;

% KAI such that relative to max 10^-10 at fkdv solution ends
% KAI=sqrt(2/delta_star)*asech(sqrt((2/delta_star)*1e-12*delta_star/2));

% KAI such that absolute min 10^-12 at fkdv solution ends
% KAI=sqrt(2/delta_star)*asech(sqrt((2/delta_star)*1e-10));

X=x{1}/pi*KAI; % X domain

% fKdVsol for no forcing
fKdVsol=delta_star/2*sech(sqrt(delta_star/2)*X).^2; % (X from -KAI to KAI)
v0=epsilon*6*s/(r*L^2)*fKdVsol*phi';

end
