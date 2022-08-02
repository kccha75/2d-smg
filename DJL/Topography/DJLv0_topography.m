% Function calculates v0 initial guess for the DJL equation
%
% Inputs:
% 
% domain.N
% domain.x
%
% DJL.epsilon - perturbation 
% DJL.L - topography length scale
% DJL.u - given wave speed
% DJL.mode - DJL solution mode
% DJL.N2 - N^2 function
%
% Outputs:
%
% v0 - DJL perturbation solution
% DJL.KAI -DJL x domain size (sufficiently large)

function [v,DJL]=DJLv0_topography(DJL,domain)

N=domain.N;
x=domain.x;

% DJL parameters
epsilon = DJL.epsilon;
alpha = DJL.alpha;
mode = DJL.mode;
N2 = DJL.N2;

% -------------------------------------------------------------------------
% Eigenvalues of phi_zz+N^2(z)*lambda*phi=0
% -------------------------------------------------------------------------

% divided both sides by z first to get -D2/z*phi=lambda*phi
D2z=4*ifct(chebdiff(fct(eye(N(2),N(2))),2)); % 2x since change in domain to [0,1] and 2x for 2nd derivative

% z domain [0,1]
z=(x{2}+1)/2;

% spectral matrix
A=-D2z./N2(z);

% Boundary conditions
A(1,:)=0;A(1,1)=1;
A(end,:)=0;A(end,end)=1;

% Find eigenvalues and eigenvectors
[phis,lambdas]=eig(A);
lambdas=diag(lambdas); % turn into vector

% Sort
[lambdas,index]=sort(lambdas,'ascend');
phis=phis(:,index);

% Disregarding the last 2 that do not satisfy BCs
phis=phis(:,3:end);
lambdas=lambdas(3:end)';

% Specific interested mode
phi=phis(:,mode);
lambda=lambdas(mode);

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

gamma=C/2*phiz0/int_phi_z_2;

% -------------------------------------------------------------------------
% Order 1 fKdVsol for forcing, picking alpha / gamma_star
% -------------------------------------------------------------------------

% Topography length
KAI=20;

% pick specific gamma
gamma_star=-0.01;

% determine delta_star from hydraulic fall plot
delta_star=1/2*(gamma_star+8);

% Solve for mu
mu=(gamma*r/(6*s^2*gamma_star))^(1/4);
% epsilon=mu^2;
    
% Solution of fkdv equation
X=x{1}/pi*KAI; % -KAI to KAI
B=2*sech(X).^2;
    
% back to original compatibility equation
A=6*s*mu^2/r*B;

% find delta
delta=delta_star*s*mu^2;

% find u
u=delta*epsilon+C;
    
DJL.u=u;
DJL.KAI=KAI;
DJL.mu=mu;
DJL.epsilon=epsilon;

v0=epsilon*A*phi';

% -------------------------------------------------------------------------
% Order epsilon solution
% -------------------------------------------------------------------------

% Take only 10 modes ...
phis=phis(:,1:10);
lambdas=lambdas(1:10);

% phi_n'
phi_z=2*ifct(chebdiff(fct(phis),1));

% % phi_n'(0)
phi_z_0=phi_z(end,:);

% int phi_N * phi_n
int1=1/2*clenshaw_curtis(phis.*phi);

% int (phi_N')^2 * phi_n'
int2=1/2*clenshaw_curtis(dphi.^2.*phi_z);

% int (phi_n')^2
int3=1/2*clenshaw_curtis(phi_z.^2);

% c_n/c_N
cn2=1./lambdas;
cN2=1/lambda;

% A_xx in x domain (KAI/mu)
A_xx=ifft(-(pi/(1/mu*KAI)*domain.k{1}).^2.*fft(A));

b=sech(X).^2;

% beta
beta=-b*cn2./cN2.*phi_z_0-A_xx*int1+A.^2*((cN2./(2*cn2)-2).*int2);

beta=beta./(cn2.*int3);

% coefficients
an=beta./(lambdas(mode)-lambdas);

% n=N case
an(:,mode)=0;

% O(e) solution
v1=an*phis';

% Back to zai coordinates
zai=epsilon^2*(v1+b*(1-z)');

% solution!
v=v0+zai;
% v=v0;
end