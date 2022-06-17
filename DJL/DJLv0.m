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

function [v,DJL]=DJLv0(DJL,domain)

N=domain.N;
x=domain.x;

u=DJL.u;
epsilon = DJL.epsilon;
mu=DJL.mu;
mode=DJL.mode;

N2=DJL.N2;

% -------------------------------------------------------------------------
% Eigenvalues of phi_zz+N^2(z)*lambda*phi=0
% -------------------------------------------------------------------------

% divided both sides by z first to get -D2/z*phi=lambda*phi
D2z=4*ifct(chebdiff(fct(eye(N(2),N(2))),2)); % 2x since change in domain to [0,1] and 2x for 2nd derivative

% z domain [0,1]
z=(x{2}+1)/2;

A=-D2z./N2(z);

% Boundary conditions
A(1,:)=0;A(1,1)=1;
A(end,:)=0;A(end,end)=1;

% Find smallest eigenvectors
[phis,lambdas]=eig(A);
lambdas=diag(lambdas); % turn into vector

% Sort
[lambdas,index]=sort(lambdas,'ascend');
phis=phis(:,index);

% disregarding the last 2 that do not satisfy BCs
phi=phis(:,mode+2);
lambda=lambdas(mode+2);

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
delta_star=delta/(s*mu^2);

% KAI such that relative to max 10^-10 at fkdv solution ends
KAI=sqrt(2/delta_star)*asech(sqrt((2/delta_star)*1e-12*delta_star/2));

% KAI such that absolute min 10^-12 at fkdv solution ends
% KAI=sqrt(2/delta_star)*asech(sqrt((2/delta_star)*1e-10));

X=x{1}/pi*KAI; % X domain

% fKdVsol for no forcing
B=delta_star/2*sech(sqrt(delta_star/2)*X).^2; % (X from -KAI to KAI)

% back to original compatibility equation
A=6*s*mu^2/r*B;

% psi_0 solution
v0=epsilon*A*phi';

DJL.KAI=KAI;

% -------------------------------------------------------------------------
% Test for order 1 solution
% -------------------------------------------------------------------------

% remove the bad solutions
phis=phis(:,3:end);
lambdas=lambdas(3:end)';

% Take only 10 modes ...
phis=phis(:,1:10);
lambdas=lambdas(1:10);

% phi_n'
phi_z=2*ifct(chebdiff(fct(phis),1));

% % phi_n'(0)
% phi_z_0=phi_z(end,:);

% int phi_N * phi_n
int1=1/2*clenshaw_curtis(phis.*phis(:,mode).^2);

% int (phi_N')^2 * phi_n'
int2=1/2*clenshaw_curtis(phi_z.*phi_z(:,mode).^2);

% int (phi_n')^2
int3=1/2*clenshaw_curtis(phi_z.^2);

% c_n/c_N
cn=1./sqrt(lambdas);
cN=1/sqrt(lambdas(mode));
cncN2=(cn/cN).^2;

% beta
beta=-ifft(-domain.k{1}.^2.*fft(A))*int1+A.^2*((1/2*1./cncN2.^2-2).*int2);
beta=beta./int3;

% coefficients
an=beta./(lambdas(mode)-lambdas);

% n=N filter
an(:,mode)=0;

% v1
v1=an*phis';

% solution!

v=v0+epsilon*v1;

end
