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

function [v,DJL]=DJLv0cordsforce(DJL,domain)

N=domain.N;
x=domain.x;

% DJL parameters
% u=DJL.u;
epsilon = DJL.epsilon;
alpha=DJL.alpha;
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
% delta=(u-C)/epsilon; % v=c+delta*epsilon delta pick here ?
gamma=C/2*phiz0/int_phi_z_2;

% in DJL coordinates
r_hat=r;
s_hat=s*epsilon/mu^2;
gamma_hat=gamma*epsilon^2/alpha;
% delta_hat=epsilon*delta;

% fKdV solution (after rescaling)
% delta_hat_star=delta_hat/(s*mu^2);

% KAI set
KAI_hat=20;

gamma_star=gamma*r/(6*s^2);
disp('gamma_star')
disp(gamma_star)

delta_star=1/2*(gamma_star+8); % from relation!
disp('delta_star')
disp(delta_star)

delta=delta_star*s;
disp('delta')
disp(delta)

u=delta*epsilon+C;
disp('u')
disp(u)

X=x{1}/pi*KAI_hat; % X domain

% fKdVsol for no forcing
B_hat=2*sech(X).^2;

% back to original compatibility equation
A_hat=6*s*epsilon/r*B_hat;

% O(1) solution
v0=A_hat*phi';

DJL.u=u;
DJL.KAI=KAI_hat;
DJL.Lx=2*KAI_hat/mu;

% -------------------------------------------------------------------------
% Order epsilon solution
% -------------------------------------------------------------------------

%
A=A_hat/epsilon;

% Take only 10 modes ...
phis=phis(:,1:20);
lambdas=lambdas(1:20);

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
A_xx=ifft(-(pi/KAI_hat*domain.k{1}).^2.*fft(A));

b=sech(X).^2;

% A_xx coefficient
a1=-int1./int3.*cN2./(cn2-cN2);

% A^2 coefficient
a2=(cN2./(2*cn2)-2).*cN2./(cn2-cN2).*int2./int3;

% b coefficient
a3=-cn2./(cn2-cN2).*phi_z_0./int3;

% an
an=epsilon^2*(A_xx*a1+A.^2*a2)+alpha*b*a3;

% n=N case
an(:,mode)=0;

% O(e) solution
v1=an*phis';

% Back to zai coordinates
zai=v1+alpha*b*(1-z)';

% solution!
v=v0+zai;
% v=v0;
end