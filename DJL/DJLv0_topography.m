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

function [v0,DJL]=DJLv0_topography(DJL,domain)

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

A=-D2z./N2(z);

% Boundary conditions
A(1,:)=0;A(1,1)=1;
A(end,:)=0;A(end,end)=1;

% Find smallest eigenvectors
[phi,lambda]=eig(A);
lambda=diag(lambda); % turn into vector

% Sort
[lambda,index]=sort(lambda);

% disregarding the last 2 that do not satisfy BCs
phi=phi(:,index(mode+2));
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

gamma=C/2*phiz0/int_phi_z_2;

% -------------------------------------------------------------------------
% fKdVsol for no forcing
% -------------------------------------------------------------------------

% Topography and solution to fKdV
KAI=20;
topography=sech(x{1}/pi*KAI).^2;

if abs(max(topography))==0
    
    % DJL parameters
    u=DJL.u;
    
    % Delta 
    delta=(u-C)/epsilon; % v=c+delta*epsilon
    
    % fKdV solution (after rescaling)
    delta_star=delta/(s*mu^2);

    % KAI such that relative to max 10^-10 at fkdv solution ends
    KAI=sqrt(2/delta_star)*asech(sqrt((2/delta_star)*1e-12*delta_star/2));
    
    % KAI such that absolute min 10^-12 at fkdv solution ends
    % KAI=sqrt(2/delta_star)*asech(sqrt((2/delta_star)*1e-10));
    
    X=x{1}/pi*KAI; % X domain
    DJL.KAI=KAI;

    fKdVsol=delta_star/2*sech(sqrt(delta_star/2)*X).^2; % (X from -KAI to KAI)
    
else

    % pick specific gamma
    gamma_star=-8;

    % Solve for mu
    mu=(gamma*r*alpha/(6*s^2*gamma_star))^(1/4);

    % determine delta_star from hydraulic fall plot
    delta_star=0;
    
    % Solution of fkdv equation
    X=x{1}/pi*KAI;
    fKdVsol=2*sech(X).^2;
    
    % find delta
    delta=delta_star*s*mu^2;

    % find u
    u=delta*epsilon+C;
    
    DJL.u=u;
    DJL.KAI=KAI;
    DJL.mu=mu;

end

v0=epsilon*6*s*mu^2/r*fKdVsol*phi';


end