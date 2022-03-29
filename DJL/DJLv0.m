% -------------------------------------------------------------------------
% Eigenvalues of phi_zz+z*lambda*phi=0
% -------------------------------------------------------------------------
% divided both sides by z first ... to get D2/z*phi=lambda*phi
A=-4*ifct(chebdiff(fct(eye(N(2),N(2))),2)); % 2x since change in domain to [0,1] and 2x for 2nd derivative

x2=(x{2}+1)/2;

A=A./x2;

A(1,:)=0;A(1,1)=1;
A(end,:)=0;A(end,end)=1;

% Find smallest eigenvectors
[eigenvector,eigenvalue]=eigs(A,3,'smallestabs');
eigenvalue=diag(eigenvalue); % turn into vector

% disregarding the 2 that do not satisfy BCs
eigenvector=eigenvector(:,3);
eigenvalue=eigenvalue(3);

C=1/sqrt(eigenvalue);

% -------------------------------------------------------------------------
% Calculate integrals (coefficients for fKdV equation)
% -------------------------------------------------------------------------
% Integrals int(phi^2)
int_phi_2=clenshaw_curtis(eigenvector.^2)/2; % clenshaw_curtis (divide 2 for domain)

% differentiate
deigenvector=2*ifct(chebdiff(fct(eigenvector),1)); % 2x since change in domain to [0,1]

% phi_z(0)
phiz0=deigenvector(end);

% integrals int(phi_z^2) and int(phi_z^3)
int_phi_z_2=clenshaw_curtis(deigenvector.^2)/2;
int_phi_z_3=clenshaw_curtis(deigenvector.^3)/2;

r=3*C/2*int_phi_z_3/int_phi_z_2;
s=C/2*int_phi_2/int_phi_z_2;
delta=(u-C)/epsilon; %v=c+delta*epsilon
gamma=C/2*phiz0/int_phi_z_2;

% fKdV solution (after rescaling)
delta_star=delta*L^2/s;

% KAI such that relative to max 10^-10 at fkdv solution ends
% KAI=sqrt(2/delta_star)*asech(sqrt((2/delta_star)*1e-12*delta_star/2));

% KAI such that absolute min 10^-12 at fkdv solution ends
KAI=sqrt(2/delta_star)*asech(sqrt((2/delta_star)*1e-10));
% KAI=100;
XX=x{1}/pi*KAI; % X domain

fKdVsol=delta_star/2*sech(sqrt(delta_star/2)*XX).^2; % (X from -KAI to KAI)
v0=epsilon*6*s/(r*L^2)*fKdVsol*eigenvector';

