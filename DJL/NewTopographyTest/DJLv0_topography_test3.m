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

% -------------------------------------------------------------------------
% PICK Delta and mu
% -------------------------------------------------------------------------
function DJL=DJLv0_topography_test3(DJL,domain)

N=domain.N;
x=domain.x;

% DJL parameters
delta = DJL.delta;
mu = DJL.mu;
mode = DJL.mode;
KAI = DJL.KAI;
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
% Order 1 fKdVsol for forcing, picking delta and mu
% -------------------------------------------------------------------------

% Domain
X=x{1}/pi*KAI; % -KAI to KAI

% Solve for other variables ...
delta_star=delta/(s*mu^2); % use in fkdv ..

if DJL.soltype==0 % 2sech^2 fKdV solution

    gamma_star=2*delta_star-8;
    B=2*sech(X).^2;

end

if DJL.soltype==1 % fKdV continuation solitary wave

    fprintf('alpha=gamma_star*%f\n',6*s^2*mu^4/(gamma*r))
    gamma_star=input('Input Gamma_star\n');

    % fKdV continuation
    [B_cont,gamma_cont,fKdV,pdefkdv,domainfkdv,optionfkdv]=fkdvsol(DJL,gamma_star,delta_star);

    % Check plot!
    plot(gamma_cont,B_cont(N(1)/2,:))
    title('fKdV continuation at chosen delta_ star')
    xlabel('\gamma')
    ylabel('B(0)')

    % Find index of root
    Bindex=findroot(gamma_cont,gamma_star);
    
    % Plot line on plot
    hold on
    xline(gamma_star)
    hold off

    pause % to check plots ...

    if isempty(Bindex)==1
        disp('No valid solution detected at gamma_star!')
        DJL.v=[];
    elseif length(Bindex)>1
        flag=input('Choose solution (bottom to top)! eg 1,2,3 etc\n');
        Bindex=Bindex(flag);
    end

    % Linear Interpolation
    B=(gamma_star-gamma_cont(Bindex))*(B_cont(:,Bindex+1)-B_cont(:,Bindex)) ...
        /(gamma_cont(Bindex+1)-gamma_cont(Bindex))+B_cont(:,Bindex);

    % Newton iteration!
    % Update variable
    fKdV.gamma=gamma_star;
    pdefkdv.f=-gamma_star*fKdV.topography(X);
    B=NewtonSolve(B,fKdV,pdefkdv,domainfkdv,optionfkdv);

end

% Solve for alpha ...
alpha=gamma_star*6*s^2*mu^4/(gamma*r);
DJL.alpha=alpha;

% back to original compatibility equation
A=6*s*mu^2/r*B;

% find u
u=delta+C;
 
DJL.mu=mu;
DJL.u=u;
DJL.KAI=KAI;
DJL.Lx=2*KAI/mu^2;

v0=A*phi';

% -------------------------------------------------------------------------
% Order epsilon solution in DJL coordinates
% -------------------------------------------------------------------------

% Take only 20 modes ...
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
A_xx=ifft(-(pi/(1/mu*KAI)*domain.k{1}).^2.*fft(A));

b=sech(X).^2; % probably use topography here ...

% A_xx coefficient
a1=-int1./int3.*cN2./(cn2-cN2);

% A^2 coefficient
a2=(cN2./(2*cn2)-2).*cN2./(cn2-cN2).*int2./int3;

% b coefficient
a3=-cn2./(cn2-cN2).*phi_z_0./int3;

% an
an=A_xx*a1+A.^2*a2+alpha*b*a3;

% n=N case
an(:,mode)=0;

% O(e) solution
v1=an*phis';

% Back to zai coordinates
zai=v1+alpha*b*(1-z)';

% solution!
v=v0+zai;
% v=v0;

DJL.v=v;

% -------------------------------------------------------------------------
% CHECKS:
% -------------------------------------------------------------------------
delta_star=delta/(s*mu^2);

gamma_star=gamma*alpha*r/(6*s^2*mu^4);
% CHECKS:
fprintf('delta_star=\n')
disp(delta_star)

fprintf('gamma_star=\n')
disp(gamma_star)

fprintf('delta=\n')
disp(delta)

fprintf('gamma=\n')
disp(gamma)

fprintf('alpha=\n')
disp(alpha)

fprintf('mu=\n')
disp(mu)

fprintf('c=\n')
disp(C)

fprintf('u=\n')
disp(u)

end