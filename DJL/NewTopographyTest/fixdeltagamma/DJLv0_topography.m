% Function calculates v0 initial guess for the DJL equation using
% perturbation methods to second order
%
% Assumes topography is alpha*sech(x)^2 form
%
% -------------------------------------------------------------------------
% NOTE: PICK alpha, mu, delta_star SOLVE FOR delta, u
% -------------------------------------------------------------------------
%
% Inputs:
% 
% domain.N - grid points
% domain.x - vector structure x{1} ,x{2}
%
% DJL.delta_star - fkdv parameter
% DJL.gamma_star - fkdv parameter
% DJL.mode - DJL solution mode
% DJL.KAI - DJL domain in nondimensionalised
% DJL.N2 - N^2 function
% option.tailtol - tail end tolerance (make sure it is asymptotic)
%
% Outputs:
% DJL.v - perturbation solution
% DJL.alpha - topography height
% DJL.mu - topography length scale
% DJL.u - wave speed
% DJL.Lx - x domain in DJL coordinates (note this is 2*KAI/mu^2)
% fKdV - structures for fkdv
% pdefkdv
% domainfkdv
% optionfkdv

function [DJL,fKdV,pdefkdv,domainfkdv,optionfkdv]=DJLv0_topography(DJL,domain,option)

N=domain.N;
x=domain.x;

% DJL parameters
delta_star=DJL.delta_star;
gamma_star=DJL.gamma_star;
mu=DJL.mu;
mode = DJL.mode;
KAI = DJL.KAI;
N2 = DJL.N2;

tailtol=option.tailtol;

% -------------------------------------------------------------------------
% Eigenvalues of phi_zz+N^2(z)*lambda*phi=0
% -------------------------------------------------------------------------

% divided both sides by z first to get -D2/z*phi=lambda*phi
D2z=4*real(ifct(chebdiff(real(fct(eye(N(2),N(2)))),2))); % 2x since change in domain to [0,1] and 2x for 2nd derivative

% z domain [0,1]
z=(x{2}+1)/2;

% Spectral matrix
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
dphi=2*real(ifct(chebdiff(real(fct(phi)),1))); % 2x since change in domain to [0,1]

% phi_z(0)
phiz0=dphi(end);

% integrals int(phi_z^2) and int(phi_z^3)
int_phi_z_2=clenshaw_curtis(dphi.^2)/2;
int_phi_z_3=clenshaw_curtis(dphi.^3)/2;

% fKdV coefficients
r=3*C/2*int_phi_z_3/int_phi_z_2;
s=C/2*int_phi_2/int_phi_z_2;

gamma=C/2*phiz0/int_phi_z_2;

% Save variables
DJL.C=C;
DJL.phi=phi;
DJL.r=r;
DJL.s=s;
DJL.gamma=gamma;

% -------------------------------------------------------------------------
% Order 1 fKdVsol for forcing, picking delta and mu
% -------------------------------------------------------------------------

% Domain
X=x{1}/pi*KAI; % -KAI to KAI

% Initialise fkdv structures
fKdV.L = 2*DJL.KAI;
fKdV.d = 3;
fKdV.topography = DJL.topography;


if DJL.soltype==0 % 2sech^2 fKdV solution

    delta_star=1/2*(gamma_star+8);
    delta=delta_star*s*mu^2;
    alpha=delta^2*6*gamma_star/(gamma*r*delta_star^2);
    B=2*sech(X).^2;

    % Check boundary!
    if abs(B(1))>tailtol || abs(B(end))>tailtol
        fprintf('B boundary larger than %d!\n',tailtol)
        return
    end

    % Save variables
    fKdV.B=B;
    fKdV.gamma=gamma_star;
    fKdV.delta=delta_star;

    % Initialise
    [domainfkdv,optionfkdv,~]=fKdVinitialise();

    % Initialise PDE
    [fKdV,pdefkdv,domainfkdv]=fKdVpdeinitialise(fKdV,domainfkdv);

end

if DJL.soltype==1 % fKdV continuation solitary wave

    delta=delta_star*s*mu^2;
    alpha=delta^2*6*gamma_star/(gamma*r*delta_star^2);
    
    % fKdV continuation
    [B_cont,gamma_cont,fKdV,pdefkdv,domainfkdv,optionfkdv]=fkdvsol(DJL,gamma_star,delta_star);

    % Check plot!
    plot(gamma_cont,B_cont(N(1)/2+1,:))
    title('fKdV continuation at chosen delta_ star')
    xlabel('\gamma')
    ylabel('B(0)')

    % Unfolded ...
    figure;
    plot(gamma_cont,trapI(B_cont,domainfkdv.dx{1}))
    title('fKdV continuation at chosen delta_ star')
    xlabel('\gamma')
    ylabel('P')

    % Find index of root
    Bindex=findroot(gamma_cont,gamma_star);
    
    % Plot line on plot
    hold on
    xline(gamma_star)
    hold off

    pause % to check plots ...

    % Pick solution type
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

    % Check boundary!
    if abs(B(1))>tailtol || abs(B(end))>tailtol
        fprintf('B boundary larger than %d!\n',tailtol)
        return
    end
    % Newton iteration!
    % Update variable
    fKdV.gamma=gamma_star;
    pdefkdv.f=-gamma_star*fKdV.topography(X);
    B=NewtonSolve(B,fKdV,pdefkdv,domainfkdv,optionfkdv);

    % Save variables
    fKdV.B=B;

end

% Solve for delta
DJL.delta=delta;
DJL.alpha=alpha;

% back to original compatibility equation
A=6*s*mu^2/r*B;

% find u
u=delta+C;

% Set DJL parameters
DJL.u=u;
DJL.Lx=2*KAI/mu;

% 1st order perturbation
v0=A*phi';

% -------------------------------------------------------------------------
% Order epsilon solution in DJL coordinates
% -------------------------------------------------------------------------

% Take modes to nyquist frequency
phis=phis(:,1:(length(z)+1)/2);
lambdas=lambdas(1:(length(z)+1)/2);

% phi_n'
phi_z=2*real(ifct(chebdiff(real(fct(phis)),1)));

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

% Order 2 solution!
v=v0+zai;

DJL.phis=phis;
DJL.a1=a1;
DJL.a2=a2;
DJL.a3=a3;
DJL.b=b;
DJL.v=v;

% -------------------------------------------------------------------------
% CHECKS:
% -------------------------------------------------------------------------

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