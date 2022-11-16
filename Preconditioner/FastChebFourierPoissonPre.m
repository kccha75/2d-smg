% Fast Fourier-Chebyshev Tau Poisson solver for -au_xx-bu_yy=f using 
% pseudospectral methods with Fourier in x and Chebyshev in y
%
% Note: a,b,c must be constants only!!!
%
% See Canuto (2006) for reference to algorithm
%
% Inputs:
%
% RHS - specify RHS, if empty use pde.f
% pde.a
% pde.b
% pde.f
% domain.N
% domain.k{1} - Fourier wave numbers 
%
%
% Outputs:
% u - solution

function u=FastChebFourierPoissonPre(RHS,pde,domain,~)

% Set domain parameters
Nx=domain.N(1);
Ny=domain.N(2);
kx=domain.k{1};

% PDE parameters
a=-1;
b=-1;
c=1;
if isempty(RHS)
    f=pde.f;
else
    f=RHS;
end

% lambda value after FFT
lambda=(c-a*kx.^2)/b;

% Set RHS of equation after FFT
RHS=fft(f)/b;

RHS=fct(transpose(RHS));
RHS=transpose(RHS);

%--------------------------------------------------------------------------
% ODD POINTS (even k, since matlab starts at 1)
%--------------------------------------------------------------------------

% lower diag vector
A1=zeros(Nx,(Ny+1)/2);
% main diag vector
A2=zeros(Nx,(Ny+1)/2);
% upper diag vector
A3=zeros(Nx,(Ny+1)/2);
% RHS
g=zeros(Nx,(Ny+1)/2);

% k odd
k_odd=(3:2:Ny-2);  

% Set diagonals
A2(:,2:end-1)=1-lambda./(2*((k_odd-1).^2-1));
A1(:,2:end-1)=cc(k_odd-3).*lambda./(4*(k_odd-1).*(k_odd-2));
A3(:,2:end-1)=lambda./(4*(k_odd-1).*k_odd);

% Set RHS
g(:,2:end-1)=cc(k_odd-3)./(4*(k_odd-1).*(k_odd-2)).*RHS(:,k_odd-2)- ...
    1./(2*((k_odd-1).^2-1)).*RHS(:,k_odd)+1./(4*(k_odd-1).*k_odd).* ...
    RHS(:,k_odd+2);

% k=0 BC
g(:,1)=0;
A2(:,1)=1;
A3(:,1)=1;

% k=N BC
A2(:,(Ny+1)/2)=1;
A1(:,(Ny+1)/2)=lambda./(4*(Ny-1)*(Ny-2));
g(:,(Ny+1)/2)=1/(4*(Ny-1)*(Ny-2))*RHS(:,Ny-2);

%--------------------------------------------------------------------------
% EVEN POINTS (odd k, since matlab starts at 1)
%--------------------------------------------------------------------------

% lower diag
B1=zeros(Nx,(Ny-1)/2);
% main diag
B2=zeros(Nx,(Ny-1)/2);
% upper diag
B3=zeros(Nx,(Ny-1)/2);
% RHS
h=zeros(Nx,(Ny+1)/2-1);

% k even
k_even=(4:2:Ny-3);

% Set diagonals
B2(:,2:end-1)=1-lambda./(2*((k_even-1).^2-1));
B1(:,2:end-1)=lambda./(4*(k_even-1).*(k_even-2));
B3(:,2:end-1)=lambda./(4*(k_even-1).*k_even);

% Set RHS
h(:,2:end-1)=1./(4*(k_even-1).*(k_even-2)).*RHS(:,k_even-2)- ...
    1./(2*((k_even-1).^2-1)).*RHS(:,k_even)+1./(4*(k_even-1).*k_even).* ...
    RHS(:,k_even+2);

% k=1 BC
h(:,1)=0;
B2(:,1)=1;
B3(:,1)=1;

% k=N
B2(:,(Ny-1)/2)=1;
B1(:,(Ny-1)/2)=lambda/(4*(Ny-1)*(Ny-2));
g(:,(Ny-1)/2)=1/(4*(Ny-1)*(Ny-2))*RHS(:,Ny-2);

% k=N-1
B2(:,(Ny-1)/2)=1;
B1(:,(Ny-1)/2)=lambda/(4*(Ny-2)*(Ny-3));
h(:,(Ny-1)/2)=1/(4*(Ny-2)*(Ny-3))*RHS(:,Ny-3);

%--------------------------------------------------------------------------
% Sort and generate sparse matrix
%--------------------------------------------------------------------------

% Index row odd (main,lower,upper)
rowA=[(1:(Ny+1)/2)';(2:(Ny+1)/2)';(1:(Ny+1)/2-1)';ones((Ny+1)/2-2,1)];

% Index columns odd 
columnA=[(1:(Ny+1)/2)';(1:(Ny+1)/2-1)';(2:(Ny+1)/2)';(3:(Ny+1)/2)'];

% Matrix entry odd 
entryA=[A2 A1(:,2:end) A3(:,1:end-1) ones(Nx,(Ny+1)/2-2)];

% Index row even (main,lower,upper)
rowB=[(1:(Ny-1)/2)';(2:(Ny-1)/2)';(1:(Ny-1)/2-1)';ones((Ny-1)/2-2,1)];

% Index columns even
columnB=[(1:(Ny-1)/2)';(1:(Ny-1)/2-1)';(2:(Ny-1)/2)';(3:(Ny-1)/2)'];

% Matrix entry even
entryB=[B2 B1(:,2:end) B3(:,1:end-1) ones(Nx,(Ny-1)/2-2)];

%--------------------------------------------------------------------------
% Solve ...
%--------------------------------------------------------------------------

u=zeros(Nx,Ny);

for j=1:Nx

    % ODD
    u1_tilde_hat=sparse(rowA,columnA,entryA(j,:),(Ny+1)/2,(Ny+1)/2)\transpose(g(j,:));
    
    % EVEN
    u2_tilde_hat=sparse(rowB,columnB,entryB(j,:),(Ny-1)/2,(Ny-1)/2)\transpose(h(j,:));

    % Combine:
    u(j,1:2:end)=u1_tilde_hat;
    u(j,2:2:end)=u2_tilde_hat;

end

% Inverse transforms
u=real(ifft(transpose(ifct(transpose(u)))));

% Calculate c_k for Chebyshev, where c_k=2 if k=0, c_k=1 if k>=1
function c=cc(k)

    c=2*(k==0)+(k>=1);

end

end