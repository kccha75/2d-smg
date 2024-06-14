% -------------------------------------------------------------------------
% Function performs conformal mapping from (x,z) with topography to (u,v)
% rectangular coordinates
%
% NOTE: uses Fourier x / Chebyshev z discretisation
% Conformal mapping done on x=(-pi,pi) and z=(0,1) domain
%
% Uses Fast spectral Fourier Cheb Poisson solver
% -------------------------------------------------------------------------

[domain,option,cont_option]=DJLinitialise_topography();

% topography parameters
alpha=0.01;
topography=@(X) sech(X).^2;
Lx=2*pi;

N=domain.N;
x=domain.x;
k=domain.k;
dx=domain.dx;

% -------------------------------------------------------------------------
% Conformal Mapping here
% -------------------------------------------------------------------------
H0=2*pi/Lx; % Calculate to keep aspect ratio correct

h = @(x) alpha*H0*topography(x*2*pi); % Bump function

% Maximum iterations
loops=100;

% Start iterating!
u=x{1};
x_old=u;
kinv=[0;1./(1i*k{1}(2:N(1)))]; % mean 0, avoid dividing by 0

for i=1:loops
    
    y_bc=h(x_old); % assume initially x=u 
    L=H0-trapI(y_bc,dx{1})/(2*pi); % New L value (see boundary integral)
    v=L/2*(x{2}+1); % To set the new domain to be [0 L]
    
    % Solve laplace's equation, poisson's after transformation
    % need to change coefficients a,b define new RHS, homogeneous BCs
    h_uu=real(ifft(-k{1}.^2.*fft(y_bc)));
    [U,V]=ndgrid(u,v); 
    
    pde.a=1;
    pde.b=(2/L)^2;
    pde.c=0;
    pde.f=-h_uu.*(1-V/L);
%     pde.f=-Lu_2d(y_bc+(H0-y_bc).*V/L,pde,domain); % Alternative ...

    % Solve here
    Y=FastChebFourierPoisson([],pde,domain,option);
    
    % Transform back to original coordinates
    y=y_bc.*(1-V/L)+H0*V/L+Y;

    % find dy/dv
    dy=2/L*real(ifct(chebdiff(real(fct(transpose(y))),1)));
    dy=dy';

    % at Bottom boundary
    dybc=dy(:,end);

    % Solve for correction e=int(dy/dv-1)du
    epsilon=real(ifft(kinv.*fft(dybc-1))); % mean 0 solution
    epsilon(1)=0; % Boundary condition

    % Find new x, x=u+e
    x_new=u+epsilon;

    % Compare old x with new x, break if tol met, loop otherwise
    if rms(x_new - x_old) < 1e-10
        fprintf('x diff is %d\n',rms(x_new - x_old))
        fprintf('Reached tolerance after %d iterations!\n',i)
        break;
    end
    
    fprintf('x diff is %d\n',rms(x_new - x_old))
    
    % Check final loop
    if i~=loops
        x_old=x_new;
    end
    
end

% x=u+e
ex=real(ifft(kinv.*fft(dy-1)));

% Output coordinates
XX=U+ex;
YY=y;

domain.YY=YY;
domain.XX=XX;
% DJL.Ly=L;

% -------------------------------------------------------------------------
% Jacobian Calculation
% -------------------------------------------------------------------------

% dx/du
dxdu=1+ifft(1i*domain.k{1}.*fft(ex));

% dx/dv
dxdv=ifct(2/L*chebdiff(real(fct(transpose(ex))),1));
dxdv=dxdv';

% dz/du
dzdu=ifft(1i*domain.k{1}.*fft(y));

% dz/dv
dzdv=ifct(2/L*chebdiff(real(fct(transpose(y))),1));
dzdv=dzdv';

% Check Cauchy-Riemann equations
surf(real(dxdu)-dzdv)
figure;
surf(real(dzdu)+dxdv)

% Jacobian calculation
jac=real(dzdv).^2+real(dzdu).^2;
