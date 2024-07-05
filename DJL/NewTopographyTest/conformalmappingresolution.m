% -------------------------------------------------------------------------
% Function performs conformal mapping from (x,z) with topography to (u,v)
% rectangular coordinates
%
% NOTE: uses Fourier x / Chebyshev z discretisation
% Conformal mapping done on x=(-pi,pi) and z=(0,1) domain
%
% Uses Fast spectral Fourier Cheb Poisson solver

% Script compares number of iterations to grid resolution
% -------------------------------------------------------------------------
clear;

N=(1:7)+5;
numit=zeros(length(N),1);
cr1=zeros(length(N),1);
cr2=zeros(length(N),1);

for jj=1:length(N)

Nx=2^N(jj);
x=2*pi*(-Nx/2:Nx/2-1)'/Nx;
kx=[0:Nx/2-1 -Nx/2 -Nx/2+1:-1]';
dx=x(2)-x(1);

Ny=2^(N(jj)-1)+1;
ky=(0:Ny-1)';
y=cos(pi*ky/(Ny-1));

[xgrid,ygrid]=ndgrid(x,y);

% domain struct for solver
domain.N(1)=Nx;
domain.N(2)=Ny;
domain.k{1}=kx;

% topography parameters
alpha=0.1;
Lx=2*pi;

% -------------------------------------------------------------------------
% Conformal Mapping here
% -------------------------------------------------------------------------
H0=2*pi/Lx; % Calculate to keep aspect ratio correct

h = @(x) alpha*H0*sech(x*pi).^2; % Bump function
% h = @(x) alpha*H0*sech(3/2*(x+pi/3)*pi).^2+alpha*H0*sech(3/2*(x-pi/3)*pi).^2; % Bump function

% Maximum iterations
loops=1000;

% Start iterating!
u=x;
x_old=u;
kinv=[0;1./(1i*kx(2:Nx))]; % mean 0, avoid dividing by 0

for i=1:loops
    
    y_bc=h(x_old); % assume initially x=u 
    L=H0-trapI(y_bc,dx)/(2*pi); % New L value (see boundary integral)
    v=L/2*(y+1); % To set the new domain to be [0 L]
    
    % Solve laplace's equation, poisson's after transformation
    % need to change coefficients a,b define new RHS, homogeneous BCs
    h_uu=real(ifft(-kx.^2.*fft(y_bc)));
    [U,V]=ndgrid(u,v); 
    
    pde.a=1;
    pde.b=(2/L)^2;
    pde.c=0;
    pde.f=-h_uu.*(1-V/L);
%     pde.f=-Lu_2d(y_bc+(H0-y_bc).*V/L,pde,domain); % Alternative ...

    % Solve here
    Y=FastChebFourierPoisson([],pde,domain,'');
    
    % Transform back to original coordinates
    yy=y_bc.*(1-V/L)+H0*V/L+Y;

    % find dy/dv
    dy=2/L*real(ifct(chebdiff(real(fct(transpose(yy))),1)));
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
YY=yy;

numit(jj)=i;

% Cauchy-Riemann
% dx/du
dxdu=1+ifft(1i*kx.*fft(ex));

% dx/dv
dxdv=ifct(2/L*chebdiff(real(fct(transpose(ex))),1));
dxdv=dxdv';

% dz/du
dzdu=ifft(1i*kx.*fft(yy));

% dz/dv
dzdv=ifct(2/L*chebdiff(real(fct(transpose(yy))),1));
dzdv=dzdv';

%L2?
cr1(jj)=sqrt(L/2*clenshaw_curtis(transpose(trapI(abs(real(dxdu)-dzdv).^2,dx))));
cr2(jj)=sqrt(L/2*clenshaw_curtis(transpose(trapI(abs(real(dzdu)+dxdv).^2,dx))));

%L1
% cr1(jj)=clenshaw_curtis(transpose(trapI(abs(real(dxdu)-dzdv),dx)));
% cr2(jj)=clenshaw_curtis(transpose(trapI(abs(real(dzdu)+dxdv),dx)));

% something
% cr1(jj)=sqrt(L/2*clenshaw_curtis(transpose(trapI(abs(dxdu-dzdv).^2,dx))));
% cr2(jj)=sqrt(L/2*clenshaw_curtis(transpose(trapI(abs(dzdu+dxdv).^2,dx))));
end

plot(N,numit)

% % -------------------------------------------------------------------------
% % Jacobian Calculation
% % -------------------------------------------------------------------------
% 
% % dx/du
% dxdu=1+ifft(1i*kx.*fft(ex));
% 
% % dx/dv
% dxdv=ifct(2/L*chebdiff(real(fct(transpose(ex))),1));
% dxdv=dxdv';
% 
% % dz/du
% dzdu=ifft(1i*kx.*fft(yy));
% 
% % dz/dv
% dzdv=ifct(2/L*chebdiff(real(fct(transpose(yy))),1));
% dzdv=dzdv';
% 
% % Check Cauchy-Riemann equations
% figure;
% surf(real(dxdu)-dzdv)
% figure;
% surf(real(dzdu)+dxdv)
% 
% % Jacobian calculation
% jac=real(dzdv).^2+real(dzdu).^2;
% 
% % Mapping plots
% figure('Position',[300 300 600 300]); fsz=15;
% 
% subplot(2,1,1)
% contour(XX,YY,xgrid,50,'Color','#0072BD');
% hold on
% contour(XX,YY,ygrid,50,'Color','#0072BD');
% plot(XX(:,1),h(XX(:,1)))
% xlabel('$x$','interpreter','latex','fontsize',fsz);
% ylabel('$y$','interpreter','latex','fontsize',fsz)
% 
% subplot(2,1,2)
% contour(xgrid,ygrid,XX,50,'Color','#0072BD');
% hold on;
% contour(xgrid,ygrid,YY,50,'Color','#0072BD');
% xlabel('$u$','interpreter','latex','fontsize',fsz);
% ylabel('$v$','interpreter','latex','fontsize',fsz)
% 
% 
