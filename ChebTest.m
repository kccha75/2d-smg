clear;

grids=5;

% initial grid
N=2^grids+1;
k=(0:N-1)';
x=cos(pi*k/(N-1));

% coarse grid
Nc=2^(grids-1)+1;
kc=(0:Nc-1)';
xc=cos(pi*kc/(Nc-1));

% fine grid
Nf=2^(grids+1)+1;
kf=(0:Nf-1)';
xf=cos(pi*kf/(Nf-1));

% function
f=sech(x).^2;

% Dirichlet boundary
f_dirich=f;
f_dirich(1)=0;f_dirich(end)=0;

% plot(x,f)


% Cheb restrict
fc=cheb_restrict(f_dirich);
% figure;
% plot(xc,fc);

% Cheb prolong
ff=cheb_prolong(f_dirich);
% figure;
% plot(xf,ff)

% Neumann

% f(1)=
a=-1/2*(-1)^(N-1)*x(1)-2*sum((-1).^(N-1+k(2:end-1))./(1+x(2:end-1)).*x(2:end-1))-(2*(N-1)^2+1)/6*x(end);
b=(2*(N-1)^2+1)/6*x(1)+2*sum((-1).^k(2:end-1)./(1-x(2:end-1)))+1/2*(-1)^(N-1)