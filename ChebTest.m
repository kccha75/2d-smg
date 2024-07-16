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
f=tanh(x);

% Dirichlet boundary
f(1)=0;f(end)=0;

plot(x,f)


% Cheb restrict
fc=cheb_restrict(f);
figure;
plot(xc,fc);

% Cheb prolong
ff=cheb_prolong(f);
figure;
plot(xf,ff)