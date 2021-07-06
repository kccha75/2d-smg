% Performs multiplication L*u using Fourier collocation methods
% L is the spectral operator for the given PDE in the form
% (au_xx+bu+cu^2)_xx+du_yy
%
% Inputs:
% v - best estimate
% pde.a
% pde.b
% pde.c
% pde.d
% domain.k - wave number
%
% Ouputs:
% A - L*u
%
function Lu=fourier_KPu_2d_dealias2(v,pde,domain)

kx=domain.k{1};
ky=domain.k{2};

a=pde.a;
b=pde.b;
c=pde.c;
d=pde.d;

[KX,KY]=ndgrid(kx,ky);

a_hat=fft2(a);
b_hat=fft2(b);
d_hat=fft2(d);
v_hat=fft2(v);
% assumed c is constant ...

% inside bracket terms in fourier
Lx=fourier_dealias_2d(a_hat,-KX.^2.*v_hat)+fourier_dealias_2d(b_hat,v_hat)+c.*fourier_dealias_2d(v_hat,v_hat);

% x terms in fourier
Lu_x=-KX.^2.*Lx;

% y term in fourier
Lu_y=-KY.^2.*fourier_dealias_2d(d_hat,v_hat);

Lu=real(ifft2(Lu_x+Lu_y));

end