% Performs multiplication L*u using Fourier collocation methods
% L is the spectral operator for the given PDE in the form
% (-u_xx+au_x+bu+cu^2)_xx+du_yy=f
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
function Lu=fourier_Ku_2d(v,pde,domain)

% kx=domain.k(:,1);
ky=domain.k(:,2);

% a=pde.a;
% b=pde.b;
c=pde.c;
d=pde.d;

% x term
Lu1=fourier_Ku_1d(v,pde,domain);

% y term
Lu2=d.*ifft(-ky.^2.*fft(v'));

Lu=real(Lu1+Lu2');

end