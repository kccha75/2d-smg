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
function Lu=fourier_Ku_2d(v,pde,domain)

kx=domain.k{1};
ky=domain.k{2};

a=pde.a;
b=pde.b;
c=pde.c;
d=pde.d;

% inside bracket terms
Lx=a.*ifft(-kx.^2.*fft(v))+b.*v+c.*v.*v;
% x terms
Lu_x=ifft(-kx.^2.*fft(Lx));

% y term
Lu_y=d.*ifft(-ky.^2.*fft(v'));

Lu=real(Lu_x+Lu_y');

end