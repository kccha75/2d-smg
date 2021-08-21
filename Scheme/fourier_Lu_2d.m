% Performs multiplication L*u using Fourier collocation methods
% L is the spectral operator for the given PDE in the form
% -(au_x)_x-(bu_y)_y+cu=f
%
% Inputs:
% v - best estimate
% pde.a
% pde.b
% pde.c
% domain.k - wave number
%
% Ouputs:
% A - L*u
%
function Lu=fourier_Lu_2d(v,pde,domain)

kx=domain.k{1};
ky=domain.k{2};

a=pde.a;
b=pde.b;
c=pde.c;

% -(au_x)_x term
Lu1=-ifft(1i.*kx.*fft(a.*ifft(1i*kx.*fft(v))));

% -(bu_y)_y term
Lu2=-ifft(1i.*ky.*fft(b'.*ifft(1i*ky.*fft(v'))));

Lu=real(Lu1+Lu2'+c.*v);

end