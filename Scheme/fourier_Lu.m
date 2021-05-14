% Performs multiplication L*u using Fourier collocation methods
% L is the spectral operator for the given PDE in the form
% -(au_x)_x-(bu_y)_y+cu=f
%
% Inputs:
% v - best estimate
% pde.a
% pde.b
% pde.c
% pde.f
% scheme.kx - wave number in x
% scheme.ky - wave number in y
%
% Ouputs:
% A - L*u
%
function Lu=fourier_Lu(v,pde,scheme)

a=pde.a;
b=pde.b;
c=pde.c;

% -(au_x)_x term
Lu1=-ifft(1i.*scheme.kx.*fft(a.*ifft(1i*scheme.kx.*fft(v))));

% -(bu_y)_y term
Lu2=-ifft(1i.*scheme.ky.*fft(b'.*ifft(1i*scheme.ky.*fft(v'))));

Lu=real(Lu1+Lu2'+c.*v);

end