% Performs multiplication L*u using Fourier collocation methods
% L is the spectral operator for the given PDE in the form
% -(au_x)_x-(bu_y)_y+cu=f
%
% Uses midpoints method for au_x and bu_y
%
% Inputs:
% v - best estimate
% pde.a
% pde.b
% pde.c
% pde.f
% scheme.KX - wave number in x
% scheme.KY - wave number in y
%
% Ouputs:
% A - L*u
%
function Lu=fourier_Lu_mid(u,pde,scheme)

a=pde.a;
b=pde.b;
c=pde.c;

[Nx,Ny]=size(u);

% Shifted midpoint at i+1/2
a_mid=real(ifft(fft(a).*exp(pi*1i*scheme.kx/Nx)));
b_mid=real(ifft(fft(b').*exp(pi*1i*scheme.ky/Ny)));

% -(au_x)_x term
Lu1=-ifft(1i*scheme.kx.*fft(a_mid.*ifft(1i*scheme.kx.*fft(u).*exp(pi*1i*scheme.kx/Nx))).*exp(-pi*1i*scheme.kx/Nx));

% -(bu_y)_y term
Lu2=-ifft(1i*scheme.ky.*fft(b_mid.*ifft(1i*scheme.ky.*fft(u').*exp(pi*1i*scheme.ky/Ny))).*exp(-pi*1i*scheme.ky/Ny));

Lu=real(Lu1+Lu2'+c.*u);

end