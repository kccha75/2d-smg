% Performs multiplication L*u using Fourier collocation methods
% L is the spectral operator for the given PDE in the form
% -u_xx+au_x+bu using
%
% Inputs:
% v - best estimate
% pde.a
% pde.b
% domain.k - wave number
% domain.N
%
% Ouputs:
% A - L*u
%

function Lu=fourier_Lu_1d(v,pde,domain)

kx=domain.k{1};

a=pde.a;
b=pde.b;

Lu=real(-ifft(-kx.^2.*fft(v))+a.*ifft(1i.*kx.*fft(v))+b.*v);

end