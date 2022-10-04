% Performs multiplication L*u using Fourier collocation methods
% L is the spectral operator for the given PDE in the form
% (-u_xx+au_x+bu+cu^2)_xx using
%
% Inputs:
% v - best estimate
% pde.a
% pde.b
% pde.c
% domain.k - wave number
% domain.N
%
% Ouputs:
% A - L*u
%

function Lu=fourier_Ku_1d(v,pde,domain)

kx=domain.k(:,1);

a=pde.a;
b=pde.b;
c=pde.c;

Ku=-ifft(-kx.^2.*fft(v))+a.*ifft(1i.*kx.*fft(v))+b.*v+c.*v.^2;

Lu=real(ifft(-kx.^2.*fft(Ku)));

end