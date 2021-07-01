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
function Lu=fourier_KPu_2d_dealias(v,pde,domain)

kx=domain.k{1};
ky=domain.k{2};

a=pde.a;
b=pde.b;
c=pde.c;
d=pde.d;

% inside bracket terms

% if length(v)==2^10
% Lx=a.*ifft(-kx.^2.*fft(v))+b.*v+c.*v.*v;
% else
Lx=a.*ifft(-kx.^2.*fft(v))+dealias_2d(b,v)+c.*dealias_2d(v,v);
% end
% x terms
Lu_x=real(ifft(-kx.^2.*fft(Lx)));

% y term
Lu_y=d.*real(ifft(-ky.^2.*fft(v')));

Lu=Lu_x+Lu_y';

end