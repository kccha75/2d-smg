% Performs multiplication L*u using Fourier collocation methods
% L is the spectral operator for the given PDE in the form
% (-u_xxxx+au_xx+bu+cu^2)_xx-u_yy using
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

function Lu=fourier_5KPu_2d_adjoint(v,pde,domain)

kx=domain.k{1};
ky=domain.k{2};

a=pde.a;
b=pde.b;
c=pde.c;


Luy=ifft(-ky.^4.*fft(v'))+a.*ifft(-ky.^2.*fft(v'))+b'.*v';

if c~=0
    Luy=Luy+c.*v.^2;
end

Lux=ifft(kx.^2.*fft(v));

Lu=real(Lux+Luy');

end