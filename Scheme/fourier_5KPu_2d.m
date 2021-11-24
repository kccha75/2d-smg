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

function Lu=fourier_5KPu_2d(v,pde,domain)

kx=domain.k{1};
ky=domain.k{2};

a=pde.a;
b=pde.b;
c=pde.c;

% Bracket terms
Lux=ifft(-kx.^4.*fft(v))+a.*ifft(-kx.^2.*fft(v))+b.*v;

% Nonlinear term
if c~=0
    Lux=Lux+c.*v.^2;
end

% Derivate of brackets
Lux=ifft(-kx.^2.*fft(Lux));

% y derivatives
Luy=ifft(ky.^2.*fft(v'));

Lu=real(Lux+Luy');

end