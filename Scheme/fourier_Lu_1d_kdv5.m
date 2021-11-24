% Performs multiplication L*u using Fourier collocation methods
% L is the spectral operator for the given PDE in the form
% -u_xxxx+au_xx+bu+cu^2 using
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

function Lu=fourier_Lu_1d_kdv5(v,pde,domain)

kx=domain.k{1};

a=pde.a;
b=pde.b;
c=pde.c;

Lu=real(-ifft(kx.^4.*fft(v))+a.*ifft(-kx.^2.*fft(v))+b.*v);

% Nonlinear term
if c~=0
    Lu=Lu+c.*v.^2;
end

end