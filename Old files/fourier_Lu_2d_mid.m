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
% domain.k
% domain.N
%
% Ouputs:
% A - L*u
%
function Lu=fourier_Lu_2d_mid(u,pde,domain)

kx=domain.k{1};
ky=domain.k{2};
Nx=domain.N(1);
Ny=domain.N(2);

a=pde.a;
b=pde.b;
c=pde.c;

% Constant coefficient check
if length(a)==1
    a_mid=a;
else   
    % Shifted midpoint at i+1/2
    a_mid=real(ifft(fft(a).*exp(pi*1i*kx/Nx)));
end

if length(b)==1
    b_mid=b;
else   
    b_mid=real(ifft(fft(b').*exp(pi*1i*ky/Ny)));
end

% -(au_x)_x term
Lu1=-ifft(1i*kx.*fft(a_mid.*ifft(1i*kx.*fft(u).*exp(pi*1i*kx/Nx))).*exp(-pi*1i*kx/Nx));

% -(bu_y)_y term
Lu2=-ifft(1i*ky.*fft(b_mid.*ifft(1i*ky.*fft(u').*exp(pi*1i*ky/Ny))).*exp(-pi*1i*ky/Ny));

Lu=real(Lu1+Lu2'+c.*u);

end