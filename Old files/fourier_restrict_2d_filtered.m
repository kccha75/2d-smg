% 2D Fourier Restriction operator assuming fine grid is 2x
% coarse grid points
%
% Assumed filtered for highest mode
%
% Inputs:
% vf - fine grid
%
% Outputs:
% vc - fine grid restriction

function vc=fourier_restrict_2d_filtered(vf)

[Nx,Ny]=size(vf);

% FFT + shift
f_hat=fftshift(fft2(vf));

% Filter high frequency
f_hat=f_hat(Nx/4+1:3*Nx/4,Ny/4+1:3*Ny/4);

% shift + iFFT
vc=0.25*real(ifft2(fftshift(f_hat)));


end
