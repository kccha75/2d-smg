% 2D Fourier Prolongation operator with x and y sweeps
% assuming fine grid is 2x coarse grid points
%
% Assumes using midpoint method with adjustment
%
% Inputs:
% vc - coarse grid (Nx * Ny size)
%
% Outputs:
% vc - coarse grid prolongation (Nx*2 * Ny*2 size)

function vf=fourier_prolong_2d_mid(vc)

% FFT + shift in x
vf=fourier_prolong_1d_x(vc);

% FFT + shift in y
vf=fourier_prolong_1d_x(vf');

vf=vf';

end

% -------------------------------------------------------------------------
% Subroutine 1d prolongation
% -------------------------------------------------------------------------

function vf=fourier_prolong_1d_x(vc)

[Nx,Ny]=size(vc);

% FFT + shift
vc_hat=fftshift(fft(vc),1);

% Filter high frequency
vc_hat(1,:)=vc_hat(1,:)/2; % Midpoint adjustment
vc_hat=[zeros(Nx/2,Ny);vc_hat;zeros(Nx/2,Ny)];

% shift + iFFT
vf=2*real(ifft(fftshift(vc_hat,1)));

end

