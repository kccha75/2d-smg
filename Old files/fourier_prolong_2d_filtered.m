% 2D Fourier Prolongation operator assuming fine grid is 2x
% coarse grid points
%
% Inputs:
% vc - coarse grid
%
% Outputs:
% vf - coarse grid prolongation

function vf=fourier_prolong_2d_filtered(vc)

[Nx,Ny]=size(vc);

% FFT + shift
vc_hat=fftshift(fft2(vc));

% Pad high frequency with 0's
vc_hat=[zeros(Nx/2,2*Ny);zeros(Nx,Ny/2) vc_hat zeros(Nx,Ny/2); zeros(Nx/2,2*Ny)];

% shift + iFFT
vf=4*real(ifft2(fftshift(vc_hat)));

end