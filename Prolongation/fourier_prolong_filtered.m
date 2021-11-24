% Fourier Prolongation operator on x direction (first component)
% assuming fine grid is 2x coarse grid points in x
%
% Assumed filtered for highest mode
%
% Inputs:
% vc - coarse grid (Nx * Ny * ... size)
%
% Outputs:
% vc - coarse grid prolongation (Nx*2 * Ny * ... size)


function vf=fourier_prolong_filtered(vc)

N=size(vc);
N(1)=N(1)/2;

% FFT + shift
vc_hat=fftshift(fft(vc),1);

% Filter high frequency
vc_hat=[zeros(N);vc_hat;zeros(N)];

% shift + iFFT
vf=2*real(ifft(fftshift(vc_hat,1)));

end