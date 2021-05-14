% 2D Fourier Restriction operator in the x direction only
% assuming fine grid is 2x coarse grid points
%
% Assumes using midpoint method with adjustment
%
% Inputs:
% vf - fine grid (Nx * Ny size)
%
% Outputs:
% vc - fine grid restriction (Nx/2 * Ny size)

function vc=fourier_restrict_1d_mid(vf)

[Nx,Ny]=size(vf);

% FFT + shift
v_hat=fftshift(fft(vf),1);

% Filter high frequency
v_hat=v_hat(Nx/4+1:3*Nx/4,:);
v_hat(1,:)=v_hat(1,:)*2; % Midpoint adjustment

% shift + iFFT
vc=0.5*real(ifft(fftshift(v_hat,1)));

end
