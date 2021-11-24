% 2D Cheb-Fourier Prolongation operator with x and y sweeps
% assuming fine grid is 2x coarse grid points
%
% Assumes using filtered method with adjustment
%
% Inputs:
% vc - coarse grid (Nx * Ny size)
%
% Outputs:
% vc - coarse grid prolongation (Nx*2 * Ny*2 size)

function vf=cheb_fourier_prolong_filtered(vc)

% Cheb prolong in x
vf=cheb_prolong(vc);

% Fourier prolong in y
vf=fourier_prolong_filtered(vf');

vf=vf';

end