% Chebyshev Prolongation operator on x direction (first component)
% assuming fine grid is 2x+1 coarse grid points in x
%
% Inputs:
% vc - coarse grid (Nx * Ny * ... size)
%
% Outputs:
% vc - coarse grid prolongation (Nx*2-1 * Ny * ... size)

function vf=cheb_prolong(vc)

N=size(vc);
N(1)=N(1)-1;

% Chebyshev transform
vc_hat=fct(vc);

% Filter high frequency and invert
vf=ifct([vc_hat;zeros(N)]);

% Apply BCs
vf(1,:)=vc(1,:);
vf(end,:)=vc(end,:);

end
