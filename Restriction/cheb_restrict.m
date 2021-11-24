% Chebyshev Restriction operator in the x direction only (first component)
% assuming fine grid is 2x coarse grid points
%
% Inputs:
% vf - fine grid (Nx * Ny * ... size)
%
% Outputs:
% vc - coarse grid (ceil(Nx/2+1) * Ny * ... size)

function vc=cheb_restrict(vf)

N=size(vf);

% Chebyshev transform
v_hat=fct(vf);

% For adjoint operator adjustment
v_hat(ceil(N(1)/2),:)=v_hat(ceil(N(1)/2),:)/2;

% Remove high frequency components and invert
vc=ifct(v_hat(1:ceil(N(1)/2),:));

% Reshape to correct size (required for 3d or higher)
N(1)=ceil(N(1)/2);
vc=reshape(vc,N);

end
