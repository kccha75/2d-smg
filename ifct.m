% Inverse Fast Chebyshev transform in the x direction only (first component)
%
% Inputs:
% v_hat - Chebyshev space (Nx * Ny * ... size)
%
% Ouputs:
% v_hat - real space (Nx * Ny * ... size)
%

function v = ifct(v_hat)

N=size(v_hat);
M=N(1)-1;

v_hat=v_hat(:,:);
v_hat = M*[v_hat(1,:)*2; v_hat(2:M,:); v_hat(end,:)*2];
v = ifft([v_hat; flipud(v_hat(2:end-1,:))]);
v = v(1:M+1,:);

% Reshape to correct size (required for 3d or higher)
v=reshape(v,N);

end