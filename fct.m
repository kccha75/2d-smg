% Fast Chebyshev transform in the x direction only (first component)
%
% Inputs:
% v - real space (Nx * Ny * ... size)
%
% Ouputs:
% v_hat - Chebyshev space (Nx * Ny * ... size)
%

function v_hat=fct(v)

N=size(v);
M=N(1)-1;

v=v(:,:);
v = [v; flipud(v(2:M,:))];
v_hat = real(fft(v))/M;
v_hat = [v_hat(1,:)/2; v_hat(2:M,:); v_hat(M+1,:)/2];

% Reshape to correct size (required for 3d or higher)
v_hat=reshape(v_hat,N);

end