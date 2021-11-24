% Chebyshev differentiation in the x direction only (first component)
%
% Inputs:
% v_hat - v in Chebyshev space (Nx * Ny * ... size)
% n - number of times to differentiate
%
% Ouputs:
% dv_hat - derivative in Chebyshev space (Nx * Ny * ... size)
%

function dv_hat = chebdiff(v_hat,n)

N=size(v_hat);
M=N(1)-1;
k=linspace(0,M,M+1)';
    
v_hat=v_hat(:,:);
dv_hat=zeros(size(v_hat));


for i=1:n
% Cheb differentiation (see Canuto)
dv_hat(M+1,:) = .0;
dv_hat(M,:) = 2*k(M+1)*v_hat(M+1,:);
    
for j = M-1:-1:2
	dv_hat(j,:) = dv_hat(j+2,:)+2*k(j+1)*v_hat(j+1,:);
end
dv_hat(1,:) = .5*dv_hat(3,:)+v_hat(2,:);
    
v_hat=dv_hat;

end

% Reshape to correct size (required for 3d or higher)
dv_hat=reshape(dv_hat,N);

end