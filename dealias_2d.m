% 2D Dealias fourier u*v using 2/3 rule (padding)
% 
% Inputs:
% u 
% v (must be same size as u)
%
% Outputs:
% w - dealiased multiplication

function w=dealias_2d(u,v)

if length(u) ==1 || length(v) == 1
    w=u.*v;
    return
end

[N,M]=size(u);

% new padding size
K1=ceil(3/2*N);
K2=ceil(3/2*M);

% index for padding
padindex1 = [1:N/2,K1-N/2+1:K1];
padindex2 = [1:M/2,K2-M/2+1:K2];

% Check if vectors are equal
if isequal(u,v)==0

    u_hat=fft2(u);
    v_hat=fft2(v);
    
    % 2/3 rules, pad with 0's
    u_hat_pad=zeros(K1,K2);
    v_hat_pad=zeros(K1,K2);
    
    u_hat_pad(padindex1,padindex2)=u_hat;
    v_hat_pad(padindex1,padindex2)=v_hat;
    
    % Multiply
    w=real(ifft2(u_hat_pad).*ifft2(v_hat_pad));
    
else
    
    u_hat=fft2(u);
    
    % 2/3 rules, pad with 0's
    u_hat_pad=zeros(K1,K2);
    
    u_hat_pad(padindex1,padindex2)=u_hat;
    
    % Multiply
    w=real(ifft2(u_hat_pad).^2);    

end

% remove pads and -N/2 mode
w_hat=fft2(w)*K1/N*K2/M;

w_hat=w_hat(padindex1,padindex2);
w_hat(N/2+1,:)=0;
w_hat(:,M/2+1)=0;

w=ifft2(w_hat);

end