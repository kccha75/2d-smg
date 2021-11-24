% Performs preconditioning of M=c-dxx(dxxxx +2*dxx-mu)-dyy, c=0.0001
% (see Yang paper)
%
% TESTING ADDING -dyy term
%
% Note: mu is set to -1.2
%
% Inputs:
% v - best estimate
% pde.f
% domain.k
%
% Ouputs:
% v - best estimate after preconditioner

function v=fourier_5KP_pre_2d_test(~,pde,domain,~)

    kx=domain.k{1};
    ky=domain.k{2};
    
    [KX,KY]=ndgrid(kx,ky);

    c=0.0001;
    mu=-1.2;
    
    v_hat=fft2(pde.f)./(-KY.^2-KX.^2.*(KX.^4-2*KX.^2-mu)-c);
    
%     % Mean 0
    v_hat(1,1)=0;
    v=ifft2(v_hat);
    
end