% Performs preconditioning of M=c-dxx(adxx+b)-dyy, c=??
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

function v=fourier_Ku_pre_2d(~,pde,domain,~)

    kx=domain.k{1};
    ky=domain.k{2};
    
    [KX,KY]=ndgrid(kx,ky);

    c=1;
    
    v_hat=fft2(pde.f)./(-KY.^2-KX.^2.*(pde.a.*KX.^2+pde.b)-c);
    v_hat(1,1)=0;
    v=real(ifft2(v_hat));
    
end