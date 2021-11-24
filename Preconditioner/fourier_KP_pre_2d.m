% Performs preconditioning of M=cc-dxx(adxx+b)-dyy
% for a given constant cc using Fourier Spectral methods
%
% Inputs:
% v - best estimate
% pde.a
% pde.b
% pde.d
% pde.f
% domain.k
%
% Ouputs:
% v - best estimate after preconditioner

function v=fourier_KP_pre_2d(~,pde,domain,~)

    kx=domain.k{1};
    ky=domain.k{2};
    
    [KX,KY]=ndgrid(kx,ky);

    cc=0.0001;
    
    v=real(ifft2(fft2(pde.f)./(cc+pde.d.*KY.^2-KX.^2.*(pde.a.*KX.^2-pde.b))));

end