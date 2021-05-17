% Performs preconditioning of M=c-nabla^2 (see Yang Nonlinear Waves)
%
% Inputs:
% v - best estimate
% pde.f
% domain.k
%
% Ouputs:
% v - best estimate after preconditioner

function v=fourier_yang_pre_2d(~,pde,domain,~)

    kx=domain.k(:,1);
    ky=domain.k(:,2);

    c=3;
    [KX,KY]=ndgrid(kx,ky);
    K2=KX.^2+KY.^2;
    v=ifft2(fft2(pde.f)./-(c+K2));
    
end