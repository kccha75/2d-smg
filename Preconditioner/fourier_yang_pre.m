% Performs preconditioning of M=c-nabla^2 (see Yang Nonlinear Waves)
%
% Inputs:
% v - best estimate
% pde.f
% scheme.kx
% scheme.ky
%
% Ouputs:
% v - best estimate after preconditioner

function v=fourier_yang_pre(v,pde,scheme,option)

    c=3;
    [KX,KY]=ndgrid(scheme.kx,scheme.ky);
    K2=KX.^2+KY.^2;
    v=ifft2(fft2(pde.f)./-(c+K2));
    
end