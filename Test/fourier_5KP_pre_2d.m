% Performs preconditioning of M=c-dxx(dxxxx +2*dxx-mu), c=0.0001
% (see Yang paper)
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

function v=fourier_5KP_pre_2d(~,pde,domain,~)

    kx=domain.k{1};

    c=0.0001;
    mu=-1.2;
    
    v=ifft(fft(pde.f)./(-kx.^2.*(kx.^4-2*kx.^2-mu)-c));
    
end