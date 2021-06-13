% Performs preconditioning for 5th order KdV (see Yang's paper)
% in the form M=-d_xxxx+ad_xx+b
%
% Inputs:
% v - best estimate
% pde.f
% domain.k
%
% Ouputs:
% v - best estimate after preconditioner

function v=fourier_yang_pre_1d_kdv5(~,pde,domain,~)

    k=domain.k{1};
    
    v=real(ifft(fft(pde.f)./(-k.^4-k.^2+pde.b)));
    
end