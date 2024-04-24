% Performs multiplication L*u
%
% L is the spectral operator for the given PDE in the form
%
% (au_xxxx+bu_xx+cu)_xx+du_yy=0

function f=Lu_kp(v,pde,domain)

k=domain.k;

f=real(ifft((-pde.a.*k{1}.^6+pde.b.*k{1}.^4).*fft(v))-ifft(k{1}.^2.*fft(pde.c.*v))-transpose(ifft(pde.d.*k{2}.^2.*fft(transpose(v)))));

end