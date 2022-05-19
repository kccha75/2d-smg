function Lspec = KPNLOpFF(u,V,K,L,kappa,N,M)
% if function is passed to MG, this could be embedded
% construct spectral operator for the LHS
%
    um = reshape(u,[N,M]);    
    Vm = reshape(V,[N,M]);    
	Dx2 = real(ifft(K.^2.*fft(um)));	
	Xop = V.*um+kappa*um.^3-Dx2;
	ut = um';
	Dy2 = real(ifft(L.^2.*fft(ut)));
	Lspec = Xop-Dy2';
	Lspec = reshape(Lspec,[N*M,1]);
end
