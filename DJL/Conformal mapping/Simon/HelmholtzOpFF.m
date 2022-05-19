function Lspec = HelmholtzOpFC(u,K,L,G,N,M)
    G = reshape(G,[N,M]);
% construct spectral operator for the LHS
    um = reshape(u,[N,M]);    
	Dx2 = real(ifft(K.^2.*fft(um)));
	ut = um';
	Dy2 = real(ifft(L.^2.*fft(ut)));
	Lspec = reshape(G,[N*M,1]).*u-reshape(Dx2+Dy2',[N*M,1]);
end
