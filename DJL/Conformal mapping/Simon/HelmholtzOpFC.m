function Lspec = HelmholtzOpFC(u,G,xa,ya,K,L,N,M)
    G = reshape(G,[N,M]);
% construct spectral operator for the LHS
    um = reshape(u,[N,M]);    
	Dx2 = real(ifft(K.^2.*fft(um)));
	ut = um';
    Dy2 = ChebBackward(ChebDeriv(ChebForward(ut),L,2));
	Lspec = reshape(G,[N*M,1]).*u-reshape(Dx2+Dy2',[N*M,1]);
    Lspec(1:N) = u(1:N);
    Lspec((M-1)*N+1:M*N) = u((M-1)*N+1:M*N);
end
