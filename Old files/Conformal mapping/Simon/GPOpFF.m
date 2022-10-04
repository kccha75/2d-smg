function Lspec = GPOpFC(u,G,X,Y,K,L,N,M)
    G = reshape(G,[N,M]);
% needs a parameter for Yop (sigsq)
% if function is passed to MG, this could be embedded
% construct spectral operator for the LHS
%
    um = reshape(u,[N,M]);    
	Dx2 = -real(ifft(K.^2.*fft(um)));
	Xop = G.*um+Dx2;
	ut = um';
	Dy2 = -real(ifft(L.^2.*fft(ut)));
	Lspec = Xop+Dy2';
	Lspec = reshape(Lspec,[N*M,1]);
end
