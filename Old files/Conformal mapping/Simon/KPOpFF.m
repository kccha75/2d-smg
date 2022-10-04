function Lspec = KPOpFC(u,G,X,Y,K,L,N,M)
    G = reshape(G,[N,M]);
% needs a parameter for Yop (sigsq)
% if function is passed to MG, this could be embedded
% construct spectral operator for the LHS
%
	global sigsq
%
    um = reshape(u,[N,M]);    
	Dx2 = real(ifft(K.^2.*fft(um)));
	wm = -G.*um+Dx2;
	Xop = real(ifft(K.^2.*fft(wm)));
	ut = um';
	Dy2 = real(ifft(L.^2.*fft(ut)));
	Lspec = Xop+sigsq*Dy2';
	Lspec = reshape(Lspec,[N*M,1]);
end
