function [fr] = FourierChebProlong(fc,N,M)
% prolongation operator
%
	fc = reshape(fc,[N,M]);
%
% Fourier
%
    gr = fft(fc);
    gp = [gr(1:N/2-1,:); zeros(N,M); gr(N/2:N,:)];
    ff = real(2*ifft(gp));
%
% Chebyshev
%
	ft = ff';
    f1 = [ChebForward(ft); zeros(M-1,2*N)];
    ff = sqrt(2)*ChebBackward(f1);
	fr = ff';
%
	fr = reshape(fr,[2*N*(2*M-1),1]);
end
