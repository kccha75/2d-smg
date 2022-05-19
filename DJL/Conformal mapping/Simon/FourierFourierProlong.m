function [fr] = FourierFourierProlong(fc,N,M)
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
% Fourier
%
	ft = ff';
    gr = fft(ft);
    gp = [gr(1:M/2-1,:); zeros(M,2*N); gr(M/2:M,:)];
    ff = real(2*ifft(gp));
	fr = ff';
%
	fr = reshape(fr,[2*N*2*M,1]);
end
