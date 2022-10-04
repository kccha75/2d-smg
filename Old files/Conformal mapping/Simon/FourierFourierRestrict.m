function [fr] = FourierFourierRestrict(ff,N,M)
% restriction operator
%
	ff = reshape(ff,[N,M]);
%
% Fourier
%
    gt = fft(ff);
    gr = [gt(1:N/4-1,:); gt(3*N/4:N,:)];
    fc = .5*ifft(gr);
%
% Fourier
%
	fc = fc';
    gt = fft(fc);
    gr = [gt(1:M/4-1,:); gt(3*M/4:M,:)];
    fr = .5*ifft(gr);
    fr = fr';
%
    fr = reshape(fr,[N/2*M/2,1]);
end
