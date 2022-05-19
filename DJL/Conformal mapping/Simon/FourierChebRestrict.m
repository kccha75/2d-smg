function [fr] = FourierChebRestrict(ff,N,M)
% restriction operator
%
	ff = reshape(ff,[N,M]);
	Nc = N/2;
	Mc = (M-1)/2+1;
%
% Fourier
%
    gt = fft(ff);
    gr = [gt(1:N/4-1,:); gt(3*N/4:N,:)];
    fc = .5*ifft(gr);
%
% Chebyshev
%
	fc = fc';
    f1 = ChebForward(fc);
    fc = ChebBackward(f1(1:Mc,:))/sqrt(2);
    fr = fc';
%
    fr = reshape(fr,[Nc*Mc,1]);
end
