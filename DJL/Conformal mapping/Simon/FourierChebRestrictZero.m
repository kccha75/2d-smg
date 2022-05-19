function [fr] = FourierChebRestrictZero(ff,N,M)
% restriction operator to account for dirichlet boundary conditions
%
	ff = reshape(ff,[N,M]);
	M1 = M-1;
	y = cos(pi*(0:M1)/M1)';
%
% Fourier
%
    gt = fft(ff);
    gr = [gt(1:N/4-1,:); gt(3*N/4:N,:)];
    fc = real(.5*ifft(gr));
    fc = real(fc);
%
% Chebyshev
%
	fc = fc';
    ind = 2:M1;
    dop1 = spdiags(2*(-1).^(ind)./(1-y(ind)),0,M1-1,M1-1);
    dop2 = spdiags(2*(-1).^(M1-ind)./(1+y(ind)),0,M1-1,M1-1);
    a = sum(dop1*fc(ind,:),1);
    b = sum(dop2*fc(ind,:),1);
    uo = ((2*M1^2+1)/6*a-.5*(-1)^M1*b)/(((2*M1^2+1)/6)^2-.25);
    un = ((2*M1^2+1)/6*b-.5*(-1)^M1*a)/(((2*M1^2+1)/6)^2-.25);
    f1 = ChebForward([uo; fc(2:M1,:); un]);
    fc = ChebBackward(f1(1:M1/2+1,:))/sqrt(2);
    fc(1:M1/2:M1/2+1,:) = .0;
    fr = fc';
%
    fr = reshape(fr,[N/2*(M1/2+1),1]);
end
