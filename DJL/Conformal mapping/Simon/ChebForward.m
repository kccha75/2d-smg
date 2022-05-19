function [fc] = ChebForward(ff)
	[N,~] = size(ff);
    fe = [ff; ff(N-1:-1:2,:)]/sqrt(2*(N-1));
    ft = real(fft(fe));
    fc = [.5*ft(1,:); ft(2:N-1,:); .5*ft(N,:)];
end
