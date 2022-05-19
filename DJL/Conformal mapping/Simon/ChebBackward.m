function [fc] = ChebBackward(ff)
	[N,~] = size(ff);
	fe = [2*ff(1,:); ff(2:N-1,:); 2*ff(N,:); ff(N-1:-1:2,:)];
    ft = real(fft(fe));
	fc = ft(1:N,:)/sqrt(2*(N-1));
end
