function [wint] = fentint(z,w,zint)

ci = complex(.0,1.);
N = length(z);
h = 2*pi/N;
K = ci*h*[0:N/2-1 0 -N/2+1:-1];
zd = ifft(K'.*fft(z));

[m,n] = size(zint);

for i = 1:m
	for j = 1:n
		wint(i,j) = 1/(2*pi*ci)*sum(w.*zd./(z-zint(i,j)+eps));
	end
end

end