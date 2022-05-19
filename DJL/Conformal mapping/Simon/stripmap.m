%
% add:
% documentation
% calculation of Jacobian = x_u^2+y_v^2
%
NI = 10;
%
N = 512;
M = 128;
H0 = 1;
h = @(x) .1*(1-tanh((x-.5)/.1).^2)+.1*(1-tanh((x+.5)/.1).^2);
% h = @(x) .5*(1-tanh(x/.25).^2);
% h = @(x) .25*exp(-(x/.25).^2);

u = -pi+2*pi*(0:N-1)/N;
K = 1i*[0:N/2-1 -N/2:-1];
Kinv = [.0 1./K(2:N)];
xi = .5*(1-cos(pi*(0:M)/M));
L1 = [0:M];
%
xb = u;
for nit = 1:NI
	yb = h(xb);
	L = H0-1/N*sum(yb);
	v = L*xi;
	LL = -2/L*L1;
	[U,V] = ndgrid(u,v);
	F = zeros(N,M+1);
	for ix = 1:N
		F(ix,:) = yb(ix)+(H0-yb(ix)).*V(ix,:)/L;
	end
	G = zeros(N,M+1);
	%
	% initialize the RHS
	%
	F = reshape(F,[N*(M+1),1]);
	nabF = -HelmholtzOpFC(F,G,u,v,K',LL,N,M+1);
	nabF = reshape(nabF,[N,M+1]);
	nabF(:,1:M:M+1) = .0;
	%
	[Y,R] = SpecMG(1,u,v,G,nabF,@HelmholtzOpFD,@HelmholtzOpFC);
	Yn = Y+reshape(F,[N,M+1]);
	Yt = Yn';
	Dy = ChebBackward(ChebDeriv(ChebForward(Yt),LL,1));
	Dx = Dy(1,:);
	Ix = mean(Dx)*(u+pi)-pi+real(ifft(Kinv.*fft(Dx)));
%    Ix = 2*pi/N*(cumsum(Dx)-Dx(1))-pi;
    fprintf('Change in x = %g \n',norm(Ix-xb));
	if norm(Ix-xb) < 1.e-8
		break;
	end
	xb = Ix;
end
%
Dx = Dy(M+1,:);
% fo = mean(Dx);
xt = mean(Dx)*(u+pi)-pi+real(ifft(Kinv.*fft(Dx)));
% xt = 2*pi/N*(cumsum(Dx)-Dx(1))-pi;
xbp = xb-u;
xtp = xt-u;
F = zeros(N,M+1);
for ix = 1:N
	F(ix,:) = xbp(ix)+(xtp(ix)-xbp(ix)).*V(ix,:)/L;
end
F = reshape(F,[N*(M+1),1]);
nabF = -HelmholtzOpFC(F,G,u,v,K',LL,N,M+1);
nabF = reshape(nabF,[N,M+1]);
nabF(:,1:M:M+1) = .0;
% X = HelmholtzOpFD(G,reshape(nabF,[N*(M+1),1]),u,v,K',LL,N,M+1);
[X,R] = SpecMG(1,u,v,G,nabF,@HelmholtzOpFD,@HelmholtzOpFC);
Xn = X+reshape(F,[N,(M+1)]);
XF = Xn+U;
YF = Yn;
%
clf
subplot(211)
contour(XF,YF,U,40); 
hold on 
contour(XF,YF,V,40)
xlabel('x')
ylabel('y')
plot(xb,h(xb))
subplot(212)
contour(U,V,XF,40); 
hold on 
contour(U,V,YF,40)
xlabel('u')
ylabel('v')
