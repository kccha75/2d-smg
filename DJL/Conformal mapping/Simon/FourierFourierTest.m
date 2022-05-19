function FourierChebTest(tcase)

nnewton = 10;

clf
switch tcase
	case 0
%
		kappa = -1;
		P0 = 2.35;
        M = 128;
        N = 256;
        Lx = 10*pi;
        Ly = 2*pi;
        x = Lx*(0:N-1)/N;
        y = Ly*(0:M-1)/M;
        ci = complex(0,1);
		K1 = 2*pi/Lx*ci*[0:N/2-1 -N/2:-1]';
		L1 = 2*pi/Ly*ci*[0:M/2-1 -M/2:-1]';;
		weights = Lx*Ly/(N*M);
        [X,Y] = ndgrid(x,y);
        xo = x(N/2);
%
% unforced KdV solitary wave
%
		pert = .5;
		noise = .0;
		kap = 5;
		Delta = -1;
		A0 = 1;
		A = A0*(1-tanh(kap*(X-xo-pert*cos(Y))).^2) ...
			+noise*randn(N,M);
		V = (1-tanh(kap*(X-xo-pert*cos(Y))).^2)-Delta;
		uf = zeros(N*M,1);
		rf = zeros(N*M,1);
		P = sqrt(sum(weights.*sum(A.^2)));
		A = P0*A/P;
		for inewton = 1:nnewton
			G = V+3*kappa*A.^2;
			F = -reshape(GPNLOpFF(A,V,K1,L1,kappa,N,M),[N,M]);
			rms = sqrt(sum(weights.*sum(F.^2)));
			[uf,rf] = SpecMG(0,x,y,G,F,@GPOpFD,@GPOpFF);
			A = A+reshape(uf,[N,M]);
			P = sqrt(sum(weights.*sum(A.^2)));
			fprintf('inewt = %g, rmse = %g, P = %g \n',inewton,rms,P);
			A = P0*A/P;
		end
        subplot(3,1,1); mesh(X,Y,reshape(uf,[N,M]));
        subplot(3,1,2); mesh(X,Y,A+reshape(uf,[N,M]));
        subplot(3,1,3); mesh(X,Y,reshape(rf,[N,M]));
    otherwise
        fprintf('Nada. \n')
end