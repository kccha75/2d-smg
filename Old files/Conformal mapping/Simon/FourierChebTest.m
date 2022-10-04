function FourierChebTest(tcase)
clf
switch tcase
    case 0
        M = 32;
        K = 0:M;
        y = cos(pi*(0:M)/M);
        theta = tanh(y);
        ff = zeros(2,M+1);
        ff(1,:) = theta;
        ff(2,:) = theta.^2;
        ft = ff';
        ftt = ChebBackward(ChebDeriv(ChebForward(ft),K,0));
        plot(y,ftt(:,1)'-ff(1,:),y,ftt(:,2)'-ff(2,:));
        legend('ft1','ft2');
	case 1
        M = 32;
        K = 0:M;
        y = cos(pi*(0:M)/M);
        theta = tanh(y);
        ff = zeros(2,M+1);
        ff(1,:) = theta;
        ff(2,:) = theta.^2;
        fd(1,:) = -2*theta.*(1.-theta.^2);
        fd(2,:) = (2.-6*theta.^2).*(1-theta.^2);
        ft = ff';
        Deriv = ChebBackward(ChebDeriv(ChebForward(ft),K,2));
        plot(y,Deriv(:,1)'-fd(1,:),y,Deriv(:,2)'-fd(2,:));
        legend('fd1','fd2');
    case 2
        M = 32;
        N = 2;
        x = 2*pi*(0:M-1)/M;
       	K = complex(0,1)*[0:M/2-1 -M/2:-1]';
        cs = cos(x);
        sn = sin(x);
        ff = zeros(M,N);
        fd = ff;
        ff(:,1) = cs.^3+ sn.^2;
        ff(:,2) = cs.^4.*sn;
        fd(:,1) = -3*cs.^2.*sn + 2*sn.*cs;
        fd(:,2) = -4.*cs.^3.*sn.^2+cs.^5;
        Dx = real(ifft(K.*fft(ff)));
        plot(x,Dx(:,1)-fd(:,1),x,Dx(:,2)-fd(:,2));
    case 3
        M = 32;
        N = 32;
        x = 2*pi*(0:N-1)/N;
       	K = complex(0,1)*[0:N/2-1 -N/2:-1]';
        L = 0:M;
        y = cos(pi*(0:M)/M);
        [X,Y] = ndgrid(x,y);
        cs = cos(X);
        sn = sin(X);
        th = tanh(Y);
        ff = (cs.^3+ sn.^2).*th;
        fr = FourierChebRestrict(ff);
        xc = 2*pi*(0:N/2-1)/(N/2);
        yc = cos(pi*(0:M/2)/(M/2));
        [Xc,Yc] = ndgrid(xc,yc);
        subplot(2,1,1); mesh(X,Y,ff)
        subplot(2,1,2); mesh(Xc,Yc,fr)
    case 4
        M = 32;
        N = 32;
        x = 2*pi*(0:N-1)/N;
       	K = complex(0,1)*[0:N/2-1 -N/2:-1]';
        L = 0:M;
        y = cos(pi*(0:M)/M);
        [X,Y] = ndgrid(x,y);
        cs = cos(X);
        sn = sin(X);
        th = tanh(Y);
        ff = (cs.^3+ sn.^2).*th;
        ff(:,1:M:M+1) = .0;
        fr = FourierChebRestrictZero(ff,y);
        xc = 2*pi*(0:N/2-1)/(N/2);
        yc = cos(pi*(0:M/2)/(M/2));
        [Xc,Yc] = ndgrid(xc,yc);
        subplot(2,1,1); mesh(X,Y,ff)
        subplot(2,1,2); mesh(Xc,Yc,fr)
     case 5
        M = 32;
        N = 32;
        x = 2*pi*(0:N-1)/N;
       	K = complex(0,1)*[0:N/2-1 -N/2:-1]';
        L = 0:M;
        y = cos(pi*(0:M)/M);
        [X,Y] = ndgrid(x,y);
        cs = cos(X);
        sn = sin(X);
        th = tanh(Y);
        ff = (cs.^3+ sn.^2).*th;
        fr = FourierChebProlong(ff);
        xc = 2*pi*(0:2*N-1)/(2*N);
        yc = cos(pi*(0:2*M)/(2*M));
        [Xc,Yc] = ndgrid(xc,yc);
        subplot(2,1,1); mesh(X,Y,ff)
        subplot(2,1,2); mesh(Xc,Yc,fr)   
     case 6
        M = 128;
        N = 128;
        x = 2*pi*(0:N-1)/N;
       	K = complex(0,1)*[0:N/2-1 -N/2:-1]';
        L = 0:M;
        y = cos(pi*(0:M)/M);
        dy = y-circshift(y,1);
        [X,Y] = ndgrid(x,y);
        G = ones(N,M+1);
        G = 100*Y;
        F = zeros(N,M+1);
        F = sin(2*X).*(1-tanh(Y).^2)+2;
        F(:,1:M:M+1) = .0;
        MM = HelmholtzOpFD(G,K,dy);
        uf = MM\reshape(F,[N*(M+1),1]);
        res = HelmholtzOpFC(uf,K,L,G)-reshape(F,[N*(M+1),1]);
        subplot(2,1,1); mesh(X,Y,reshape(uf,[N,M+1]));
        subplot(2,1,2); mesh(X,Y,reshape(res,[N,M+1]));
    	weights = 2*pi^2/(N*M)*sqrt(1-y.^2);
        rms = sqrt(sum(weights.*sum(res.^2,1)));
        disp(rms);
     case 7
        M = 128;
        N = 128;
        x = 2*pi*(0:N-1)/N;
        y = cos(pi*(0:M)/M);
        [X,Y] = ndgrid(x,y);
        G = ones(N,M+1);
        G = 100*Y;
        F = zeros(N,M+1);
        F = sin(2*X).*(1-tanh(Y).^2).^2+2;
        F(:,1:M:M+1) = .0;
        [uf,rf] = FourierChebMG(1,x,y,G,F,@HelmholtzOpFD,@HelmholtzOpFC);
        subplot(2,1,1); mesh(X,Y,reshape(uf,[N,M+1]));
        subplot(2,1,2); mesh(X,Y,reshape(rf,[N,M+1]));
    otherwise
        fprintf('Nada. \n')
end
end