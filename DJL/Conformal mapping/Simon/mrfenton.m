function [w Ck] = mrfenton(z,zd,zk,phi,phid,psi,psid,iflag)
%
N = length(z);
%
for m = 1:N;
	Omega(:,m) = zd./(z-z(m)+eps);
	Omega(m,m) = .0;
end
%
alpha = real(Omega);beta = imag(Omega);
%
OmegaC = zd./(z-zk);
SO = sum(OmegaC);
%
Ak = real(SO);Bk = imag(SO);
%
alphaC = real(OmegaC);betaC = imag(OmegaC);
%
vp = phi;
up = psi;
wp = complex(up,vp);
err = 1.;
%
while err > 1.e-6
%
for m = 1:N;
  if (iflag(m) == 1)
     phi(m) = (psid(m)+sum(alpha(:,m).*(psi-psi(m))+beta(:,m).*phi))...
         /sum(beta(:,m));
  elseif (iflag(m) == 0)
     psi(m) = (-phid(m)-sum(alpha(:,m).*(phi-phi(m))-beta(:,m).*psi))...
         /sum(beta(:,m)); 
  end
  Ck = sum((Ak*betaC-Bk*alphaC).*phi+(Ak*alphaC+Bk*betaC).*psi)...
      ./(Ak.^2 + Bk.^2);
  psi(N/2+1:3*N/4) = Ck;
end
%
w = complex(phi,psi);
err = norm(w-wp,1);
% if err > 1.e2
%     disp('error is too big. Run the program again')
%     return
% end
disp(err);
wp = w;
end
w = complex(phi,psi);
end
