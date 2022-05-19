function u = HelmholtzOpFD(G,F,xa,ya,K,L,nc,mc)
    G = reshape(G,[nc,mc]);
% construct FD approximation to operator for LHS
    ex = ones(nc,1);
    K1 = imag(K(2));
    D2fdx = (nc*K1/(2*pi))^2*spdiags([ex ex -2*ex ex ex], ...
        [-nc+1 -1 0 1 nc-1],nc,nc);
	dy = ya-circshift(ya,1);
    dy = reshape(dy,[mc,1]);
    dym = .5*(dy+circshift(dy,-1));
    dm1 = 1./(dy.*dym);
    dp1 = 1./(circshift(dy,-1).*dym);
    D2fdy = spdiags([circshift(dm1,-1),-(dm1+dp1),circshift(dp1,1)], ...
        [-1,0,1],mc,mc);
    D2fdy(1,:) = [1. zeros(1,mc-1)];
    D2fdy(mc,:) = [zeros(1,mc-1) 1.];
    G(:,1:mc-1:mc) = .0;
    M = -kron(eye(mc),D2fdx)-kron(D2fdy,eye(nc)) ...
        +spdiags(reshape(G,[mc*nc,1]),0,nc*mc,nc*mc);
	u = M\F;
end
