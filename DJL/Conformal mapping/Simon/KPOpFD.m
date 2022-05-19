function xout = KPOpFD(G,F,x,y,K,L,nc,mc)
%
% needs a parameter for Yop (sigsq)
% if function is passed to MG, this could be embedded
%
	global sigsq
%
    G = reshape(G,[nc,mc]);
% construct FD approximation to operator for LHS
    ex = ones(nc,1);
    K1 = imag(K(2));
    D2fdx = (nc*K1/(2*pi))^2*spdiags([ex ex -2*ex ex ex], ...
        [-nc+1 -1 0 1 nc-1],nc,nc);
	LD2fx = kron(eye(mc),D2fdx);
	Xop = LD2fx*(-spdiags(reshape(G,[mc*nc,1]),0,nc*mc,nc*mc)+LD2fx);
	ey = ones(mc,1);
    L1 = imag(L(2));
    D2fdy = (mc*L1/(2*pi))^2*spdiags([ey ey -2*ey ey ey], ...
        [-mc+1 -1 0 1 mc-1],mc,mc);
	Yop = kron(D2fdy,eye(nc));
    M = Xop+sigsq*Yop;
    xout = M\F;
end
