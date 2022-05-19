%
% SMG function to solve
% L{u} = -nabla^2 u + g*u = f , (1)
% on periodic channel with homogeneous Dirichlet boundary conditions
%
% fourier is x direction [0,2*pi), and cheb is y direction [-1,1]
% assume for the moment that x and y are normalized
%
function [uf,rf] = FourierChebMG(x,y,g,f)
%
% test:
% Newton-multigrid,
% FAS,
% nonlinear wave problems (how to do autonomous systems),
% underdetermined systems, e.g., Poisson's equation
%
% todo:
% set up iterators and solvers at top of script
% set up parameters
% two-dimensional
% Chebyshev
% set up of problem
% want a structure which specifies sequence to follow
%
	[N,M] = size(g); % number of points in the finest domain
	if length(x) ~= N || length(y) ~= M || sum(size(f) ~= size(g))
		fprintf('Array sizes incompatible. \n');
	end	
%
% normalize variables so that zeta = [1,-1] and eta = Lx/Ly*[0,1),
% where Ly = .5*(y(M)-y(1)),
% y = y(1)+Ly*(1-zeta), x = abs(Ly)*eta
%
	Ly = (y(1)-y(M))/2.;
	Lx = N*(x(2)-x(1));
	K1 = 2*pi*Ly/Lx;
	f = Ly^2*f;
	g = Ly^2*g;
	ci = complex(0,1);
	zeta = cos(pi*(0:M-1)/(M-1));
	eta = Lx/Ly*(0:N-1)/N;
	weights = Lx*pi/(Ly*N*(M-1))*sqrt(1-zeta.^2);
%
% parameters for multigrid
%
	nlevels = 4; % number of levels
    nmgit = 5; % number of iterations
    tolerance = 1.e-10; % desired accuracy
%
%%%%%%%%%%%%%%%% setup %%%%%%%%%%%%%%%%%%%
%
    nxgrid = (N./2.^(0:nlevels-1))'; % array with number of points in x direction
    nygrid = ((M-1)./2.^(0:nlevels-1)+1)'; % array with number of points in y direction
    ngrid = nxgrid.*nygrid;
	dtlength = sum(ngrid); % length of dataset to store solution and residuals
	dxlength = sum(nxgrid); % length of dataset to x, dx and k
	dylength = sum(nygrid); % length of dataset to y and l
	strt = zeros(nlevels,3); % starting points for storage on each domain
	strt(1,:) = 1; % starting point on top domain
	for ilevel = 2:nlevels 
		strt(ilevel,1) = strt(ilevel-1,1)+ngrid(ilevel-1); % starting point on lower domains
		strt(ilevel,2) = strt(ilevel-1,2)+nxgrid(ilevel-1); % starting point on lower domains
		strt(ilevel,3) = strt(ilevel-1,3)+nygrid(ilevel-1); % starting point on lower domains
	end
	fin(:,1) = strt(:,1)+ngrid-1; % finishing point for each domain
	fin(:,2) = strt(:,2)+nxgrid-1; % finishing point for each domain
	fin(:,3) = strt(:,3)+nygrid-1; % finishing point for each domain
	u = zeros(dtlength,1); % solution
	fa = zeros(dtlength,1); % RHS
	ga = zeros(dtlength,1); % function on LHS of (1)
	r = zeros(dtlength,1); % residuals
	xa = zeros(dxlength,1); % array of x values
	K = zeros(dxlength,1); % array of x wavenumbers
	ya = zeros(dylength,1); % array of y values
	dy = zeros(dylength,1); % array of dy values
	L = zeros(dylength,1); % y wavenumbers
%
    clevel = 1;
    nc = nxgrid(clevel);
    mc = nygrid(clevel);
    cind = strt(clevel,1):fin(clevel,1);
    cindx = strt(clevel,2):fin(clevel,2);
    cindy = strt(clevel,3):fin(clevel,3);
	fa(cind) = reshape(f,[nc*mc,1]); % initialize RHS on top domain
	ga(cind) = reshape(g,[nc*mc,1]); % initialize LHS function on top domain
	xa(cindx) = eta;
	K(cindx) = K1*ci*[0:nc/2-1 -nc/2:-1]';
	ya(cindy) = zeta;
	dy(cindy) = ya(cindy)-circshift(ya(cindy),1);
	L(cindy) = [0:mc-1];
%
% do the same for the lower levels, except f is not passed between domains,
% rather the residual is passed
% use restriction to pass g down to the lower domains
%
	for clevel = 2:nlevels
		steps = 2^(clevel-1);
		nc = nxgrid(clevel);
		mc = nygrid(clevel);
		cind = strt(clevel,1):fin(clevel,1);
		cindx = strt(clevel,2):fin(clevel,2);
		cindy = strt(clevel,3):fin(clevel,3);
		ncp = nxgrid(clevel-1);
		mcp = nygrid(clevel-1);
		pind = strt(clevel-1,1):fin(clevel-1,1);
		pindx = strt(clevel-1,2):fin(clevel-1,2);
		pindy = strt(clevel-1,3):fin(clevel-1,3);
		ga(cind) = FourierChebRestrict(ga(pind),ncp,mcp);
		xa(cindx) = eta(1:steps:end);
		K(cindx) = K1*ci*[0:nc/2-1 -nc/2:-1]';
		ya(cindy) = zeta(1:steps:end);
		dy(cindy) = ya(cindy)-circshift(ya(cindy),1);
		L(cindy) = [0:mc-1];
    end
%
%%%%%%%%%%%%%%%% calculations %%%%%%%%%%%%%%%%%%%
%
    clevel = 1;
    nc = nxgrid(clevel);
    mc = nygrid(clevel);
	cind = strt(clevel,1):fin(clevel,1);
	cindx = strt(clevel,2):fin(clevel,2);
	cindy = strt(clevel,3):fin(clevel,3);
% setup periodic FD approximation M to L
    MM = HelmholtzOpFD(ga(cind),K(cindx),dy(cindy),nc,mc);
% % function to calculate residual at top domain
    Lop = @(uc) HelmholtzOpFC(uc,K(cindx),L(cindy),ga(cind),nc,mc);
% use FD to obtain an initial approximation to the solution
	xin = MM\fa(cind);
	r(cind) = Lop(xin)-fa(cind);
	rem = reshape(r(cind),[nc,mc]);
	rms = sqrt(sum(weights.*sum(rem.^2)));
	fprintf('rmse = %g \n',rms)
% use pCG as relaxation to improve solution and calculate residual
	[u(cind),r(cind)] = mpcg(Lop,fa(cind),xin,MM,2);
% do multigrid v-cycles
% should set this as a parameter and have a tolerance parameter to check after each 
% v-cycle
	for iv = 1:nmgit
% down the snakes
		for clevel = 2:nlevels
            nc = nxgrid(clevel);
            mc = nygrid(clevel);
            cind = strt(clevel,1):fin(clevel,1);
            cindx = strt(clevel,2):fin(clevel,2);
            cindy = strt(clevel,3):fin(clevel,3);
            ncp = nxgrid(clevel-1);
            mcp = nygrid(clevel-1);
            pind = strt(clevel-1,1):fin(clevel-1,1);
            pindy = strt(clevel-1,3):fin(clevel-1,3);
% rhs is residual from previous level
			fa(cind) = FourierChebRestrictZero(r(pind),ya(pindy),ncp,mcp);
% setup periodic FD approximation M to L
            MM = HelmholtzOpFD(ga(cind),K(cindx),dy(cindy),nc,mc);
% function to calculate residual at this level
			Lop = @(uc) HelmholtzOpFC(uc,K(cindx),L(cindy),ga(cind),nc,mc);
% use FD to obtain an initial approximation to the solution
			xin = MM\fa(cind);
% if we are at an intermediate level do pCG relaxation
% if we are at the bottom level use pCG to solve the problem
% should have parameter in pCG to iterate to solution
			if (clevel < nlevels)
				[u(cind),r(cind)] = mpcg(Lop,fa(cind),xin,MM,2);
			else
				[u(cind),r(cind)] = mpcg(Lop,fa(cind),xin,MM,20);
			end
		end
% up the ladders
		for clevel = nlevels-1:-1:1
            nc = nxgrid(clevel);
            mc = nygrid(clevel);
            cind = strt(clevel,1):fin(clevel,1);
            cindx = strt(clevel,2):fin(clevel,2);
            cindy = strt(clevel,3):fin(clevel,3);
            ncp = nxgrid(clevel+1);
            mcp = nygrid(clevel+1);
            pind = strt(clevel+1,1):fin(clevel+1,1);
% prolong approximate solution from previous level
			xout = FourierChebProlong(u(pind),ncp,mcp);
% add this to the solution at the current level
			u(cind) = u(cind)+xout;
% setup periodic FD approximation M to L
            MM = HelmholtzOpFD(ga(cind),K(cindx),dy(cindy),nc,mc);
% function to calculate residual at this level
			Lop = @(uc) HelmholtzOpFC(uc,K(cindx),L(cindy),ga(cind),nc,mc);
% initialize the solution
			xin = u(cind);
% do some pCG sweeps
			[u(cind),r(cind)] = mpcg(Lop,fa(cind),xin,MM,2);
        end
        rem = reshape(r(cind),[nc,mc]);
        rms = sqrt(sum(weights.*sum(rem.^2)));
        fprintf('rmse = %g \n',rms)
        if rms < tolerance
			break;
		end
	end
% calculate the final solution and residual
	uf = reshape(real(u(strt(1):fin(1))),[nc,mc]);
	rf = reshape(real(r(strt(1):fin(1))),[nc,mc]);
end




