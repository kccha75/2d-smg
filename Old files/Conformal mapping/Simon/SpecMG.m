%
% Last updated: 20210902
%
% 2-dim SMG function to solve equations of the form
% L{u} = -nabla^2 u + g*u = f , (1)
% on torus [0,Lx) x [0,Ly), periodic strip or rectangle
% Fourier on periodic directions, Chebyshev on finite intervals
% boundary conditions can be periodic, homogeneous Dirichlet or Neumann
%
% theoretically, this should be able to do 1-dim problems on 
% periodic and finite intervals through appropriate configuration
%
% Dirichlet or Neumann bounday conditions should only appear in the 
% pcond and specop functions, hence should be no difference
% could maybe even handle Robin boundary conditions
%
% uses PCG for relaxations and solution on coarsest grid
% pcond returns an approximate solution to the problem
% specop returns the LHS of (1), including boundary conditions
%
% need to test whether PCG is appropriate, and have alternative solvers
%
% framework could be used for 2-dim FDMG, where preconditioners are ADI
%
% todo:
% * framework for rectangle, and 1-dim problems
%   * for rectangle need to modify ChebRestrictZero so that works for 
%     intervals other than y in [-1,1]
%   * 1-dim problems, make nygrid ones of length nlevels and 
%     transfer operators just are identities
% * standardize arguments for operators, and how to pass parameters,
%   probably global is the easiest way
% * test whether passing of operator functions work
% * write 2-dim operations in terms of 1-dim operations, but test efficiency
% * add to git
% * test suite
% * combine with boundary integral method for conformal mapping to handle
%   topography
% * handling undetermined systems, probably only required in pcond and specop
% * initial approximation, either output of pcond or a previous solution
% * parameter to handle structure of cycles and stopping criteria
% * self-adjoint forms of nabla^2, e.g., spherical
%
function [uf,rf] = SpecMG(topology,x,y,g,f,pcond,specop)
%
	ci = complex(0,1);
%
% check that inputs are correct
%
	[N,M] = size(g); % number of points in the finest domain
	if length(x) ~= N || length(y) ~= M || sum(size(f) ~= size(g))
		fprintf('Array sizes incompatible. \n');
		whos x y f g
		return;
	end	
%
% parameters for multigrid
%
	nlevels = 4; % number of levels
    nmgit = 1; % number of iterations
    tolerance = 1.e-12; % desired accuracy
    nrelax = 2; % number of iterations for relaxation
    nsolve = 20; % max number of iterations for solving on bottom grid
%
% setup functions dependent on topology
%
% functions for restricting wavenumber
%
    function kc = FourierWvnoRestrict(kk,nin)
		kc = [kk(1:nin/4-1,:); kk(3*nin/4:nin,:)];
	end
    function kc = ChebWvnoRestrict(kk,nin)
		kc = [kk(1:(nin-1)/2+1,:)];
	end
%
	switch topology
		case 0 
% torus
			grestrict = @FourierFourierRestrict;
			frestrict = @FourierFourierRestrict;
			rprolong = @FourierFourierProlong;
			krestrict = @FourierWvnoRestrict;
			lrestrict = @FourierWvnoRestrict;
%
			Lx = N*(x(2)-x(1));
			Ly = M*(y(2)-y(1));
			K1 = 2*pi/Lx*ci*[0:N/2-1 -N/2:-1]';
			L1 = 2*pi/Ly*ci*[0:M/2-1 -M/2:-1]';;
			eta = Lx*(0:N-1)/N;
			zeta = Ly*(0:M-1)/M;
			weights = Lx*Ly/(N*M);
			nxgrid = (N./2.^(0:nlevels-1))'; % array with number of points in x direction
			nygrid = (M./2.^(0:nlevels-1))'; % array with number of points in y direction
		case 1
% periodic strip
			grestrict = @FourierChebRestrict;
			frestrict = @FourierChebRestrictZero;
			rprolong = @FourierChebProlong;
			krestrict = @FourierWvnoRestrict;
			lrestrict = @ChebWvnoRestrict;
%
			Ly = (y(1)-y(M))/2.;
			Lx = N*(x(2)-x(1));
			K1 = 2*pi/Lx*ci*[0:N/2-1 -N/2:-1]';
			L1 = [0:M-1];
			f = Ly^2*f;
			g = Ly^2*g;
			zeta = cos(pi*(0:M-1)/(M-1));
			eta = Lx/Ly*(0:N-1)/N;
			weights = sqrt(1-zeta.^2);
			rem = ones(N,M);
			wnorm = sum(weights.*sum(rem));
			weights = weights/wnorm;
			nxgrid = (N./2.^(0:nlevels-1))'; % array with number of points in x direction
			nygrid = ((M-1)./2.^(0:nlevels-1)+1)'; % array with number of points in y direction
		case 2
% rectangle
% need to adapt zero restriction so that it takes into account intervals other
% than [-1,1]	
			grestrict = @ChebChebRestrict;
			frestrict = @ChebChebRestrictZero;
            frestrict = @ChebChebRestrict;
			rprolong = @ChebChebProlong;
			krestrict = @ChebWnoRestrict;
			lrestrict = @ChebWnoRestrict;
		case 3
% periodic interval
			nxgrid = (N./2.^(0:nlevels-1))'; 
			nygrid = ones(nlevels,1); 

	end
%
%%%%%%%%%%%%%%%% setup %%%%%%%%%%%%%%%%%%%
%
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
	K = zeros(dxlength,1); % array of x wavenumbers
	L = zeros(dylength,1); % y wavenumbers
	xa = zeros(dxlength,1); % array of x wavenumbers
	ya = zeros(dylength,1); % y wavenumbers
%
    clevel = 1;
    nc = nxgrid(clevel);
    mc = nygrid(clevel);
    cind = strt(clevel,1):fin(clevel,1);
    cindx = strt(clevel,2):fin(clevel,2);
    cindy = strt(clevel,3):fin(clevel,3);
	fa(cind) = reshape(f,[nc*mc,1]); % initialize RHS on top domain
	ga(cind) = reshape(g,[nc*mc,1]); % initialize LHS function on top domain
	K(cindx) = K1;
	L(cindy) = L1;
	xa(cindx) = eta;
	ya(cindy) = zeta;
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
		ga(cind) = grestrict(ga(pind),ncp,mcp);
		K(cindx) = krestrict(K(pindx),ncp);
		L(cindy) = lrestrict(L(pindy),mcp);
		xa(cindx) = eta(1:steps:end);
		ya(cindy) = zeta(1:steps:end);
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
    MM = @(fin) pcond(ga(cind),fin,xa(cindx),ya(cindy),K(cindx),L(cindy),nc,mc);
% function to calculate residual at top domain
    Lop = @(uc) specop(uc,ga(cind),xa(cindx),ya(cindy),K(cindx),L(cindy),nc,mc);
% use FD to obtain an initial approximation to the solution
	xin = MM(fa(cind));
	r(cind) = Lop(xin)-fa(cind);
	rem = reshape(r(cind),[nc,mc]);
	rmsp = sqrt(sum(weights.*sum(rem.^2)));
	uf = xin;
	rf = rem;
	fprintf('n = %g, rmse = %g \n',0,rmsp)
	if rmsp < tolerance
		return;
	end
% use pCG as relaxation to improve solution and calculate residual
	[u(cind),r(cind)] = mpcg(Lop,fa(cind),xin,MM,nrelax);
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
			fa(cind) = frestrict(r(pind),ncp,mcp);
% setup periodic FD approximation M to L
			MM = @(fin) pcond(ga(cind),fin,xa(cindx),ya(cindy),K(cindx),L(cindy),nc,mc);
% function to calculate residual at top domain
			Lop = @(uc) specop(uc,ga(cind),xa(cindx),ya(cindy),K(cindx),L(cindy),nc,mc);
% use FD to obtain an initial approximation to the solution
			xin = MM(fa(cind));
% if we are at an intermediate level do pCG relaxation
% if we are at the bottom level use pCG to solve the problem
% should have parameter in pCG to iterate to solution
			if (clevel < nlevels)
				[u(cind),r(cind)] = mpcg(Lop,fa(cind),xin,MM,nrelax);
			else
				[u(cind),r(cind)] = mpcg(Lop,fa(cind),xin,MM,nsolve);
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
			xout = rprolong(u(pind),ncp,mcp);
% add this to the solution at the current level
			u(cind) = u(cind)+xout;
% setup periodic FD approximation M to L
			MM = @(fin) pcond(ga(cind),fin,xa(cindx),ya(cindy),K(cindx),L(cindy),nc,mc);
% function to calculate residual at top domain
			Lop = @(uc) specop(uc,ga(cind),xa(cindx),ya(cindy),K(cindx),L(cindy),nc,mc);
% initialize the solution
			xin = u(cind);
% do some pCG sweeps
			[u(cind),r(cind)] = mpcg(Lop,fa(cind),xin,MM,nrelax);
        end
        rem = reshape(r(cind),[nc,mc]);
        rms = sqrt(sum(weights.*sum(rem.^2)));
		fprintf('n = %g, rmse = %g \n',iv,rms)
		uf = reshape(real(u(strt(1):fin(1))),[nc,mc]);
		rf = reshape(real(r(strt(1):fin(1))),[nc,mc]);
        if rms < tolerance
			break;
		end
		if rms > rmsp
			return;
		end
		rmsp = rms;
	end
end




