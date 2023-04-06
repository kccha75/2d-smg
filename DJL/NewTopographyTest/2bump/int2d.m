function intans=int2d(v,domain,DJL)

intans=v.*1./domain.jac;
intans=trapI(intans.^2,DJL.Lx/domain.N(1));
intans=permute(intans,[2,1,3]);
intans=clenshaw_curtis(intans)/2*DJL.Ly/domain.H;

end