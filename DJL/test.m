% x=x{1};
f=exp(x);

N=2^6+1;
k=(0:N-1)';

Nf=2^7+1;
Nff=2^8+1;
Nfff=2^9+1;
Nffff=2^10+1;
Nfffff=2^11+1;

Nc=2^5+1;
Ncc=2^4+1;
Nccc=2^3+1;
Ncccc=2^2+1;
Nccccc=2^1+1;

kf=(0:Nf-1)';
kff=(0:Nff-1)';
kfff=(0:Nfff-1)';
kffff=(0:Nffff-1)';
kfffff=(0:Nfffff-1)';

kc=(0:Nc-1)';
kcc=(0:Ncc-1)';
kccc=(0:Nccc-1)';
kcccc=(0:Ncccc-1)';
kccccc=(0:Nccccc-1)';

xf=cos(pi*kf/(Nf-1));
xff=cos(pi*kff/(Nff-1));
xfff=cos(pi*kfff/(Nfff-1));
xffff=cos(pi*kffff/(Nffff-1));
xfffff=cos(pi*kfffff/(Nfffff-1));

xc=cos(pi*kc/(Nc-1));
xcc=cos(pi*kcc/(Ncc-1));
xccc=cos(pi*kccc/(Nccc-1));
xcccc=cos(pi*kcccc/(Ncccc-1));
xccccc=cos(pi*kccccc/(Nccccc-1));


fc=cheb_restrict(f);
fcc=cheb_restrict(fc);
fccc=cheb_restrict(fcc);
fcccc=cheb_restrict(fccc);
fccccc=cheb_restrict(fcccc);

ff=cheb_prolong(f);
fff=cheb_prolong(ff);
ffff=cheb_prolong(fff);
fffff=cheb_prolong(ffff);
ffffff=cheb_prolong(fffff);

df=ifct(chebdiff(fct(f),1));

dff=ifct(chebdiff(fct(ff),1));
dfff=ifct(chebdiff(fct(fff),1));
dffff=ifct(chebdiff(fct(ffff),1));
dfffff=ifct(chebdiff(fct(fffff),1));
dffffff=ifct(chebdiff(fct(ffffff),1));

dfc=ifct(chebdiff(fct(fc),1));
dfcc=ifct(chebdiff(fct(fcc),1));
dfccc=ifct(chebdiff(fct(fccc),1));
dfcccc=ifct(chebdiff(fct(fcccc),1));
dfccccc=ifct(chebdiff(fct(fccccc),1));

f(1)-exp(1)

disp('coarse boundary error')
fc(1)-exp(1)
fcc(1)-exp(1)
fccc(1)-exp(1)
fcccc(1)-exp(1)
fccccc(1)-exp(1)

disp('fine boundary error')
ff(1)-exp(1)
fff(1)-exp(1)
ffff(1)-exp(1)
fffff(1)-exp(1)
ffffff(1)-exp(1)

disp('coarse boundary derivative error')
dfc(1)-exp(1)
dfcc(1)-exp(1)
dfccc(1)-exp(1)
dfcccc(1)-exp(1)
dfccccc(1)-exp(1)

disp('fine boundary derivative error')
dff(1)-exp(1)
dfff(1)-exp(1)
dffff(1)-exp(1)
dfffff(1)-exp(1)
dffffff(1)-exp(1)