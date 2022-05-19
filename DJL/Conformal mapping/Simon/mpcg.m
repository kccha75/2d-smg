function [xout,rem] = mpcg(Lm,fam,xin,Mm,nrelax)
% does a few iterations of standard preconditioned conjugate gradient
% need option to iterate to specific tolerance or just make nrelax large
% and bail when tolerance satisfied
	tolerance = 1.e-12;
    rem = fam-Lm(xin);
    zm = Mm(rem);
    pm = zm;
    gmp = sum(rem.*zm);
    for n = 1:nrelax
        Lp = Lm(pm);
        omega = sum(rem.*zm)/sum(pm.*Lp);
        xin = xin+omega*pm;
        rem = rem-omega*Lp;
        if norm(rem) < tolerance
			break;
        end
        zm = Mm(rem);
        gm = sum(rem.*zm);
        pm = zm+(gm/gmp)*pm;
        gmp = gm;
    end
    xout = xin;
end