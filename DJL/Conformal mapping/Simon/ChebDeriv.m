function [fc] = ChebDeriv(ff,k,nd)
	[N,~] = size(ff);
    fc = ff;
    for id = 1:nd
        fc(N,:) = .0;
        fc(N-1,:) = 2*k(N)*ff(N,:);
        for j = N-2:-1:2
            fc(j,:) = fc(j+2,:)+2*k(j+1)*ff(j+1,:);
        end
        fc(1,:) = .5*fc(3,:)+k(2)*ff(2,:);
        ff = fc;
    end
end
