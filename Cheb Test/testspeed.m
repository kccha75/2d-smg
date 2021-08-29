    N = 2^9;
    k = 0:N-1;
    x = cos(pi*k/(N-1))';
    
    theta = tanh(x);
    u = theta+3*theta.^2;
    ud = (1+6*theta).*(1-theta.^2);
    
    udd=6*sech(x).^4-2*tanh(x).*(6*tanh(x)+1).*sech(x).^2;
    
   z = fct(u);
    uud = bct2(ctd2(z,k));
    
    
    tic
    for i=1:10000
        du=fchd2(u);
    end
    
   toc
   
   tic
   for i=1:10000
       du2=ifct(chebdiff(fct(u),2));
   end
   toc
    
    
function [fc] = fct2(ff)
    N = length(ff)-1;
    ff(1:N:N+1) = ff(1:N:N+1)/sqrt(2);
    fc = dct(ff,'Type',1);
    fc(1:N:N+1) = fc(1:N:N+1)/sqrt(2);
end

function [fc] = bct2(ff)
    N = length(ff)-1;
    ff(1:N:N+1) = ff(1:N:N+1)*sqrt(2);
    fc = dct(ff,'Type',1);
    fc(1:N:N+1) = fc(1:N:N+1)*sqrt(2);
end

function [fc] = ctd2(ff,k)
    N = length(ff)-1;
    fc(N+1) = .0;
    fc(N) = 2*k(N+1)*ff(N+1);
    for j = N-1:-1:2
        fc(j) = fc(j+2)+2*k(j+1)*ff(j+1);
    end
    fc(1) = .5*fc(3)+ff(2);
end