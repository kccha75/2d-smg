function cttest
    N = 32;
    k = 0:N;
    x = cos(pi*k/N);
    
    theta = tanh(x);
    u = theta+3*theta.^2;
    ud = (1+6*theta).*(1-theta.^2);
    z = fct(u);
    figure(1); semilogy(k, abs(z));
    uu = bct(z);
    figure(2); plot(x,u-uu);
    uud = bct(ctd(z,k));
    figure(3); plot(x,ud,'r',x,uud,'g');
end

function [fc] = fct(ff)
    N = length(ff)-1;
    ff(1:N:N+1) = ff(1:N:N+1)/sqrt(2);
    fc = dct(ff,'Type',1);
    fc(1:N:N+1) = fc(1:N:N+1)/sqrt(2);
end

function [fc] = bct(ff)
    N = length(ff)-1;
    ff(1:N:N+1) = ff(1:N:N+1)*sqrt(2);
    fc = dct(ff,'Type',1);
    fc(1:N:N+1) = fc(1:N:N+1)*sqrt(2);
end

function [fc] = ctd(ff,k)
    N = length(ff)-1;
    
    % why ?
    fc(N+1) = .0; % okkkk
    % implies this line?
    fc(N) = 2*k(N+1)*ff(N+1);
    
    
    % ok i get this part ...
    for j = N-1:-1:2
        fc(j) = fc(j+2)+2*k(j+1)*ff(j+1);
    end
    fc(1) = .5*fc(3)+ff(2);
end
    