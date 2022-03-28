du=0.005;

global u
U(1)=u;
V(1)=clenshaw_curtis2(2*trapint(v(1:N(1)/2+1,:).^2,x{1}(1:N(1)/2+1))'/pi*KAI*L/mu);

for ii=1:1
    fprintf('Step %d\n',ii)
    u=u+du;
    
    v0=v;
    DJLsolve

    if converge==false
        fprintf('OH NOOOO')
        break
    end

    U(ii+1)=u;
    V(ii+1)=clenshaw_curtis2(2*trapint(v(1:N(1)/2+1,:).^2,x{1}(1:N(1)/2+1))'/pi*KAI*L/mu);

    % check dv<1 requirement
    dv=2*ifct(chebdiff(fct(v'),1));
    if max(dv(:))>1
        break
    end
end

plot(U,V)
xlabel('u');ylabel('Momentum');title('mode 1 DJL')

figure
contour(X2,Y2,Y2-v,100)
title("C=" + U(end))

figure;
plot(X2,Y2-v)
title("C=" + U(end))

% check dv<1 requirement
dv=2*ifct(chebdiff(fct(v'),1));
max(dv(:))
min(dv(:))