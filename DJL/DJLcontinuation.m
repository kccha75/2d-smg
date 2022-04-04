N=domain.N;
x=domain.x;

du=0.0001;

u=pde.u;

U(1)=u;
V(1)=clenshaw_curtis2(2*trapint(v(1:N(1)/2+1,:).^2,x{1}(1:N(1)/2+1))'/pi*KAI*L/mu);

converge=true;

for ii=1:200
    fprintf('Step %d\n',ii)
    u=u+du;
    DJL.u=u;
    v0=v;
    [pde,domain]=DJL_pde_initialise(DJL,domain);
    v=DJLsolve(v0,pde,domain,option);

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

