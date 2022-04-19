% Function performs natural parameter continuation
%
% Input:
% v - initial solution at initial parameter
% u - initial parameter
% s - initial step size
% ds - step size
% dsmax
% dsmin
% DJL
% domain
% option
%
% Output:
% v - solution vector at each parameter value 
% u - parameter vector

% function [v,lambda]=naturalparametercontinuation(v,u,du,domain)

N=domain.N;
x=domain.x;
dx=domain.dx;

mu=DJL.mu;
L=DJL.L;
u=pde.u;
KAI=DJL.KAI;


U(1)=u;
V(1)=clenshaw_curtis(2*trapI(v.^2,dx{1})'/pi*KAI*L/mu);

converge=true;

for ii=1:250
    fprintf('Step %d\n',ii)
    u=u+du;
    DJL.u=u;
    v0=v;
    [pde,domain]=DJL_pde_initialise(DJL,domain);
    v=NewtonSolve(v0,pde,domain,option);

    if converge==false
        fprintf('OH NOOOO')
        break
    end

    U(ii+1)=u;
    V(ii+1)=clenshaw_curtis(2*trapI(v.^2,dx{1})'/pi*KAI*L/mu);

    % check dv<1 requirement
    dv=2*ifct(chebdiff(fct(v'),1));
    if max(dv(:))>1
        break
    end
end
