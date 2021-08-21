% Function uses Bi Conjugate Gradient to solve system Au=f
%
% Iterates until tolerance condition OR maximum iterations reached
%
% Inputs:
% v - best estimate
% pde - structure consisting of pde coefficients
% scheme.k - wave number in x
% option.operator - operator used to find Lu
% option.numit - number of iterations
% option.preconditioner - preconditioner used
% option.prenumit - number of preconditioner iterations
%
% Ouputs:
% v - solution after tol / maxit reached
% r - residual for solution
%
% Optional display messages can be commented out or left in

function [v,r]=bicg(v,pde,domain,option)
f=pde.f;
maxit=1000;

% Initial residual
r=pde.f-option.operator(v,pde,domain);
r_hat=r;

if rms(r)<option.tol
    return
end

% Check if preconditioner available
if ~isempty(option.preconditioner) && option.prenumit~=0
    
    pde.f=r;
    option.numit=option.prenumit;
    d=option.preconditioner(zeros(size(v)),pde,domain,option);
    d_hat=d;
    
else
    
    d=r;
    d_hat=d;
    
end

delta=sum(sum(r_hat.*d));

for i=1:maxit
    
    q=option.operator(d,pde,domain);
    
    % tranpose here
    pde_trans.a=pde.b;
    pde_trans.b=pde.a;
    pde_trans.c=pde.c;
    domain_trans=domain;
    domain_trans.k(:,1)=domain.k(:,2);
    domain_trans.k(:,2)=domain.k(:,1);
    
    q_hat=option.operator(d_hat,pde_trans,domain_trans);
    alpha=delta/sum(sum(d_hat.*q));
    v=v+alpha*d; % estimate of new solution
    
    % Recalculate every 100 iteration to remove roundoff errors
    if mod(i,100)==0
        
        r=f-option.operator(v,pde,domain);
        r_hat=r;
        
    else
        
        r=r-alpha*q;
        r_hat=r_hat-alpha*q_hat;
        
    end
    
    % Tolerance check
    if rms(r)<option.tol
         
        fprintf('Conjugate Gradient Converged after %d iterations!\n',i);
        break
        
    end
    
    % update for next guess
    % Check if preconditioner available
    if ~isempty(option.preconditioner) && option.prenumit~=0
    
        pde.f=r;
        option.numit=option.prenumit;
        s=option.preconditioner(zeros(size(v)),pde,domain,option);
        pde.f=r_hat;
        s_hat=option.preconditioner(zeros(size(v)),pde,domain,option);
        
    else
    
        s=r;
        s_hat=r_hat;
    
    end
     
    deltanew=sum(sum(r_hat.*s));
    beta=deltanew/delta;
    
    d=s+beta*d;
    d_hat=s_hat+beta*d_hat;
    
    delta=deltanew;
    
end

if i==maxit
    
    fprintf('Bicg did not converge after %d iterations!\n',i);
    fprintf('Residual %d\n',rms(rms(r)));
    
end

end