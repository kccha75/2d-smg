% Function uses Preconditioned Conjugate Gradient to solve system Au=f
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

function [v,r]=cg(v,pde,domain,option)
f=pde.f;
maxit=10000;

% Initial residual
r=pde.f-option.operator(v,pde,domain);

if rms(r(:))<option.tol
    fprintf('Initial guess already below tolerance!\n');
    return
end

% Check if preconditioner available
if ~isempty(option.preconditioner) && option.prenumit~=0
    
    pde.f=r;
    option.numit=option.prenumit;
    d=option.preconditioner(zeros(size(v)),pde,domain,option);
    
else
    
    d=r;
    
end

delta=sum(sum(r.*d));

for i=1:maxit
    
    q=option.operator(d,pde,domain);
    alpha=delta/sum(sum(d.*q));
    v=v+alpha*d; % estimate of new solution
    
    % Recalculate every 100 iteration to remove roundoff errors
    if mod(i,100)==0
        
        r=f-option.operator(v,pde,domain);
        
    else
        
        r=r-alpha*q;
        
    end
    
    % Tolerance check
    if rms(r(:))<option.tol
         
        fprintf('Conjugate Gradient Converged after %d iterations!\n',i);
        break
        
    end
    
    % update for next guess
    % Check if preconditioner available
    if ~isempty(option.preconditioner) && option.prenumit~=0
    
        pde.f=r;
        option.numit=option.prenumit;
        s=option.preconditioner(zeros(size(v)),pde,domain,option);
    
    else
    
        s=r;
    
    end
     
    deltanew=sum(sum(r.*s));
    beta=deltanew/delta;
    d=s+beta*d;
    delta=deltanew;
    
end

if i==maxit
    
    fprintf('Conjugate Gradient did not converge after %d iterations!\n',i);
    fprintf('Residual %d\n',rms(r(:)));
    
end

end