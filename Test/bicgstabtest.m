% Function uses Preconditioned BiConjugate Gradient Stability to solve system Au=f
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

function [v,r]=bicgstabtest(v,pde,domain,option)
f=pde.f;
maxit=2000;

% Initial residual
r=pde.f-option.operator(v,pde,domain);
p=r;
d=r;

if rms(r)<option.tol
    return
end

% % Check if preconditioner available
% if ~isempty(option.preconditioner) && option.prenumit~=0
%     
%     pde.f=r;
%     option.numit=option.prenumit;
%     p=option.preconditioner(zeros(size(v)),pde,domain,option);
%     
% else
%     
%     p=r;
%     
% end

delta=dot(r(:),d(:));

for i=1:maxit
    
    q=option.operator(p,pde,domain);
    alpha=delta/dot(d(:),q(:));
    s=r-alpha*q;
    
    z=option.operator(s,pde,domain);% A*s
    
    w=dot(z(:),s(:))/dot(z(:),z(:));
    v=v+alpha*p+w*s; % estimate of new solution
    
    % Recalculate every 100 iteration to remove roundoff errors
    if mod(i,100)==0
        
        r=f-option.operator(v,pde,domain);
        
    else
        
        r=s-w*z;
        
    end
    
    % Tolerance check
    if rms(r)<option.tol
         
        fprintf('Bi Conjugate Gradient STAB Converged after %d iterations!\n',i);
        break
        
    end
    
%     % update for next guess
%     % Check if preconditioner available
%     if ~isempty(option.preconditioner) && option.prenumit~=0
%     
%         pde.f=r;
%         option.numit=option.prenumit;
%         s=option.preconditioner(zeros(size(v)),pde,domain,option);
%     
%     else
%     
%         s=r;
%     
%     end
     
    deltanew=dot(r(:),d(:));
    beta=(alpha/w)*deltanew/delta;
    p=r+beta*(p-w*q);
    delta=deltanew;
    
end

if i==maxit
    
    fprintf('Bi Conjugate Gradient STAB did not converge after %d iterations!\n',i);
    fprintf('Residual %d\n',rms(rms(r)));
    
end

end