% Performs Minimum Residual Richardson smoother on Au=f
%
% Inputs:
% v - best estimate
% pde - structure consisting of pde coefficients
% domain.k - wave number
% option.operator - operator used to find Lu
% option.numit - number of iterations
% option.preconditioner - preconditioner used
% option.prenumit - number of preconditioner iterations
%
% Ouputs:
% v - best estimate after relaxation
% r - residual after relaxation
%
% See Boyd / Canuto for reference to algorithm

function [v,r] = MRR(v,pde,domain,option)

numit=option.numit;

% Find residual only if numit=0
if numit==0
    r=pde.f-option.operator(v,pde,domain);
    return;
end

r=pde.f-option.operator(v,pde,domain);

if r==0
    return
end

% Check if preconditioner available
if ~isempty(option.preconditioner) && option.prenumit~=0
    
    pde.f=r;
    option.numit=option.prenumit;
    z=option.preconditioner(zeros(size(v)),pde,domain,option);
    
else
    
    z=r;
    
end

for i=1:numit
    
    Az=option.operator(z,pde,domain);
    tau=sum(sum(r.*Az))/sum(sum(Az.*Az));
    v=v+tau*z;
    r=r-tau*Az;
    
    % Preconditioner step only if not at last iteration
    if i~=numit
        
        if ~isempty(option.preconditioner) && option.prenumit~=0
            
            pde.f=r;
            option.numit=option.prenumit;
            z=option.preconditioner(zeros(size(v)),pde,domain,option);
            
        else
            
            z=r;
            
        end
    end
    
end

end