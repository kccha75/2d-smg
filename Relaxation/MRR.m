% Performs Minimum Residual Richardson smoother on Lv=f
%
% Inputs:
% v - best estimate
% pde.a
% pde.b
% pde.c
% pde.f
% scheme.kx - wave number in x
% scheme.ky - wave number in y
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

function [v,r] = MRR(v,pde,scheme,option)

numit=option.numit;

% Find residual only if numit=0
if numit==0
    r=pde.f-option.operator(v,pde,scheme);
    return;
end

r=pde.f-option.operator(v,pde,scheme);

if r==0
    return
end

% Check if preconditioner available
if ~isempty(option.preconditioner) && option.prenumit~=0
    
    pde.f=r;
    option.numit=option.prenumit;
    z=option.preconditioner(zeros(size(v)),pde,scheme,option);
    
else
    
    z=r;
    
end

for i=1:numit
    
    Az=option.operator(z,pde,scheme);
    tau=sum(sum(r.*Az))/sum(sum(Az.*Az));
    v=v+tau*z;
    r=r-tau*Az;
    
    % Preconditioner step only if not at last iteration
    if i~=numit
        
        if ~isempty(option.preconditioner) && option.prenumit~=0
            
            pde.f=r;
            option.numit=option.prenumit;
            z=option.preconditioner(zeros(size(v)),pde,scheme,option);
            
        else
            
            z=r;
            
        end
    end
    
end

end