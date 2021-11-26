% Function to set structures for multigrid
%
% Inputs:
% v0 - initial guess
% pde - struct of pde parameters
% domain.N 
% domain.k
% domain.dx 
% option.discretisation - discretisation flags for each direction
% option.grids - total number of grids in multigrid used
%
% Outputs:
% pde - struct of pde parameters for each grid in multigrid
% domain - struct of domain parameters for each grid in multigrid

function [pde,domain]=setstructures(v0,pde,domain,option)

% Find dimension of problem
if isvector(v0)
    d=1;
else
    d=ndims(v0);
end

% Loop to set parameters for each grid level
for i=2:option.grids
    
    % Loop dimensions
    for j=1:d
        
        domain(i).N(j)=ceil(domain(i-1).N(j)/2);
        
        switch option.discretisation(j)
            
            case 1  
            % Fourier
            domain(i).k{j}=[domain(i-1).k{j}(1:domain(i-1).N(j)/4);domain(i-1).k{j}(3*domain(i-1).N(j)/4+1);domain(i-1).k{j}(3*domain(i-1).N(j)/4+2:end)];
            domain(i).dx{j}=domain(i-1).dx{j}/2;
                
            case 2
            % Cheb
            domain(i).k{j}=domain(i-1).k{j}(1:domain(i).N(j));
            domain(i).dx{j}=domain(i-1).dx{j}(2:2:end)+domain(i-1).dx{j}(1:2:end-1);
            
        end
        
        % Ensure zeros(domain(option.grids).N); gives a vector
        if d==1
            domain(i).N(2)=1;
        end

    end
    
    % Loop though all pde coefficients and restrict to coarse grid
    fn=fieldnames(pde);
    for ii=1:numel(fn)
        
        % Constant coefficient check
        if length(pde(i-1).(fn{ii}))==1
            pde(i).(fn{ii})=pde(i-1).(fn{ii});
        else
            pde(i).(fn{ii})=option.restriction(pde(i-1).(fn{ii}));
        end

    end
    
end

end