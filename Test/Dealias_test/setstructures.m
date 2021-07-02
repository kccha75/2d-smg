% Function to set structures

function [jacobian,domain]=setstructures(v0,pde,jacobian,domain,option)


% Find dimension of problem
if isvector(v0)
    d=1;
else
    d=ndims(v0);
end

solution(1).v=v0;

% Loop to set parameters for each grid level

for i=2:option.grids

    for j=1:d
        domain(i).L(j)=domain(1).L(j); % L does not change between grids
        domain(i).N(j)=domain(i-1).N(j)/2;
        domain(i).k{j}=[domain(i-1).k{j}(1:domain(i-1).N(j)/4);domain(i-1).k{j}(3*domain(i-1).N(j)/4+1);domain(i-1).k{j}(3*domain(i-1).N(j)/4+2:end)];
        domain(i).dx(j)=domain(1).L(j)/domain(i).N(j);
        
        % Ensure zeros(domain(option.grids).N); gives a vector
        if d==1
            domain(i).N(2)=1;
        end
        
        solution(i).v=option.restriction(solution(i-1).v);
        
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

    % define jacobian here??
    jacobian(i)=option.jacobian_dealias(solution(i).v,pde(i),domain(i));
    
end

end