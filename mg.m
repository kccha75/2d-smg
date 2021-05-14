% Function uses multigrid to solve given pde
function [v,r]=mg(v0,pde,domain,scheme,option)

% -------------------------------------------------------------------------
% Set up structures
% -------------------------------------------------------------------------

solution(1).v=v0;

% Loop to set parameters for each grid level
for i=2:option.grids

    domain(i).Nx=domain(i-1).Nx/2;
    domain(i).Ny=domain(i-1).Ny/2;
    
    % unfiltered
    scheme(i).kx=[scheme(i-1).kx(1:domain(i-1).Nx/4);scheme(i-1).kx(3*domain(i-1).Nx/4+1);scheme(i-1).kx(3*domain(i-1).Nx/4+2:end)];
    scheme(i).ky=[scheme(i-1).ky(1:domain(i-1).Ny/4);scheme(i-1).ky(3*domain(i-1).Ny/4+1);scheme(i-1).ky(3*domain(i-1).Ny/4+2:end)];
    
    % dx dy
    scheme(i).dx=domain(1).Lx/domain(i).Nx;
    scheme(i).dy=domain(1).Ly/domain(i).Ny;
    
    % Check if coefficients constant
    if length(pde(1).a)==1
        pde(i).a=pde(i-1).a;
    else
        pde(i).a=option.restriction(pde(i-1).a);
    end
    
    if length(pde(1).b)==1
        pde(i).b=pde(i-1).b;
    else
        pde(i).b=option.restriction(pde(i-1).b);
    end
    
    if length(pde(1).c)==1
        pde(i).c=pde(i-1).c;
    else
        pde(i).c=option.restriction(pde(i-1).c);
    end   

    % step down RHS (for FMG only)
    if option.solver == "FMG"
        
        pde(i).f=option.restriction(pde(i-1).f);
        
    end
    
end

% -------------------------------------------------------------------------
% Solution method here
% -------------------------------------------------------------------------

% Perform initial relaxations on finest grid?
option.initialrelaxations=true;

switch option.solver
    
    case 'V-cycle'
        
        % Set finest grid and V-cycle grids
        option.initialgrid=1;
        
        % Solve using V-cycle
        [solution.v,solution.r]=vcycle(v0,pde,domain,scheme,option);
    
    case 'FMG'
             
        % Coarse grid initial guess (zeros)
        solution(option.grids).v=zeros(domain(option.grids).Nx,domain(option.grids).Ny);
        
        % Coarse grid solver
        solution(option.grids).v=option.coarsegridsolver(solution(option.grids).v,pde(option.grids),scheme(option.grids),option);
        
        % FMG options at coarsest grid
        option.initialgrid=option.grids;
        option.grids=1;
        
        % Save num_vcycle parameter
        num_vcycles=option.num_vcycles;
        option.num_vcycles=1;
        
        j=option.finestgrid-option.coarsestgrid+1;
        
        for i=j:-1:2
            
            % Increment for next FMG cycle
            option.initialgrid=option.initialgrid-1;
            option.grids=option.grids+1;
            
            % Step up solution
            solution(i-1).v=option.prolongation(solution(i).v);
            
            if i~=2
                
                % Do one V-cycle
                [solution(i-1).v,solution(i-1).r]=vcycle(solution(i-1).v,pde,domain,scheme,option);
                
            else
                
            	% Do 1 V-cycles + num_vcycles V-cycles after FMG
                option.num_vcycles=num_vcycles+1;
                [solution(option.initialgrid).v,solution(option.initialgrid).r]=vcycle(solution(option.initialgrid).v,pde,domain,scheme,option);
            
            end
            
        end
        

    otherwise
        
        disp('Invalid Solver!\n')
        return
        
end

% Print solution
v=solution(1).v;
r=solution(1).r;

end

% -------------------------------------------------------------------------
% Subroutine V-cycle
% -------------------------------------------------------------------------

function [v,r]=vcycle(v,pde,domain,scheme,option)

solution(option.initialgrid).v=v;

% Down cycle iterations
option.numit=option.Nd;

% Initial iterations
[solution(option.initialgrid).v,solution(option.initialgrid).r]=option.relaxation(solution(option.initialgrid).v,pde(option.initialgrid),scheme(option.initialgrid),option);
    

% loop through vcycles
for p=1:option.num_vcycles

    % Stepping down
    for i=option.initialgrid:option.initialgrid+option.grids-2
        
        % Correction vs FAS scheme?
        switch option.mgscheme
            
            case 'Correction'
                
                % Step down
                pde(i+1).f=option.restriction(solution(i).r);
                % clear v from previous loop
                solution(i+1).v=zeros(domain(i+1).Nx,domain(i+1).Ny);
                
            case 'FAS'

                % Step down solution
                solution(i+1).v=option.restriction(solution(i).v);
                % Calculate RHS
                pde(i+1).f=option.operator(solution(i+1).v,pde(i+1),scheme(i+1))+option.restriction(solution(i).r);
                
            otherwise
                
                disp('Invalid mgscheme selected!\n')
                return
                        
        end
        
        % If at coarsest level, solve exactly
        if i==option.initialgrid+option.grids-2
            
            solution(i+1).v=option.coarsegridsolver(solution(i+1).v,pde(i+1),scheme(i+1),option);
            
        else
            
            % Iterate
            [solution(i+1).v,solution(i+1).r]=option.relaxation(solution(i+1).v,pde(i+1),scheme(i+1),option);

        end
        
    end

    % Stepping up
    
    % Up cycle iterations
    option.numit=option.Nu;
    
    for i=option.initialgrid+option.grids-1:-1:option.initialgrid+1
        
        % Correction vs FAS scheme?
        switch option.mgscheme
            
            case 'Correction'
                
                % Step up and update guess
                solution(i-1).v=solution(i-1).v+option.prolongation(solution(i).v);
                               
            case 'FAS'
                
                % Step up and update guess
                solution(i-1).v=solution(i-1).v+option.prolongation(solution(i).v-option.restriction(solution(i-1).v));
                
        end
        
         % Iterate
         [solution(i-1).v,solution(i-1).r]=option.relaxation(solution(i-1).v,pde(i-1),scheme(i-1),option);
                        
    end

    disp(rms(rms(solution(option.initialgrid).r))) % optional display rms residual
    
    if rms(rms(solution(option.initialgrid).r))<option.tol
        
        fprintf('SMG converged after %d V-cycles!\n',p);
        break
        
    end
    
end

v=solution(option.initialgrid).v;
r=solution(option.initialgrid).r;

end