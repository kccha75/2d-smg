% Function uses multigrid to solve given pde
%
% Inputs:
% v0 - initial guess
% pde - structure consisting of pde coefficients
%
% domain.N - number of points
% domain.k - spectral wave numbers
% domain.dx - step size
%
% option.num_vcycles - number of vcycles performed
% option.tol - tolerance set (solution / iteration)
% option.Nd - number of iterations on the down cycle
% option.Nu - number of iterations on the up cycle
% option.solver - multigrid solver used 'V-cycle' or 'FMG'
% option.mgscheme - multigrid scheme used 'Correction' or 'FAS'
% option.operator - operator function to calculate L*u
% option.coarsegridsolver - coarse grid solver function
% option.restriction - restriction operator for pde coefficients
% option.restriction_residual - restriction operator for residual
% option.restriction - restriction function
% option.prolongation - prolongation function
% option.preconditioner - preconditioner function
% option.prenumit - number of precondition iterations
% option.finestgrid - finest grid used given by 2^finestgrid (or +1)
% option.coarsestgrid - coarsest grid used given by 2^coarsestgrid (or +1)
% option.grids - numbers of grids used in mg (counting from coarsest grid)
%
% Outputs:
% v - best estimate after multigrid
% r - residual of v

% Notes: assumes constant dx dy (for now)
function [v,r]=mg(v0,pde,domain,option)

% Solution structure here
solution(1).v=v0;

% Set multigrid structures
[pde,domain]=setstructures(v0,pde,domain,option);

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
        [solution.v,solution.r]=vcycle(v0,pde,domain,option);
    
    case 'FMG'
             
        % Coarse grid initial guess (zeros)
        solution(option.grids).v=zeros(domain(option.grids).N);
        
        % Coarse grid solver
        solution(option.grids).v=option.coarsegridsolver(solution(option.grids).v,pde(option.grids),domain(option.grids),option);
        
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
                [solution(i-1).v,solution(i-1).r]=vcycle(solution(i-1).v,pde,domain,option);
                
            else
                
            	% Do 1 V-cycles + num_vcycles V-cycles after FMG
                option.num_vcycles=num_vcycles+1;
                [solution(option.initialgrid).v,solution(option.initialgrid).r]=vcycle(solution(option.initialgrid).v,pde,domain,option);
            
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

function [v,r]=vcycle(v,pde,domain,option)

solution(option.initialgrid).v=v;

% Down cycle iterations
option.numit=option.Nd;

% Initial iterations
[solution(option.initialgrid).v,solution(option.initialgrid).r]=option.relaxation(solution(option.initialgrid).v,pde(option.initialgrid),domain(option.initialgrid),option);
    

% loop through vcycles
for p=1:option.num_vcycles

    % Stepping down
    for i=option.initialgrid:option.initialgrid+option.grids-2
        
        % Correction vs FAS scheme?
        switch option.mgscheme
            
            case 'Correction'
                
                % Step down
                pde(i+1).f=option.restriction_residual(solution(i).r);
                % clear v from previous loop
                solution(i+1).v=zeros(domain(i+1).N);
                
            case 'FAS'

                % Step down solution
                solution(i+1).v=option.restriction_residual(solution(i).v);
                % Calculate RHS
                pde(i+1).f=option.operator(solution(i+1).v,pde(i+1),domain(i+1))+option.restriction_residual(solution(i).r);
                
            otherwise
                
                disp('Invalid mgscheme selected!\n')
                return
                        
        end
        
        % If at coarsest level, solve exactly
        if i==option.initialgrid+option.grids-2
            
            solution(i+1).v=option.coarsegridsolver(solution(i+1).v,pde(i+1),domain(i+1),option);
            
        else
            
            % Iterate
            [solution(i+1).v,solution(i+1).r]=option.relaxation(solution(i+1).v,pde(i+1),domain(i+1),option);

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
                solution(i-1).v=solution(i-1).v+option.prolongation(solution(i).v-option.restriction_residual(solution(i-1).v));
                
        end

         % Iterate
         [solution(i-1).v,solution(i-1).r]=option.relaxation(solution(i-1).v,pde(i-1),domain(i-1),option);
                        
    end

    disp(rms(solution(option.initialgrid).r(:))) % optional display rms residual
    
    if rms(solution(option.initialgrid).r(:))<option.tol
        
        fprintf('SMG converged after %d V-cycles!\n',p);
        break
        
    end
    
end

v=solution(option.initialgrid).v;
r=solution(option.initialgrid).r;

end