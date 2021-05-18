% Function calculates num_vcycles loops of vcycles using choice of
% relaxation

% Vcycle will attempt to solve exactly at a specified coarse grid level
% using matlab backslash
% This avoids cycling down coarser grids

% Inputs:
% relaxation - function input for choice of relaxation
% initialgrid - finest grid index (1 = finest grid, 2 = next finest etc)
% grids - number of grids to cycle through
% cellN - cell of grid points
% cellk - cell of wave number
% tol - error tolerance required
% Nd - number of iterations per grid down cycle
% Nu - number of iterations per grid up cycle
% cella - cell of function a(x) in -u_xx+au_x+bu=f
% cellb - cell of function b(x) in -u_xx+au_x+bu=f
% cellf - cell of RHS of Au=f
% cellv - cell of best estimate of solution (or initial estimate)

% Ouputs:
% v - best guess after vcycle iterations
% r - residual of best guess

% For use in FMG ideally

% -------------------------------------------------------------------------

function [v,r]=vcycle_pre_tol(pre_relaxation,initialgrid,grids,cellN,cellk, ...
    tol,Nd,Nu,cella,cellb,cellf,cellH,cellv)

% Maximum vcycles before breaking
num_vcycles_max=10;

% Initial iterations
cellv{initialgrid}=pre_relaxation(cellv{initialgrid},cellk{initialgrid}, ...
    cella{initialgrid},cellb{initialgrid},cellf{initialgrid},cellH{initialgrid},Nd);

% Find initial residual
r=findR(cellf{initialgrid},cellv{initialgrid},cellk{initialgrid}, ...
    cella{initialgrid},cellb{initialgrid});

% loop through vcycles
for p=1:num_vcycles_max

    % Stepping down
    for i=initialgrid:initialgrid+grids-2
    
        % Step down
        cellf{i+1}=Rmg(r,cellN{i+1});
        % clear v from previous loop
        cellv{i+1}=zeros(cellN{i+1},1);
        
        % If at coarsest level, solve exactly
        if i==initialgrid+grids-2
            cellv{i+1}=LA(eye(cellN{i+1},cellN{i+1}),cellk{i+1}, ...
                cella{i+1},cellb{i+1})\cellf{i+1};
        else
            % Iterate
            cellv{i+1}=pre_relaxation(cellv{i+1},cellk{i+1},cella{i+1}, ...
                cellb{i+1},cellf{i+1},cellH{i+1},Nd);
            
            % Find residual
            r=findR(cellf{i+1},cellv{i+1},cellk{i+1},cella{i+1}, ...
                cellb{i+1});
        end
        
    end

    % Stepping up
    for i=initialgrid+grids-1:-1:initialgrid+1
    
        % Step up and update guess
        cellv{i-1}=cellv{i-1}+Pmg(cellv{i},cellN{i-1});
        % Iterate
        cellv{i-1}=pre_relaxation(cellv{i-1},cellk{i-1},cella{i-1}, ...
            cellb{i-1},cellf{i-1},cellH{i-1},Nu);
        
    end

    % calculate v and r for next loop (or output)
    v=cellv{initialgrid};
    r=findR(cellf{initialgrid},cellv{initialgrid},cellk{initialgrid}, ...
        cella{initialgrid},cellb{initialgrid});
    
    if rms(r)<= tol
%         fprintf('Converged after %d v-cycle\n',p)
        break
    end
    
end

if p==num_vcycles_max
    
    fprintf('Did not converge to required tolerance after %d vcycles\n',p)
    fprintf('Residual = %d\n',rms(r))
    
end

end