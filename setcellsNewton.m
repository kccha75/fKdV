% Function to set cells for arrays of different lengths
% cell{1} - finest grid
% cell{2} - next finest grid etc etc

% Inputs:
% vcyclegrids - grids to vcycle through
% cellbnew - cell of function new b(x) in -u_xx+au_x+bu=f for Newton
% bnew - new b(x) function for Newton method (finest grid)

% Outputs:
% cellbnew - cell updated with step down variants

function cellbnew=setcellsNewton(vcyclegrid,cellbnew,bnew)

% Set fine grid point parameters to cell 1
cellbnew{1}=bnew;

% Loop to set parameters for coarse grids
for i=2:vcyclegrid
   
    % step down
    cellbnew{i}=cellbnew{i-1}(1:2:end);
%     cellbnew{i}=Rmg(cellbnew{i-1},length(cellbnew{i-1})/2);
    % step down RHS (for FMG only)
    % cellf{i}=Rmg(cellf{i-1},cellN{i});
  
end

end