% Function to set cells for arrays of different lengths
% cell{1} - finest grid
% cell{2} - next finest grid etc etc
%
% Inputs:
% vcyclegrids - grids to vcycle through
% cellN - cell of grid points
% cellRHS1 - cell of RHS1 old
% cellRHS2 - cell of RHS2 old
% cellRHS3 - cell of RHS3 old
% RHS1 - RHS1 function
% RHS2 - RHS2 function
% RHS3 - RHS3 function
%
% Outputs:
% cellRHS1 - cell of RHS1 updated
% cellRHS2 - cell of RHS2 updated
% cellRHS3 - cell of RHS3 updated
% 
% Application to Pseudo Arclength continuation solving
% (see Nonlinear oceanography for solution technique)

function [cellRHS1,cellRHS2,cellRHS3]=setcellspseudo2(vcyclegrid,cellN,cellRHS1,cellRHS2,cellRHS3,RHS1,RHS2,RHS3)

% Set fine grid point parameters to cell 1
cellRHS1{1}=RHS1;
cellRHS2{1}=RHS2;
cellRHS3{1}=RHS3;

% Loop to set parameters for coarse grids
for i=2:vcyclegrid
   
    % step down RHS (for FMG only)
    cellRHS1{i}=Rmg(cellRHS1{i-1},cellN{i});
    cellRHS2{i}=Rmg(cellRHS2{i-1},cellN{i});
    cellRHS3{i}=Rmg(cellRHS3{i-1},cellN{i});
  
end

end