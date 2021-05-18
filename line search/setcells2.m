% Function to set cells for arrays of different lengths
% cell{1} - finest grid
% cell{2} - next finest grid etc etc

% Inputs:
% vcyclegrids - grids to vcycle through
% N - grid points (finest grid)
% k - wave number (finest grid)
% a - a(x) function (finest grid)
% b - b(x) function (finest grid)
% f - RHS of Au=f (finest grid)
% v0 - initial guess (finest grid)

% Outputs:
% cellN - cell of grid points
% cellk - cell of wave number
% cella - cell of function a(x) in -u_xx+au_x+bu+cu^2=f
% cellb - cell of function b(x) in -u_xx+au_x+bu+cu^2=f
% cellc - cell of function b(x) in -u_xx+au_x+bu+cu^2=f
% cellRHS - cell of RHS of Au=f
% cellv - cell of best estimate of solution (or initial estimate)

function [cellN,cellk,cella,cellb,cellc,cellRHS,cellv]=setcells2(vcyclegrid,N,k,a,b,c,RHS,v0)

% Define cells to store arrays of different lengths
cellN=cell(vcyclegrid,1);
cellk=cell(vcyclegrid,1);
cellRHS=cell(vcyclegrid,1);
cella=cell(vcyclegrid,1);
cellb=cell(vcyclegrid,1);
cellc=cell(vcyclegrid,1);
cellv=cell(vcyclegrid,1);

% Set fine grid point parameters to cell 1
cellN{1}=N;
cellk{1}=k;
cellRHS{1}=RHS;
cella{1}=a;
cellb{1}=b;
cellc{1}=c;
cellv{1}=v0;

% Loop to set parameters for coarse grids
for i=2:vcyclegrid
   
    % step down
    cellN{i}=cellN{i-1}/2;
    cellk{i}=[cellk{i-1}(1:cellN{i-1}/4) 0 cellk{i-1}(3*cellN{i-1}/4+2:end)];
    cella{i}=cella{i-1}(1:2:end);
    cellb{i}=cellb{i-1}(1:2:end);
    cellc{i}=cellc{i-1}(1:2:end);
    
    % step down RHS (for FMG only)
    % cellRHS{i}=Rmg(cellRHS{i-1},cellN{i});
  
end

end