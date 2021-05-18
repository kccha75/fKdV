% Function performs Full Multigrid from coarest grid to finest grid using
% vcycle function

% FOR NON-PRECONDITIONING ONLY

% Inputs:
% relaxation - function input for choice of relaxation
% vcyclegrid - number of grids to cycle through
% cellN - cell of grid points
% cellk - cell of wave number
% num_vcycles - number of vcycles performed 
% Nd - number of iterations per grid down cycle
% Nu - number of iterations per grid up cycle
% cella - cell of function a(x) in -u_xx+au_x+bu=f
% cellb - cell of function b(x) in -u_xx+au_x+bu=f
% cellf - cell of RHS of Au=f
% cellv - cell of best estimate of solution (or initial estimate)

% Ouputs:
% cellv - best guess after vcycle iterations
% r - residual of best guess

function [cellv,r]=FMG(relaxation,vcyclegrid,cellN,cellk,num_vcycles,Nd,Nu,cella,cellb,cellf,cellv)

% Initial coarse grid guess
cellv{vcyclegrid}=zeros(cellN{vcyclegrid},1);

% start FMG
j=2;
for i=vcyclegrid:-1:2
    
    % step up
    cellv{i-1}=Pmg(cellv{i},cellN{i-1});
    % vcycle
    [cellv{i-1},r]=vcycle(relaxation,i-1,j,cellN,cellk,num_vcycles,Nd,Nu,cella,cellb,cellf,cellv);
    j=j+1;
    
end

end