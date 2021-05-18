% Function to set cells for arrays of different lengths
% cell{1} - finest grid
% cell{2} - next finest grid etc etc

% Inputs:
% vcyclegrids - grids to vcycle through
% H - Preconditioning matrix in finest grid
% L - domain length
% cellN - cell of grid points
% cella - cell of function a(x) in -u_xx+au_x+bu+cu^2=f
% cellb - cell of function b(x) in -u_xx+au_x+bu+cu^2=f

% Outputs:
% cellH - cell of preconditioner updated with step down variants
function cellH=setcellsH(vcyclegrid,H,L,cellN,cella,cellb)

cellH=cell(vcyclegrid,1);
cellH{1}=H;

for i=2:vcyclegrid
    
    cellH{i}=Hfd(L,cellN{i},cella{i},cellb{i});
    
end