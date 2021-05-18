% Function uses Newton vcycle iterations to solve equation 

% Preconditioned relaxations only!

% Inputs:
% pre_relaxation - iterative scheme used (see main)
% vcyclegrid - number of grids used in v-cycle
% cellN - cell of grid points
% cellk - cell of wave number
% num_vcycles - number of vcycles performed
% Nd - number of iterations per grid down cycle
% Nu - number of iterations per grid up cycle
% cella - cell of function a(x) in -u_xx+au_x+bu+cu^2=f
% cellb - cell of function b(x) in -u_xx+au_x+bu+cu^2=f
% cellc - cell of function c(x) in -u_xx+au_x+bu+cu^2=f
% cellRHS - cell of RHS of Au=f
% cellv - cell of v, best estimate of solution (or initial estimate)

% Ouputs:
% cellv - best guess after vcycle iterations (cellv{1})

function cellv=Newton_vcycle_pre(gamma,pre_relaxation,vcyclegrid,cellN,cellk, ...
    num_vcycles,Nd,Nu,cella,cellb,cellc,cellRHS,preconditioner,cellv)

global L

% Save base variable value
b=cellb{1};
RHS=cellRHS{1};

% Update variable wrt parameter
cellRHS{1}=gamma*RHS;

% New b(x) function in Newton
bnew=cellb{1}+2*cellc{1}.*cellv{1};
H=preconditioner(L,cellN{1},cella{1},bnew);

% Error guess (keep at 0)
e0=zeros(cellN{1},1);

% Initialise cells
cellbnew=cell(vcyclegrid,1);
cellf=cell(vcyclegrid,1);
celle=cell(vcyclegrid,1);

for i=1:20
    
    % Initial RHS of linear equation
    F=cellRHS{1}-NA(cellv{1},cellk{1},cella{1},cellb{1},cellc{1});
    
    % Set cell for bnew and step downs
    cellbnew=setcellsNewton(vcyclegrid,cellbnew,bnew);
    cellH=setcellsH(vcyclegrid,H,L,cellN,cella,cellbnew);
    
    % Update finest grid for f and e for iteration
    cellf{1}=F;
    celle{1}=e0;
    
    r=rms(F);
    fprintf('Residual Newton = %d\n',r)
    if r<=1e-10
        fprintf('Converged after %d Newton Iterations \n',i-1)
        break
    end
    
    % Solve linear equation
    e=vcycle_pre(pre_relaxation,1,vcyclegrid,cellN,cellk,num_vcycles,Nd,Nu, ...
        cella,cellbnew,cellf,cellH,celle);
    e=Pcg(celle{1},cellk{1},cella{1},cellbnew{1},cellf{1},cellH{1});
    % Update correction
    cellv{1}=cellv{1}+e; 
    
    % Update vector b for next iteration
    bnew=cellb{1}+2*cellc{1}.*cellv{1};
    H=preconditioner(L,cellN{1},cella{1},cellbnew{1});
    
end

if i==20
    
    fprintf('Did not converge to required tolerance after %d Newton Iterations\n',i)
    
end

end