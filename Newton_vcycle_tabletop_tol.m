% Function uses Newton vcycle iterations to solve equations
%
% Given a gamma value, function will solve for u and lambda satisfying
% F(u,lambda,gamma)=0 and G(u,lambda,gamma)=0
%
% Preconditioned relaxations only!
%
% Inputs:
% lambda - initial lambda guess in  -u_xx+au_x+b*lambda*u+cu^2=gamma*f
% gamma - fixed gamma value in      -u_xx+au_x+b*lambda*u+cu^2=gamma*f
% pre_relaxation - iterative scheme used (see main)
% vcyclegrid - number of grids used in v-cycle
% cellN - cell of grid points
% cellk - cell of wave number
% num_vcycles - number of vcycles performed
% Nd - number of iterations per grid down cycle
% Nu - number of iterations per grid up cycle
% cella - cell of function a(x) in -u_xx+au_x+b*lambda*u+cu^2=gamma*f
% cellb - cell of function b(x) in -u_xx+au_x+b*lambda*u+cu^2=gamma*f
% cellc - cell of function c(x) in -u_xx+au_x+b*lambda*u+cu^2=gamma*f
% cellRHS - cell of RHS of Au=f
% cellv - cell of v, best estimate of solution (or initial estimate)
%
% Ouputs:
% cellv - best solution guess after iterations (cellv{1})
% lambda - best lambda guess after iterations

function [cellv,lambda]=Newton_vcycle_tabletop_tol(lambda,gamma,x,pre_relaxation,vcyclegrid,cellN,cellk, ...
    tol,Nd,Nu,cella,cellb,cellc,cellRHS,preconditioner,cellv)

global L

% Save base variable value
b=cellb{1};
RHS=cellRHS{1};

% Set cells for z1,z2,RHS1,RHS2 (for multigrid solver)
cellz1=cell(vcyclegrid,1);
cellz2=cell(vcyclegrid,1);
cellRHS1=cell(vcyclegrid,1);
cellRHS2=cell(vcyclegrid,1);
cellbnew=cell(vcyclegrid,1);

cellz1{1}=zeros(cellN{1},1);
cellz2{1}=cellz1{1};

% Update variable wrt parameter
cellb{1}=lambda*b;
cellRHS{1}=gamma*RHS;

% New b(x) function in Newton
bnew=cellb{1}+2*cellc{1}.*cellv{1};
H=preconditioner(L,cellN{1},cella{1},bnew);

% G_u calculation (no loop dependence)
RHSx=real(ifft(1i*cellk{1}'.*fft(RHS)));
cn=zeros(cellN{1},1);
cn(2:cellN{1}/2)=L/cellN{1};
cn(1)=1/2*L/cellN{1};
cn(cellN{1}/2+1)=1/2*L/cellN{1};
G_u=gamma*RHSx.*cn;
    
for i=1:20
    
    % F and F derivatives
    F=NA(cellv{1},cellk{1},cella{1},cellb{1},cellc{1})-cellRHS{1};
    F_lambda=sign(b).*cellv{1};
    
    % G calculation
    G=1/54*lambda^3+trapint(gamma*RHSx(1:cellN{1}/2+1).*cellv{1}(1:cellN{1}/2+1),x(1:cellN{1}/2+1));

    % G_lambda
    G_lambda=1/18*lambda^2;
            
    % Check convergence condition
    if rms(F)<=1e-8 && rms(G)<=1e-8
        fprintf('Converged after %d Newton Iterations\n',i-1)
    	break % End if converged to sufficient accuracy
    end
    disp(rms(F));disp(rms(G));
    % Initial RHS of linear equation
    RHS1=-F;
    RHS2=F_lambda;
       
    % Set cell for bnew and step downs
    cellbnew=setcellsNewton(vcyclegrid,cellbnew,bnew);
    cellH=setcellsH(vcyclegrid,H,L,cellN,cella,cellbnew);
    
    % Update cell RHS with new values
    [cellRHS1,cellRHS2]=setcellspseudo(vcyclegrid,cellN,cellRHS1, ...
        cellRHS2,RHS1,RHS2); 
    
    % Solve F_u*z1=-F and F_u*z2=F_lambda (F_u is jacobian)
%     [cellz1{1},r1]=vcycle_pre_tol(pre_relaxation,1,vcyclegrid,cellN,cellk,tol,Nd,Nu,cella,cellbnew,cellRHS1,cellH,cellz1);
%     [cellz2{1},r2]=vcycle_pre_tol(pre_relaxation,1,vcyclegrid,cellN,cellk,tol,Nd,Nu,cella,cellbnew,cellRHS2,cellH,cellz2);
    [cellz1{1},r1]=Pcg(cellz1{1},cellk{1},cella{1},cellbnew{1},cellRHS1{1},H);
    [cellz2{1},r2]=Pcg(cellz2{1},cellk{1},cella{1},cellbnew{1},cellRHS2{1},H);
%     cellz1{1}=LA(eye(cellN{1},cellN{1}),cellk{1},cella{1},cellbnew{1})\cellRHS1{1};
%     cellz2{1}=LA(eye(cellN{1},cellN{1}),cellk{1},cella{1},cellbnew{1})\cellRHS2{1};    
    
    % Solving for delta_v and delta_lambda
    delta_lambda=(-G-dot(G_u,cellz1{1}))/(G_lambda-dot(G_u,cellz2{1}));
    delta_v=cellz1{1}-cellz2{1}*delta_lambda;
    
    % Update variable for next Newton iteration
    cellv{1}=cellv{1}+delta_v;
    lambda=lambda+delta_lambda; 
    
    % Update variable wrt parameter
    cellb{1}=lambda*b;    
    
    % Update vector b for next iteration
    bnew=cellb{1}+2*cellc{1}.*cellv{1};
    H=preconditioner(L,cellN{1},cella{1},cellbnew{1}); 
    
end

if i==20
    
    fprintf('Did not converge to required tolerance after %d Newton Iterations\n',i)
    
end

end