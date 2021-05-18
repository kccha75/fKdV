% Function uses Pseudo arclength continuation to trace out solution family

% Non-preconditioned relaxations only!

% Uses fixed step length ds
% Continues until steps or if Newton does not converge (after 20
% iterations)

% Algorithm follows that of solution technique in nonlinear oceanography

% Inputs:
% v - initial solution 
% lambda - initial lambda
% dv - initial tangent of v (or approximation of)
% dlambda - initial tangent of lambda (or approximation of)
% ds - step size (constant)
% steps - number of steps to take
% x - vector of x values
% relaxation - iterative scheme used (see main)
% vcyclegrid - number of grids used in v-cycle
% cellN - cell of grid points
% cellk - cell of wave number
% Nd - number of iterations per grid down cycle
% Nu - number of iterations per grid up cycle
% cella - cell of function a(x) in -u_xx+au_x+bu+cu^2=f
% cellb - cell of function b(x) in -u_xx+au_x+bu+cu^2=f
% cellc - cell of function c(x) in -u_xx+au_x+bu+cu^2=f

% Outputs:
% v - vector of solution curves
% lambda - vector of lambda values

function [v,lambda]=pseudoarclength(v,lambda,dv,dlambda,ds,steps,... 
    x,relaxation,vcyclegrid,cellN,cellk,Nd,Nu,cella,cellb,cellc,cellRHS)
     
% Save base variable value
b=cellb{1};

% Set cells for z1,z2,RHS1,RHS2 (for FMG solver)
cellz1=cell(vcyclegrid,1);
cellz2=cell(vcyclegrid,1);
cellRHS1=cell(vcyclegrid,1);
cellRHS2=cell(vcyclegrid,1);
cellbnew=cell(vcyclegrid,1);

cellz1{1}=zeros(cellN{1},1);
cellz2{1}=cellz1{1};

for j=1:steps
    
    % Tangent predictor
    v(:,j+1)=v(:,j)+ds*dv;
    lambda(j+1)=lambda(j)+ds*dlambda;
    
    % Update variable wrt parameter
    cellb{1}=lambda(j+1)*b;
    
    % New b(x) function in Newton
    bnew=cellb{1}+2*cellc{1}.*v(:,j+1);

    % Correction loop using Newton iterations
    for i=1:20
        
        % Initial RHS of linear equation
        RHS1=-(NA(v(:,j+1),cellk{1},cella{1},cellb{1},cellc{1})- ...
            cellRHS{1});
        RHS2=findF_lambda(cella{1},cellb{1},cellc{1},v(:,j+1),x,lambda(j+1)); 
        
        % Set cell for bnew and step downs
        cellbnew=setcellsNewton(vcyclegrid,cellbnew,bnew);

        % Update cell RHS with new values
        [cellRHS1,cellRHS2]=setcellspseudo(vcyclegrid,cellN,cellRHS1, ...
            cellRHS2,RHS1,RHS2); 
        
        % FMG to solve F_x*z1=-F
        [cellz1,r1]=FMG(relaxation,vcyclegrid,cellN,cellk,1,Nd,Nu,cella,cellbnew,cellRHS1,cellz1);
%         [cellz1{1},r1]=vcycle(relaxation,1,vcyclegrid,cellN,cellk,10,Nd,Nu,cella,cellbnew,cellRHS1,cellz1);
%         cellz1{1}=LA(eye(cellN{1},cellN{1}),cellk{1},cella{1},cellbnew{1})\cellRHS1{1}; % Exact solver
%         [cellz1{1},r1]=Pcg(zeros(cellN{1},1),cellk{1},cella{1},cellbnew{1},cellRHS1{1}); % Pcg
%         r1=cellRHS1{1}-LA(cellz1{1},cellk{1},cella{1},cellbnew{1});
%          disp(rms(r1))

        % FMG to solve F_x*z2=F_lambda
        [cellz2,r2]=FMG(relaxation,vcyclegrid,cellN,cellk,1,Nd,Nu,cella,cellbnew,cellRHS2,cellz2);
%         [cellz2{1},r2]=vcycle(relaxation,1,vcyclegrid,cellN,cellk,10,Nd,Nu,cella,cellbnew,cellRHS2,cellz2);
%         cellz2{1}=LA(eye(cellN{1},cellN{1}),cellk{1},cella{1},cellbnew{1})\cellRHS2{1}; % Exact solver
%         [cellz2{1},r2]=Pcg(zeros(cellN{1},1),cellk{1},cella{1},cellbnew{1},cellRHS2{1}); % Pcg
%         r2=cellRHS2{1}-LA(cellz2{1},cellk{1},cella{1},cellbnew{1});
%          disp(rms(r2))

        % Solving for delta_v and delta_lambda (see nonlinear oceanography)
        delta_lambda=(ds-dot(dv,(v(:,j+1)-v(:,j)))-dlambda*(lambda(j+1)- ...
            lambda(j))-dot(dv,cellz1{1}))/(dlambda-dot(dv,cellz2{1}));
        delta_v=cellz1{1}-delta_lambda*cellz2{1};    
        
        % if sqrt(rms(dot([delta_v;delta_lambda],[delta_v;delta_lambda])))<=1e-10
%         if abs(rms(delta_v))<=1e-10 && abs(delta_lambda)<=1e-10 %
        if rms(NA(v(:,j+1),cellk{1},cella{1},cellb{1},cellc{1})-cellRHS{1})<=1e-10
            fprintf('Converged after %d Newton Iterations step = %d\n',i,j)
            break % End if converged to sufficient accuracy
        end
%         disp(abs(rms(delta_v)));disp(abs(delta_lambda));
        % Update correction
        v(:,j+1)=v(:,j+1)+delta_v;
        lambda(j+1)=lambda(j+1)+delta_lambda;
        
        % Update variable wrt parameter
        cellb{1}=lambda(j+1)*b;
        
        % Update vector b for next iteration
        bnew=cellb{1}+2*cellc{1}.*v(:,j+1);
        
    end

    if i==20
        fprintf('Did not converge to required tolerance after %d Newton Iterations at step %d\n',i,j)
        % Do not save latest vectors if not converged
        v(:,j+1)=[];
        lambda(:,j+1)=[];
        return % end function
    end

    % Solve for new direction
%     RHS1=zeros(cellN{1},1);
%     RHS2=findF_lambda(x,lambda(j+1));
    
    % Update cell RHS with new values
%     [cellRHS1,cellRHS2]=setcellspseudo(vcyclegrid,cellN,cellRHS1, ...
%         cellRHS2,RHS1,RHS2); 
        
    % FMG to solve F_x*z1=0 SOLUTION IS JUST 0!
%     [cellz1,r1]=FMG(relaxation,vcyclegrid,cellN,cellk,10,Nd,Nu,cella,cellbnew,cellRHS1,cellz1);
%     cellz1{1}=LA(eye(cellN{1},cellN{1}),cellk{1},cella{1},cellbnew{1})\RHS1; % Exact solver
%     [cellz1{1},r1]=Pcg(zeros(cellN{1},1),cellk{1},cella{1},cellbnew{1},RHS1); % Pcg
%     r1=RHS1-LA(cellz1{1},cellk{1},cella{1},cellbnew{1});
%     disp(rms(r1))

    % FMG to solve F_x*z2=F_lambda
%     [cellz2,r2]=FMG(relaxation,vcyclegrid,cellN,cellk,10,Nd,Nu,cella,cellbnew,cellRHS2,cellz2);
%     cellz2{1}=LA(eye(cellN{1},cellN{1}),cellk{1},cella{1},cellbnew{1})\RHS2; % Exact solver
%     [cellz2{1},r2]=Pcg(zeros(cellN{1},1),cellk{1},cella{1},cellbnew{1},RHS2); % Pcg
%     r2=RHS2-LA(cellz2{1},cellk{1},cella{1},cellbnew{1});
%     disp(rms(r2))

    % Solving for dlambda and dv
    dlambda=1/(dlambda-dot(dv,cellz2{1}));
    dv=-dlambda*cellz2{1};
    
    % Normalise
    mag=sqrt(dot(dv,dv)+dlambda^2);
    dlambda=dlambda/mag;
    dv=dv/mag;
    
end